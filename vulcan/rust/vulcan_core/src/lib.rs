//! vulcan_core: Vulcan L1 proteome-index primitives (DESIGN.md §4 L1, §5).
//!
//! Reuses sage-core (MIT, Michael Lazear; https://github.com/lazear/sage, commit 99407db6 = v0.14.7):
//!   * FASTA parsing           -> sage_core::fasta::Fasta
//!   * enzymatic digestion     -> sage_core::enzyme::{Enzyme, EnzymeParameters}
//!   * modification expansion  -> sage_core::peptide::Peptide::apply (variable then static mods)
//!   * reference ion series    -> sage_core::ion_series::IonSeries (tests only, `sage_reference`)
//!
//! Vulcan-specific (not in Sage): I/L-collapsed peptide keys, peptidoform pairs whose decoy side is a
//! *given* permutation of the target (pseudo-reverse with the C-term fixed, or a hash-seeded shuffle;
//! the permutation is chosen in Python, `vulcan/index/decoys.py`), f64 masses, var-mod position masks
//! and hash-assigned u64 ids (xxhash64 of the canonical peptidoform string; no storage order and no
//! label enters an id).
use std::collections::HashMap;
use std::hash::Hasher;
use std::str::FromStr;
use std::sync::Arc;

use numpy::{IntoPyArray, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rayon::prelude::*;

use sage_core::enzyme::{Digest, Enzyme, EnzymeParameters, Position};
use sage_core::fasta::Fasta;
use sage_core::ion_series::{IonSeries, Kind};
use sage_core::mass::{monoisotopic as sage_mono, H2O as SAGE_H2O};
use sage_core::modification::ModificationSpecificity;
use sage_core::peptide::Peptide;

/// Monoisotopic residue masses (Unimod, f64). Index = ASCII - b'A'. 0.0 = not a standard residue.
pub const RESIDUE_F64: [f64; 26] = [
    71.037113805,  // A
    0.0,           // B
    103.009184505, // C
    115.026943065, // D
    129.042593135, // E
    147.068413945, // F
    57.021463721,  // G
    137.058911875, // H
    113.084064015, // I
    0.0,           // J
    128.094963050, // K
    113.084064015, // L
    131.040484645, // M
    114.042927470, // N
    0.0,           // O (pyrrolysine, excluded)
    97.052763875,  // P
    128.058577540, // Q
    156.101111050, // R
    87.032028435,  // S
    101.047678505, // T
    0.0,           // U (selenocysteine, excluded)
    99.068413945,  // V
    186.079312980, // W
    0.0,           // X
    163.063328575, // Y
    0.0,           // Z
];
pub const H2O_F64: f64 = 18.010564684;

fn is_standard(s: &str) -> bool {
    s.bytes().all(|c| c.is_ascii_uppercase() && RESIDUE_F64[(c - b'A') as usize] > 0.0)
}

fn il(s: &str) -> String {
    s.replace('I', "L")
}

fn parse_mods(mods: &[(String, f64)]) -> PyResult<Vec<(ModificationSpecificity, f32)>> {
    mods.iter()
        .map(|(k, v)| {
            ModificationSpecificity::from_str(k)
                .map(|m| (m, *v as f32))
                .map_err(|e| PyValueError::new_err(format!("bad modification spec {k}: {e:?}")))
        })
        .collect()
}

/// Digest a FASTA (text) with Sage's enzyme code.
///
/// Returns a dict with
///   keys: list[str]  unique I/L-collapsed peptide sequences (standard residues only), sorted
///   missed: u8[n]    missed cleavages (Sage's count)
///   prot_off: i64[n+1], prot_idx: u32[...]  protein indices of each key (CSR)
///   accessions: list[str] protein accessions in FASTA order (first whitespace token of the header)
///   n_nonstandard: number of distinct digests dropped for non-standard residues
#[pyfunction]
#[pyo3(signature = (fasta_text, missed_cleavages, min_len, max_len, cleave_at, restrict, c_terminal))]
fn digest_fasta<'py>(
    py: Python<'py>,
    fasta_text: String,
    missed_cleavages: u8,
    min_len: usize,
    max_len: usize,
    cleave_at: &str,
    restrict: Option<char>,
    c_terminal: bool,
) -> PyResult<Bound<'py, PyDict>> {
    // generate_decoys=true + tag "rev_": entries whose accession contains "rev_" are skipped; Vulcan
    // makes its own decoys at peptide level.
    let fasta = Fasta::parse(fasta_text, "rev_", true);
    let accessions: Vec<String> = fasta.targets.iter().map(|(a, _)| a.as_ref().clone()).collect();
    let acc_idx: HashMap<String, u32> =
        accessions.iter().enumerate().map(|(i, a)| (a.clone(), i as u32)).collect();
    let enzyme = EnzymeParameters {
        missed_cleavages,
        min_len,
        max_len,
        enyzme: Enzyme::new(cleave_at, restrict, c_terminal, false),
    };
    let digests: Vec<Digest> = py.allow_threads(|| fasta.digest(&enzyme));

    let mut map: HashMap<String, (u8, Vec<u32>)> = HashMap::with_capacity(digests.len());
    let mut nonstd: std::collections::HashSet<String> = Default::default();
    for d in digests.iter() {
        if !is_standard(&d.sequence) {
            nonstd.insert(d.sequence.clone());
            continue;
        }
        let k = il(&d.sequence);
        let p = acc_idx[d.protein.as_ref()];
        let e = map.entry(k).or_insert((d.missed_cleavages, Vec::new()));
        e.0 = e.0.min(d.missed_cleavages);
        e.1.push(p);
    }
    let mut keys: Vec<String> = map.keys().cloned().collect();
    keys.sort_unstable();
    let mut missed = Vec::with_capacity(keys.len());
    let mut prot_off: Vec<i64> = Vec::with_capacity(keys.len() + 1);
    let mut prot_idx: Vec<u32> = Vec::new();
    prot_off.push(0);
    for k in keys.iter() {
        let (m, ps) = map.get_mut(k).unwrap();
        ps.sort_unstable();
        ps.dedup();
        missed.push(*m);
        prot_idx.extend_from_slice(ps);
        prot_off.push(prot_idx.len() as i64);
    }
    let out = PyDict::new_bound(py);
    out.set_item("keys", keys)?;
    out.set_item("missed", missed.into_pyarray_bound(py))?;
    out.set_item("prot_off", prot_off.into_pyarray_bound(py))?;
    out.set_item("prot_idx", prot_idx.into_pyarray_bound(py))?;
    out.set_item("accessions", accessions)?;
    out.set_item("n_nonstandard", nonstd.len())?;
    Ok(out)
}

struct ModCfg {
    statics: HashMap<ModificationSpecificity, f32>,
    vars: Vec<(ModificationSpecificity, f32)>,
    /// per residue (A..Z): static delta (f64), var delta (f64, 0 = none)
    static_res: [f64; 26],
    var_res: [f64; 26],
    /// variable peptide-N-term mod delta (Sage `^`), 0 = none.  Carried in mask bit 63 (NTERM_BIT).
    nterm_var: f64,
    max_var: usize,
}

/// Mask bit for a variable peptide-N-terminal modification (M4).  Residue bits are 0..len-1 (len <= 50),
/// so bit 63 never collides.  The decoy keeps the bit (its own N-terminus), so pairs stay mass-identical.
pub const NTERM_BIT: u64 = 1u64 << 63;

fn mod_cfg(static_mods: &[(String, f64)], var_mods: &[(String, f64)], max_var: usize) -> PyResult<ModCfg> {
    let st = parse_mods(static_mods)?;
    let vr = parse_mods(var_mods)?;
    let mut static_res = [0.0f64; 26];
    let mut var_res = [0.0f64; 26];
    let mut nterm_var = 0.0f64;
    for ((spec, _), (_, m)) in st.iter().zip(static_mods.iter()) {
        if let ModificationSpecificity::Residue(r) = spec {
            static_res[(r - b'A') as usize] = *m;
        } else {
            return Err(PyValueError::new_err("only residue static mods are supported in the M2 MVP"));
        }
    }
    for ((spec, _), (_, m)) in vr.iter().zip(var_mods.iter()) {
        if let ModificationSpecificity::Residue(r) = spec {
            let i = (r - b'A') as usize;
            if var_res[i] != 0.0 || static_res[i] != 0.0 {
                return Err(PyValueError::new_err(
                    "M2 MVP: at most one variable mod per residue, and not on a statically modified residue",
                ));
            }
            var_res[i] = *m;
        } else if let ModificationSpecificity::PeptideN(None) = spec {
            if nterm_var != 0.0 {
                return Err(PyValueError::new_err("at most one variable peptide-N-term mod"));
            }
            nterm_var = *m;
        } else {
            return Err(PyValueError::new_err(
                "only residue and peptide-N-term (^) variable mods are supported (protein-N-term needs digest positions)",
            ));
        }
    }
    Ok(ModCfg { statics: st.into_iter().collect(), vars: vr, static_res, var_res, nterm_var, max_var })
}

fn sage_peptide(seq: &str) -> Peptide {
    let d = Digest {
        decoy: false,
        semi_enzymatic: false,
        sequence: seq.to_string(),
        protein: Arc::new(String::new()),
        missed_cleavages: 0,
        position: Position::Internal,
    };
    Peptide::try_from(d).expect("standard residues only")
}

/// Canonical peptidoform string: residues with "[+delta]" after variable-modified residues.
fn canonical(seq: &[u8], mask: u64, cfg: &ModCfg) -> String {
    let mut s = String::with_capacity(seq.len() + 12);
    if mask & NTERM_BIT != 0 {
        s.push_str(&format!("[{:+.4}]-", cfg.nterm_var));
    }
    for (k, &c) in seq.iter().enumerate() {
        s.push(c as char);
        if mask >> k & 1 == 1 {
            s.push_str(&format!("[{:+.4}]", cfg.var_res[(c - b'A') as usize]));
        }
    }
    s
}

fn xxh64(s: &str, seed: u64) -> u64 {
    let mut h = twox_hash::XxHash64::with_seed(seed);
    h.write(s.as_bytes());
    h.finish()
}

/// Sage's own f32 monoisotopic mass of a peptidoform (peptide.rs: H2O + residues in sequence order,
/// then + the f32 sum of per-residue modification masses), so precursor windows match Sage bit-for-bit.
fn sage_mass32(seq: &[u8], mask: u64, cfg: &ModCfg) -> f32 {
    let mut m: f32 = SAGE_H2O;
    for &c in seq.iter() {
        m += sage_mono(c);
    }
    let mut ms: f32 = 0.0;
    for (k, &c) in seq.iter().enumerate() {
        let i = (c - b'A') as usize;
        let d = if mask >> k & 1 == 1 { cfg.var_res[i] as f32 } else { cfg.static_res[i] as f32 };
        ms += d;
    }
    if mask & NTERM_BIT != 0 {
        // Sage: monoisotopic += (sum of residue mods) + nterm   (peptide.rs modification_mass)
        ms += cfg.nterm_var as f32;
    }
    m + ms
}

fn form_mass(seq: &[u8], mask: u64, cfg: &ModCfg) -> f64 {
    let mut m = H2O_F64;
    for (k, &c) in seq.iter().enumerate() {
        let i = (c - b'A') as usize;
        m += RESIDUE_F64[i] + cfg.static_res[i];
        if mask >> k & 1 == 1 {
            m += cfg.var_res[i];
        }
    }
    if mask & NTERM_BIT != 0 {
        m += cfg.nterm_var;
    }
    m
}

/// Expand target/decoy *pairs* of sequences into paired peptidoforms.
///
/// fwd[i] and dec[i] are a forward (target, entrapment or foreign) sequence and its decoy, with
/// dec[i][j] == fwd[i][perm_i[j]]. perm_i is perm_flat[perm_off[i]..perm_off[i+1]]; an empty slice
/// means the pseudo-reverse with the C-term fixed (perm[j] = n-2-j for j < n-1, perm[n-1] = n-1).
/// Variable-mod configurations are enumerated on the forward side with Sage's `Peptide::apply` and
/// carried to the decoy through the permutation, so each forward peptidoform has exactly one decoy
/// peptidoform of identical mass (the pair). Forms outside [min_mass, max_mass] are dropped (both sides).
///
/// Returns dict of numpy arrays, one row per pair of peptidoforms:
///   seq (u32, index i), mask_f, mask_d (u64), mass (f64), n_var (u8), id_f, id_d (u64),
///   mass32_f, mass32_d (f32: Sage's own f32 mass of each side, sequence-order summation)
#[pyfunction]
#[pyo3(signature = (fwd, dec, perm_flat, perm_off, static_mods, var_mods, max_var, min_mass, max_mass, id_seed))]
#[allow(clippy::too_many_arguments)]
fn expand_pairs<'py>(
    py: Python<'py>,
    fwd: Vec<String>,
    dec: Vec<String>,
    perm_flat: PyReadonlyArray1<'py, u8>,
    perm_off: PyReadonlyArray1<'py, i64>,
    static_mods: Vec<(String, f64)>,
    var_mods: Vec<(String, f64)>,
    max_var: usize,
    min_mass: f64,
    max_mass: f64,
    id_seed: u64,
) -> PyResult<Bound<'py, PyDict>> {
    if fwd.len() != dec.len() {
        return Err(PyValueError::new_err("fwd and dec differ in length"));
    }
    let cfg = mod_cfg(&static_mods, &var_mods, max_var)?;
    let perm_flat = perm_flat.as_slice()?.to_vec();
    let perm_off = perm_off.as_slice()?.to_vec();
    if perm_off.len() != fwd.len() + 1 {
        return Err(PyValueError::new_err("perm_off must have len(fwd)+1 entries"));
    }
    for (f, d) in fwd.iter().zip(dec.iter()) {
        if f.len() != d.len() || f.len() > 64 || !is_standard(f) || !is_standard(d) {
            return Err(PyValueError::new_err(format!("bad pair {f} / {d}")));
        }
    }
    type Row = (u32, u64, u64, f64, u8, u64, u64, f32, f32);
    let rows: Result<Vec<Vec<Row>>, String> = py.allow_threads(|| {
        fwd.par_iter()
            .zip(dec.par_iter())
            .enumerate()
            .map(|(i, (f, d))| {
                let fb = f.as_bytes();
                let db = d.as_bytes();
                let n = fb.len();
                let (a, b) = (perm_off[i] as usize, perm_off[i + 1] as usize);
                let perm: Vec<usize> = if a == b {
                    (0..n).map(|j| if j + 1 < n { n - 2 - j } else { n - 1 }).collect()
                } else {
                    perm_flat[a..b].iter().map(|&x| x as usize).collect()
                };
                if perm.len() != n || !(0..n).all(|j| perm[j] < n && db[j] == fb[perm[j]]) {
                    return Err(format!("pair {i}: decoy is not the given permutation of the forward sequence"));
                }
                let forms = sage_peptide(f).apply(&cfg.vars, &cfg.statics, cfg.max_var);
                let mut out = Vec::with_capacity(forms.len());
                for p in forms.iter() {
                    // variable-mod mask: residues carrying their (unique) variable mod delta
                    let mut mask = 0u64;
                    for (k, (&c, &m)) in p.sequence.iter().zip(p.modifications.iter()).enumerate() {
                        let v = cfg.var_res[(c - b'A') as usize];
                        if v != 0.0 && (m as f64 - v).abs() < 1e-3 {
                            mask |= 1 << k;
                        }
                    }
                    if p.nterm.is_some() {
                        mask |= NTERM_BIT;
                    }
                    let mass = form_mass(fb, mask, &cfg);
                    if mass < min_mass || mass > max_mass {
                        continue;
                    }
                    let mut mask_d = mask & NTERM_BIT;
                    for j in 0..n {
                        if mask >> perm[j] & 1 == 1 {
                            mask_d |= 1 << j;
                        }
                    }
                    let id_f = xxh64(&canonical(fb, mask, &cfg), id_seed);
                    let id_d = xxh64(&canonical(db, mask_d, &cfg), id_seed);
                    let m32f = sage_mass32(fb, mask, &cfg);
                    let m32d = sage_mass32(db, mask_d, &cfg);
                    out.push((i as u32, mask, mask_d, mass, mask.count_ones() as u8, id_f, id_d, m32f, m32d));
                }
                Ok(out)
            })
            .collect()
    });
    let rows = rows.map_err(PyValueError::new_err)?;
    let n: usize = rows.iter().map(|r| r.len()).sum();
    let (mut m32f, mut m32d): (Vec<f32>, Vec<f32>) = (Vec::with_capacity(n), Vec::with_capacity(n));
    let (mut seq, mut mf, mut md, mut mass, mut nv, mut idf, mut idd) = (
        Vec::with_capacity(n),
        Vec::with_capacity(n),
        Vec::with_capacity(n),
        Vec::with_capacity(n),
        Vec::with_capacity(n),
        Vec::with_capacity(n),
        Vec::with_capacity(n),
    );
    for r in rows.into_iter().flatten() {
        seq.push(r.0);
        mf.push(r.1);
        md.push(r.2);
        mass.push(r.3);
        nv.push(r.4);
        idf.push(r.5);
        idd.push(r.6);
        m32f.push(r.7);
        m32d.push(r.8);
    }
    let out = PyDict::new_bound(py);
    out.set_item("seq", seq.into_pyarray_bound(py))?;
    out.set_item("mask_f", mf.into_pyarray_bound(py))?;
    out.set_item("mask_d", md.into_pyarray_bound(py))?;
    out.set_item("mass", mass.into_pyarray_bound(py))?;
    out.set_item("n_var", nv.into_pyarray_bound(py))?;
    out.set_item("id_f", idf.into_pyarray_bound(py))?;
    out.set_item("id_d", idd.into_pyarray_bound(py))?;
    out.set_item("mass32_f", m32f.into_pyarray_bound(py))?;
    out.set_item("mass32_d", m32d.into_pyarray_bound(py))?;
    Ok(out)
}

/// Canonical id of one peptidoform (same hash as `expand_pairs`), for tests and reporting.
#[pyfunction]
#[pyo3(signature = (seq, mask, var_mods, id_seed))]
fn peptidoform_id(seq: &str, mask: u64, var_mods: Vec<(String, f64)>, id_seed: u64) -> PyResult<(String, u64)> {
    let cfg = mod_cfg(&[], &var_mods, 0)?;
    let s = canonical(seq.as_bytes(), mask, &cfg);
    let h = xxh64(&s, id_seed);
    Ok((s, h))
}

/// Sage reference for tests: for each sequence, every peptidoform from Sage's `Peptide::apply`
/// with Sage's own f32 monoisotopic mass and Sage's b/y neutral fragment masses (IonSeries).
/// Returns list of (seq index, sage display string, var mask, sage mass, b list, y list).
#[pyfunction]
#[pyo3(signature = (seqs, static_mods, var_mods, max_var))]
#[allow(clippy::type_complexity)]
fn sage_reference(
    seqs: Vec<String>,
    static_mods: Vec<(String, f64)>,
    var_mods: Vec<(String, f64)>,
    max_var: usize,
) -> PyResult<Vec<(usize, String, u64, f32, Vec<f32>, Vec<f32>)>> {
    let cfg = mod_cfg(&static_mods, &var_mods, max_var)?;
    let mut out = Vec::new();
    for (i, s) in seqs.iter().enumerate() {
        for p in sage_peptide(s).apply(&cfg.vars, &cfg.statics, cfg.max_var) {
            let mut mask = 0u64;
            for (k, (&c, &m)) in p.sequence.iter().zip(p.modifications.iter()).enumerate() {
                let v = cfg.var_res[(c - b'A') as usize];
                if v != 0.0 && (m as f64 - v).abs() < 1e-3 {
                    mask |= 1 << k;
                }
            }
            if p.nterm.is_some() {
                mask |= NTERM_BIT;
            }
            let b: Vec<f32> = IonSeries::new(&p, Kind::B).map(|x| x.monoisotopic_mass).collect();
            let y: Vec<f32> = IonSeries::new(&p, Kind::Y).map(|x| x.monoisotopic_mass).collect();
            out.push((i, p.to_string(), mask, p.monoisotopic, b, y));
        }
    }
    Ok(out)
}

#[pymodule]
fn vulcan_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(digest_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(expand_pairs, m)?)?;
    m.add_function(wrap_pyfunction!(peptidoform_id, m)?)?;
    m.add_function(wrap_pyfunction!(sage_reference, m)?)?;
    m.add("SAGE_COMMIT", "99407db6e3754b31a9b88b7316a0aee67293c93f")?;
    m.add("RESIDUE_F64", RESIDUE_F64.to_vec())?;
    m.add("H2O_F64", H2O_F64)?;
    m.add("NTERM_BIT", NTERM_BIT)?;
    Ok(())
}
