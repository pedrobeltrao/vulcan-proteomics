#!/usr/bin/env python3
"""Fetch reference FASTAs into data/fasta/ (gitignored). Stdlib only.

UniProt entries come from the REST stream endpoint (format=fasta, canonical sequences only,
no isoforms); the X-UniProt-Release / X-UniProt-Release-Date response headers are recorded.
Human one-protein-per-gene comes from the reference_proteomes FTP (release from RELEASE.metalink).
Contaminants: Frankenfield et al. 2022 universal contaminant library, pinned to a git commit.
Idempotent: existing files are kept (UniProt 'current' moves on; re-fetching would change the
database) unless --force.  Writes data/manifests/fasta.tsv (provenance) and fasta.sha256.
"""
import argparse, gzip, hashlib, os, re, sys, time, urllib.parse, urllib.request, datetime

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "data", "fasta")
MAN = os.path.join(ROOT, "data", "manifests")
UA = {"User-Agent": "fetch-fasta/1.0 (" + os.environ.get("VULCAN_CONTACT_EMAIL", "noreply@example.com") + " research)"}
REST = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&compressed=false&query="
RP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
HAO = ("https://raw.githubusercontent.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/"
       "bddf6da9334f4300caed75f7a87c096dd96b6bca/Universal%20protein%20contaminant%20FASTA/0602_Universal%20Contaminants.fasta")

# name, kind, spec, description, license
ITEMS = [
    ("human_sp_canonical", "rest", "(proteome:UP000005640) AND (reviewed:true)", "Human UP000005640, Swiss-Prot reviewed canonical", "CC BY 4.0"),
    ("human_proteome_canonical", "rest", "proteome:UP000005640", "Human UP000005640 all entries (Swiss-Prot+TrEMBL), canonical", "CC BY 4.0"),
    ("human_opg", "ftp", RP + "Eukaryota/UP000005640/UP000005640_9606.fasta.gz", "Human UP000005640 one-protein-per-gene (reference_proteomes FTP)", "CC BY 4.0"),
    ("pfu_UP000001013", "rest", "proteome:UP000001013", "Pyrococcus furiosus DSM 3638", "CC BY 4.0"),
    ("ecoli_K12_UP000000625", "rest", "proteome:UP000000625", "Escherichia coli K12", "CC BY 4.0"),
    ("yeast_UP000002311", "rest", "proteome:UP000002311", "Saccharomyces cerevisiae S288C", "CC BY 4.0"),
    ("arath_UP000006548", "rest", "proteome:UP000006548", "Arabidopsis thaliana (entrapment foreign proteome)", "CC BY 4.0"),
    ("contaminants_frankenfield2022", "url", HAO, "Universal protein contaminants, Frankenfield et al. J Proteome Res 2022 (HaoGroup-ProtContLib @bddf6da)", "see repo (no LICENSE file; paper CC BY)"),
]

def get(url, retries=5):
    for i in range(retries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300)
        except Exception as e:
            print("  retry", i + 1, e, file=sys.stderr); time.sleep(4 * (i + 1))
    raise SystemExit("failed: " + url)

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True); os.makedirs(MAN, exist_ok=True)
    prov_p = os.path.join(MAN, "fasta.tsv")
    prov = {}
    if os.path.exists(prov_p):
        for l in open(prov_p).read().splitlines()[1:]:
            f = l.split("\t"); prov[f[0]] = f
    for name, kind, spec, desc, lic in ITEMS:
        path = os.path.join(OUT, name + ".fasta")
        if os.path.exists(path) and name in prov and not a.force:
            print("[skip]", name); continue
        print("[get ]", name, flush=True)
        rel = reldate = ""
        if kind == "rest":
            url = REST + urllib.parse.quote(spec)
            r = get(url); rel = r.headers.get("X-UniProt-Release", ""); reldate = r.headers.get("X-UniProt-Release-Date", "")
            data = r.read()
        elif kind == "ftp":
            url = spec; data = gzip.decompress(get(url).read())
            m = get(spec.rsplit("/", 1)[0] + "/RELEASE.metalink").read().decode()
            rel = re.search(r"<version>([^<]+)</version>", m).group(1)
        else:
            url = spec; data = get(url).read(); rel = "git:bddf6da9334f4300caed75f7a87c096dd96b6bca"
        with open(path + ".tmp", "wb") as f: f.write(data)
        os.replace(path + ".tmp", path)
        n = data.count(b"\n>") + (1 if data.startswith(b">") else 0)
        prov[name] = [name, kind, url, rel, reldate, datetime.date.today().isoformat(), str(n), str(len(data)), sha(path), lic, desc]
    with open(prov_p, "w") as f:
        f.write("name\tkind\turl\trelease\trelease_date\tretrieved\tn_sequences\tbytes\tsha256\tlicense\tdescription\n")
        for it in ITEMS:
            if it[0] in prov: f.write("\t".join(prov[it[0]]) + "\n")
    with open(os.path.join(MAN, "fasta.sha256"), "w") as f:
        for it in ITEMS:
            if it[0] in prov: f.write(f"{prov[it[0]][8]}  fasta/{it[0]}.fasta\n")
    print("wrote", prov_p)

if __name__ == "__main__":
    main()
