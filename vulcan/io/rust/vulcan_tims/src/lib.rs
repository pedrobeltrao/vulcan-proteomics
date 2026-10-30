//! Vulcan L0: numpy-returning Bruker TDF frame reader (wraps timsrust 0.4).
//! `TdfFrames(path).read(start, stop)` decodes frames [start, stop) and returns flat numpy
//! arrays (no per-peak Python objects, unlike timsrust_pyo3 which returns Python lists).
use numpy::IntoPyArray;
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use timsrust::readers::FrameReader;
use timsrust::MSLevel;

#[pyclass]
struct TdfFrames {
    reader: FrameReader,
}

#[pymethods]
impl TdfFrames {
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        FrameReader::new(path)
            .map(|r| Self { reader: r })
            .map_err(|e| PyIOError::new_err(format!("{:?}", e)))
    }

    fn __len__(&self) -> usize {
        self.reader.len()
    }

    /// Decode frames [start, stop). Returns a dict of numpy arrays:
    /// frame_index (u64, TDF frame id), rt_s (f64), ms_level (u8: 1, 2, 0=unknown),
    /// event_offsets (u64, n+1), scan_offsets (u64, concatenated per-frame scan offsets, each
    /// relative to its frame), scan_ptr (u64, n+1 into scan_offsets), tof (u32), intensity (u32),
    /// quad_ptr (u64, n+1), quad_scan_start/end (u64), quad_mz/width/ce (f64).
    fn read(&self, py: Python<'_>, start: usize, stop: usize) -> PyResult<PyObject> {
        let stop = stop.min(self.reader.len());
        let n = stop.saturating_sub(start);
        let mut frame_index: Vec<u64> = Vec::with_capacity(n);
        let mut rt: Vec<f64> = Vec::with_capacity(n);
        let mut lvl: Vec<u8> = Vec::with_capacity(n);
        let mut ev_off: Vec<u64> = Vec::with_capacity(n + 1);
        let mut scan_ptr: Vec<u64> = Vec::with_capacity(n + 1);
        let mut quad_ptr: Vec<u64> = Vec::with_capacity(n + 1);
        let mut scan_offsets: Vec<u64> = Vec::new();
        let mut tof: Vec<u32> = Vec::new();
        let mut inten: Vec<u32> = Vec::new();
        let (mut qs, mut qe, mut qmz, mut qw, mut qce) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        ev_off.push(0);
        scan_ptr.push(0);
        quad_ptr.push(0);
        for i in start..stop {
            let f = self.reader.get(i).map_err(|e| PyIOError::new_err(format!("frame {}: {:?}", i, e)))?;
            frame_index.push(f.index as u64);
            rt.push(f.rt_in_seconds);
            lvl.push(match f.ms_level {
                MSLevel::MS1 => 1,
                MSLevel::MS2 => 2,
                _ => 0,
            });
            tof.extend_from_slice(&f.tof_indices);
            inten.extend_from_slice(&f.intensities);
            ev_off.push(tof.len() as u64);
            scan_offsets.extend(f.scan_offsets.iter().map(|&x| x as u64));
            scan_ptr.push(scan_offsets.len() as u64);
            if f.ms_level == MSLevel::MS2 {
                let q = &f.quadrupole_settings;
                qs.extend(q.scan_starts.iter().map(|&x| x as u64));
                qe.extend(q.scan_ends.iter().map(|&x| x as u64));
                qmz.extend_from_slice(&q.isolation_mz);
                qw.extend_from_slice(&q.isolation_width);
                qce.extend_from_slice(&q.collision_energy);
            }
            quad_ptr.push(qs.len() as u64);
        }
        let d = PyDict::new_bound(py);
        d.set_item("frame_index", frame_index.into_pyarray_bound(py))?;
        d.set_item("rt_s", rt.into_pyarray_bound(py))?;
        d.set_item("ms_level", lvl.into_pyarray_bound(py))?;
        d.set_item("event_offsets", ev_off.into_pyarray_bound(py))?;
        d.set_item("scan_offsets", scan_offsets.into_pyarray_bound(py))?;
        d.set_item("scan_ptr", scan_ptr.into_pyarray_bound(py))?;
        d.set_item("tof", tof.into_pyarray_bound(py))?;
        d.set_item("intensity", inten.into_pyarray_bound(py))?;
        d.set_item("quad_ptr", quad_ptr.into_pyarray_bound(py))?;
        d.set_item("quad_scan_start", qs.into_pyarray_bound(py))?;
        d.set_item("quad_scan_end", qe.into_pyarray_bound(py))?;
        d.set_item("quad_mz", qmz.into_pyarray_bound(py))?;
        d.set_item("quad_width", qw.into_pyarray_bound(py))?;
        d.set_item("quad_ce", qce.into_pyarray_bound(py))?;
        Ok(d.into())
    }
}

#[pymodule]
fn vulcan_tims(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<TdfFrames>()?;
    Ok(())
}
