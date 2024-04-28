
use std::sync::{Arc, RwLock};

use pyo3::prelude::*;
use pyo3::exceptions::{PyValueError, PySystemError};
use pyo3::pyclass::CompareOp;
use pyo3::types::{PyFloat, PyType, PyBytes};

use mass_fragment_index::fragment::{Fragment, FragmentSeries};
use mass_fragment_index::index::SearchIndex;
use mass_fragment_index::interval::Interval;
use mass_fragment_index::parent::{Peptide, Spectrum};
use mass_fragment_index::sort::{SortType, ParentID, MassType, Tolerance};
use mass_fragment_index::peak::DeconvolutedPeak;

use mass_fragment_index::storage;

use rmp_serde;


#[pyclass(name="Fragment", module="pyfragment_index")]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[pyo3(text_signature = "(mass, parent_id, ordinal, series)")]
pub struct PyFragment(Fragment);

#[pymethods]
impl PyFragment {
    #[new]
    fn new(mass: MassType, parent_id: ParentID, ordinal: u16, series: &str) -> PyResult<Self> {
        let series = match series {
            "b" => FragmentSeries::b,
            "y" => FragmentSeries::y,
            "c" => FragmentSeries::c,
            "z" => FragmentSeries::z,
            "a" => FragmentSeries::a,
            "x" => FragmentSeries::x,
            _ => return Err(PyValueError::new_err("Unknown fragment series")),
        };
        Ok(Self(Fragment::new(mass, parent_id, series, ordinal)))
    }

    #[getter]
    fn mass(&self) -> MassType {
        self.0.mass
    }

    #[getter]
    fn parent_id(&self) -> ParentID {
        self.0.parent_id
    }

    #[getter]
    fn ordinal(&self) -> u16 {
        self.0.ordinal
    }

    #[getter]
    fn series(&self) -> PyResult<String> {
        let series_label = match self.0.series {
            FragmentSeries::b => "b",
            FragmentSeries::y => "y",
            FragmentSeries::c => "c",
            FragmentSeries::z => "z",
            FragmentSeries::a => "a",
            FragmentSeries::x => "x",
            _ => return Err(PyValueError::new_err("Unknown fragment series")),
        }
        .to_string();
        Ok(series_label)
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "PyFragment({}, {}, {}, {})",
            self.mass(),
            self.parent_id(),
            self.ordinal(),
            self.series()?
        )
        .to_string())
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp, py: Python<'_>) -> PyObject {
        match op {
            CompareOp::Eq => self.0.eq(other.as_ref()).into_py(py),
            CompareOp::Ne => self.0.ne(other.as_ref()).into_py(py),
            _ => py.NotImplemented()
        }
    }
}

impl Into<Fragment> for PyFragment {
    fn into(self) -> Fragment {
        self.0
    }
}

impl From<Fragment> for PyFragment {
    fn from(value: Fragment) -> Self {
        Self(value)
    }
}

impl AsRef<Fragment> for PyFragment {
    fn as_ref(&self) -> &Fragment {
        &self.0
    }
}

#[pyclass(name="Interval", module="pyfragment_index")]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[pyo3(text_signature = "(start, end)")]
pub struct PyInterval(Interval);

impl From<Interval> for PyInterval {
    fn from(value: Interval) -> Self {
        Self(value)
    }
}

impl Into<Interval> for PyInterval {
    fn into(self) -> Interval {
        self.0
    }
}

impl AsRef<Interval> for PyInterval {
    fn as_ref(&self) -> &Interval {
        &self.0
    }
}

#[pymethods]
impl PyInterval {
    #[new]
    pub fn new(start: usize, end: usize) -> Self {
        Interval::new(start, end).into()
    }

    #[getter]
    fn start(&self) -> usize {
        self.0.start
    }

    #[getter]
    fn end(&self) -> usize {
        self.0.end
    }

    fn __contains__(&self, point: usize) -> bool {
        self.contains(point)
    }

    pub fn contains(&self, point: usize) -> bool {
        self.0.contains(point)
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp, py: Python<'_>) -> PyObject {
        match op {
            CompareOp::Eq => self.0.eq(other.as_ref()).into_py(py),
            CompareOp::Ne => self.0.ne(other.as_ref()).into_py(py),
            _ => py.NotImplemented()
        }
    }

    fn __repr__(&self) -> String {
        format!("PyInterval({}, {})", self.start(), self.end()).to_string()
    }
}

#[pyclass(name="Peptide", module="pyfragment_index")]
#[derive(Debug, Default, Clone, PartialEq)]
#[pyo3(text_signature = "(mass, id, protein_id, start_position, sequence)")]
pub struct PyPeptide(Peptide);

impl AsRef<Peptide> for PyPeptide {
    fn as_ref(&self) -> &Peptide {
        &self.0
    }
}

impl From<Peptide> for PyPeptide {
    fn from(value: Peptide) -> Self {
        Self(value)
    }
}

impl Into<Peptide> for PyPeptide {
    fn into(self) -> Peptide {
        self.0
    }
}

#[pymethods]
impl PyPeptide {
    #[new]
    fn new(mass: MassType, id: ParentID, protein_id: ParentID, start_position: u16, sequence: String) -> Self {
        Self(Peptide::new(mass, id, protein_id, start_position, sequence))
    }

    #[staticmethod]
    #[pyo3(signature = (proforma, id=0, protein_id=0, start_position=0))]
    fn from_proforma(
        proforma: &str,
        id: ParentID,
        protein_id: ParentID,
        start_position: u16,
    ) -> PyResult<Self> {
        Python::with_gil(|py| {
            let proforma_module = py.import("pyteomics.proforma")?;
            let proforma_type = proforma_module.getattr("ProForma")?.downcast::<PyType>()?;
            let inst = proforma_type.call_method1("parse", (proforma,))?;
            let mass_obj = inst.getattr("mass")?.downcast::<PyFloat>()?;
            let mass = mass_obj.extract()?;
            let inner = Peptide::new(mass, id, protein_id, start_position, proforma.to_string());
            Ok(inner.into())
        })
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "PyPeptide({}, {}, {}, {})",
            self.mass(),
            self.id(),
            self.protein_id(),
            self.start_position()
        )
        .to_string())
    }

    #[getter]
    fn mass(&self) -> MassType {
        self.0.mass
    }

    #[getter]
    fn id(&self) -> ParentID {
        self.0.id
    }

    #[setter]
    fn set_id(&mut self, value: ParentID) {
        self.0.id = value;
    }

    #[getter]
    fn protein_id(&self) -> ParentID {
        self.0.protein_id
    }

    #[setter]
    fn set_protein_id(&mut self, value: ParentID) {
        self.0.protein_id = value
    }

    #[getter]
    fn start_position(&self) -> u16 {
        self.0.start_position
    }

    #[setter]
    fn set_start_position(&mut self, value: u16) {
        self.0.start_position = value
    }

    #[getter]
    fn sequence(&self) -> String {
        self.0.sequence.clone()
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp, py: Python<'_>) -> PyObject {
        match op {
            CompareOp::Eq => self.0.eq(other.as_ref()).into_py(py),
            CompareOp::Ne => self.0.ne(other.as_ref()).into_py(py),
            _ => py.NotImplemented()
        }
    }
}


#[pyclass(name="Spectrum", module="pyfragment_index")]
#[derive(Debug, Default, Clone, PartialEq)]
#[pyo3(text_signature = "(scan_number, source_file_id, precursor_mass, precursor_charge)")]
pub struct PySpectrum(Spectrum);

impl AsRef<Spectrum> for PySpectrum {
    fn as_ref(&self) -> &Spectrum {
        &self.0
    }
}

impl From<Spectrum> for PySpectrum {
    fn from(value: Spectrum) -> Self {
        Self(value)
    }
}

impl Into<Spectrum> for PySpectrum {
    fn into(self) -> Spectrum {
        self.0
    }
}

#[pymethods]
impl PySpectrum {
    #[new]
    fn new(scan_number: ParentID,  source_file_id: ParentID, precursor_mass: MassType, precursor_charge: i32) -> PySpectrum {
        Spectrum {
            scan_number,
            source_file_id,
            precursor_mass,
            precursor_charge,
            sort_id: u32::MAX
        }.into()
    }

    #[getter]
    fn precursor_mass(&self) -> MassType {
        self.0.precursor_mass
    }

    #[getter]
    fn scan_number(&self) -> ParentID {
        self.0.scan_number
    }

    #[setter]
    fn set_scan_number(&mut self, value: ParentID) {
        self.0.scan_number = value;
    }

    #[getter]
    fn source_file_id(&self) -> ParentID {
        self.0.source_file_id
    }

    #[setter]
    fn set_source_file_id(&mut self, value: ParentID) {
        self.0.source_file_id = value
    }

    #[getter]
    fn precursor_charge(&self) -> i32 {
        self.0.precursor_charge
    }

    #[setter]
    fn set_precursor_charge(&mut self, value: i32) {
        self.0.precursor_charge = value
    }

    fn __repr__(&self) -> String {
        format!(
            "Spectrum({}, {}, {}, {})",
            self.0.scan_number,
            self.0.source_file_id,
            self.0.precursor_mass,
            self.0.precursor_charge
        ).to_string()
    }
}


#[pyclass(name="Peak", module="pyfragment_index")]
#[derive(Debug, Default, Clone, PartialEq)]
#[pyo3(text_signature = "(mass, intensity, charge, scan_ref)")]
pub struct PyPeak(DeconvolutedPeak);

impl AsRef<DeconvolutedPeak> for PyPeak {
    fn as_ref(&self) -> &DeconvolutedPeak {
        &self.0
    }
}

impl From<DeconvolutedPeak> for PyPeak {
    fn from(value: DeconvolutedPeak) -> Self {
        Self(value)
    }
}

impl Into<DeconvolutedPeak> for PyPeak {
    fn into(self) -> DeconvolutedPeak {
        self.0
    }
}


#[pymethods]
impl PyPeak {

    #[new]
    fn new(mass: MassType, intensity: f32, charge: i16, scan_ref: ParentID) -> PyPeak {
        DeconvolutedPeak {
            mass, intensity, charge, scan_ref
        }.into()
    }

    #[getter]
    fn get_mass(&self) -> MassType {
        self.0.mass
    }

    #[getter]
    fn get_intensity(&self) -> f32 {
        self.0.intensity
    }

    #[getter]
    fn get_charge(&self) -> i16 {
        self.0.charge
    }

    #[getter]
    fn get_scan_ref(&self) -> ParentID {
        self.0.scan_ref
    }

    #[setter]
    fn set_intensity(&mut self, value: f32) {
        self.0.intensity = value;
    }

    #[setter]
    fn set_charge(&mut self, value: i16) {
        self.0.charge = value;
    }

    #[setter]
    fn set_scan_ref(&mut self, value: ParentID) {
        self.0.scan_ref = value;
    }

    fn __repr__(&self) -> String {
        format!("Peak({}, {}, {}, {})", self.0.mass, self.0.intensity, self.0.charge, self.0.scan_ref).to_string()
    }

}


#[pyclass(name="PeptideFragmentIndex", module="pyfragment_index")]
#[derive(Debug)]
pub struct PyPeptideFragmentIndex(SearchIndex<Fragment, Peptide>);

#[pymethods]
impl PyPeptideFragmentIndex {
    #[new]
    fn new(bins_per_dalton: u32, maxfragment_size: MassType) -> PyResult<Self> {
        Ok(Self(SearchIndex::empty(bins_per_dalton, maxfragment_size)))
    }

    fn add_parent(&mut self, parent: PyPeptide) {
        self.0.add_parent(parent.into())
    }

    fn add(&mut self, fragment: PyFragment) {
        self.0.add(fragment.into());
    }

    fn parents_for(&self, mass: MassType, error_tolerance: MassType) -> PyInterval {
        let error_tolerance = Tolerance::PPM(error_tolerance);
        self.0.parents_for(mass, error_tolerance).into()
    }

    fn parents_for_range(&self, low: MassType, high: MassType, error_tolerance: MassType) -> PyInterval {
        let error_tolerance = Tolerance::PPM(error_tolerance);
        self.0.parents_for_range(low, high, error_tolerance).into()
    }

    fn __repr__(&self) -> String {
        format!(
            "PyPeptideFragmentIndex({} fragments, {} peptides, sorted={})",
            self.0.bins.iter().map(|b| b.len()).sum::<usize>(),
            self.0.parents.len(),
            matches!(self.0.sort_type, SortType::ByParentId)
        )
        .to_string()
    }

    fn sort(&mut self) {
        self.0.sort(SortType::ByParentId)
    }

    fn search(&self, query: MassType, error_tolerance: MassType, parentinterval: Option<PyInterval>) -> Vec<PyFragment> {
        let error_tolerance = Tolerance::PPM(error_tolerance);
        let parentinterval: Option<Interval> = match parentinterval {
            Some(iv) => {
                Some(iv.as_ref().clone())
            },
            None => None
        };

        let searcher = self.0.search(query, error_tolerance, parentinterval);
        searcher.into_iter().map(|f| f.into()).collect()
    }

    pub fn __getstate__<'py>(&self, py: Python<'py>) -> PyResult<&'py PyBytes> {
        match rmp_serde::to_vec(&self.0) {
            Ok(payload) => {
                match zstd::bulk::compress(&payload, 0) {
                    Ok(payload) => Ok(PyBytes::new(py, &payload)),
                    Err(err) => Err(PyValueError::new_err(err.to_string()))
                }
            },
            Err(err) => {
                Err(PyValueError::new_err(err.to_string()))
            }
        }

    }

    pub fn __setstate__(&mut self, state: &PyBytes) -> PyResult<()> {
        let buf = state.as_bytes();
        match zstd::Decoder::new(std::io::BufReader::new(buf)) {
            Ok(decoder) => {
                match rmp_serde::from_read(decoder) {
                    Ok(inst) => {
                        self.0 = inst;
                        Ok(())
                    },
                    Err(err) => {
                        Err(PyValueError::new_err(err.to_string()))
                    }
                }
            },
            Err(err) => {
                Err(PyValueError::new_err(err.to_string()))
            }
        }
    }

    pub fn __getnewargs__(&self) -> PyResult<(u32, MassType)> {
        Ok((self.0.bins_per_dalton, self.0.max_item_mass))
    }

}

impl AsRef<SearchIndex<Fragment, Peptide>> for PyPeptideFragmentIndex {
    fn as_ref(&self) -> &SearchIndex<Fragment, Peptide> {
        &self.0
    }
}


#[pyclass(name="PeakIndex", module="pyfragment_index")]
#[derive(Debug)]
pub struct PyPeakIndex(Arc<RwLock<SearchIndex<DeconvolutedPeak, Spectrum>>>);

#[pymethods]
impl PyPeakIndex {
    #[new]
    fn new(bins_per_dalton: u32, maxfragment_size: MassType) -> PyResult<Self> {
        Ok(Self(Arc::new(RwLock::new(SearchIndex::empty(bins_per_dalton, maxfragment_size)))))
    }

    pub fn __getstate__<'py>(&self, py: Python<'py>) -> PyResult<&'py PyBytes> {
        match rmp_serde::to_vec(&self.0.as_ref()) {
            Ok(payload) => {
                match zstd::bulk::compress(&payload, 0) {
                    Ok(payload) => Ok(PyBytes::new(py, &payload)),
                    Err(err) => Err(PyValueError::new_err(err.to_string()))
                }
            },
            Err(err) => {
                Err(PyValueError::new_err(err.to_string()))
            }
        }

    }

    pub fn __setstate__(&mut self, state: &PyBytes) -> PyResult<()> {
        let buf = state.as_bytes();
        match zstd::Decoder::new(std::io::BufReader::new(buf)) {
            Ok(decoder) => {
                match rmp_serde::from_read(decoder) {
                    Ok(inst) => {
                        self.0 = Arc::new(inst);
                        Ok(())
                    },
                    Err(err) => {
                        Err(PyValueError::new_err(err.to_string()))
                    }
                }
            },
            Err(err) => {
                Err(PyValueError::new_err(err.to_string()))
            }
        }
    }

    pub fn __getnewargs__(&self) -> PyResult<(u32, MassType)> {
        Ok((self.0.read().unwrap().bins_per_dalton, self.0.read().unwrap().max_item_mass))
    }

    fn add_parent(&mut self, parent: PySpectrum) {
        self.0.write().unwrap().add_parent(parent.into())
    }

    fn add(&mut self, peak: PyPeak) {
        self.0.write().unwrap().add(peak.into());
    }

    fn parents_for(&self, mass: MassType, error_tolerance: MassType) -> PyInterval {
        let error_tolerance = Tolerance::PPM(error_tolerance);
        self.0.read().unwrap().parents_for(mass, error_tolerance).into()
    }

    fn parents_for_range(&self, low: MassType, high: MassType, error_tolerance: MassType) -> PyInterval {
        let error_tolerance = Tolerance::PPM(error_tolerance);
        self.0.read().unwrap().parents_for_range(low, high, error_tolerance).into()
    }

    fn __repr__(&self) -> String {
        format!(
            "PeakIndex({} peaks, {} spectra, sorted={})",
            self.0.read().unwrap().bins.iter().map(|b| b.len()).sum::<usize>(),
            self.0.read().unwrap().parents.len(),
            matches!(self.0.read().unwrap().sort_type, SortType::ByParentId)
        )
        .to_string()
    }

    fn sort(&mut self) {
        self.0.write().unwrap().sort(SortType::ByParentId)
    }

    fn __iter__(&self) -> PeakIndexIter {
        PeakIndexIter::new(self.0.clone())
    }

    fn search(&self, query: MassType, error_tolerance: MassType, parent_interval: Option<PyInterval>) -> PyResult<Vec<PyPeak>> {
        let error_tolerance = Tolerance::PPM(error_tolerance);
        let parent_interval: Option<Interval> = match parent_interval {
            Some(iv) => {
                Some(iv.as_ref().clone())
            },
            None => None
        };

        match self.0.read() {
            Ok(obj) => {
                let searcher = obj.search(query, error_tolerance, parent_interval);
                Ok(searcher.into_iter().map(|f| f.into()).collect())
            },
            Err(err) => Err(PySystemError::new_err(format!("Lock contention: {}", err))),
        }
    }
}


#[pyclass]
pub struct PeakIndexIter {
    search_index: Arc<RwLock<SearchIndex<DeconvolutedPeak, Spectrum>>>,
    bin_index: usize,
    position_in_bin: usize
}

impl PeakIndexIter {
    fn new(search_index: Arc<RwLock<SearchIndex<DeconvolutedPeak, Spectrum>>>) -> Self {
        PeakIndexIter {
            search_index: search_index,
            bin_index: 0,
            position_in_bin: 0
        }
    }
}

#[pymethods]
impl PeakIndexIter {
    fn __len__(&self) -> usize {
        self.search_index.read().map(|index| {
            index.bins.iter().map(|b| b.len()).sum::<usize>()
        }).unwrap()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<PyPeak> {
        let val = self.search_index.read().map(|index| {
            if self.bin_index >= index.bins.len() {
                None
            } else {
                let bin = &index.bins[self.bin_index];
                if self.position_in_bin < bin.len() {
                    let value = Some(bin[self.position_in_bin].into());
                    self.position_in_bin += 1;
                    value
                } else {
                    self.position_in_bin = 0;
                    while self.bin_index < index.bins.len() {
                        self.bin_index += 1;
                        let bin = &index.bins[self.bin_index];
                        if bin.len() > 0 {
                            let value = Some(bin[self.position_in_bin].into());
                            self.position_in_bin += 1;
                            return value
                        }
                    }
                    return None
                }
            }
        });
        val.unwrap_or_default()
    }
}


#[pyfunction]
fn py_read_peak_index_pq(path: &str) -> PyResult<PyPeakIndex> {
    let index = storage::read_peak_index(&path)?;
    Ok(PyPeakIndex(Arc::new(RwLock::new(index))))
}

/// A Python module implemented in Rust.
#[pymodule]
fn pyfragment_index(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyInterval>()?;

    m.add_class::<PyFragment>()?;
    m.add_class::<PyPeptide>()?;
    m.add_class::<PyPeptideFragmentIndex>()?;

    m.add_class::<PyPeak>()?;
    m.add_class::<PySpectrum>()?;
    m.add_class::<PyPeakIndex>()?;

    m.add_function(wrap_pyfunction!(py_read_peak_index_pq, m)?)?;
    Ok(())
}
