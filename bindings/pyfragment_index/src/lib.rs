
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::pyclass::CompareOp;
use pyo3::types::{PyFloat, PyType};

use fragment_index_rs::fragment::{Fragment as _Fragment, FragmentSeries};
use fragment_index_rs::index::{SearchIndex};
use fragment_index_rs::interval::Interval as _Interval;
use fragment_index_rs::parent::Peptide as _Peptide;
use fragment_index_rs::sort::{SortType, ParentID};


#[pyclass]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Fragment(_Fragment);

#[pymethods]
impl Fragment {
    #[new]
    fn init(mass: f32, parent_id: ParentID, ordinal: u16, series: &str) -> PyResult<Self> {
        let series = match series {
            "b" => FragmentSeries::b,
            "y" => FragmentSeries::y,
            "c" => FragmentSeries::c,
            "z" => FragmentSeries::z,
            "a" => FragmentSeries::a,
            "x" => FragmentSeries::x,
            _ => return Err(PyValueError::new_err("Unknown fragment series")),
        };
        Ok(Self(_Fragment::new(mass, parent_id, series, ordinal)))
    }

    #[getter]
    fn mass(&self) -> f32 {
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

impl Into<_Fragment> for Fragment {
    fn into(self) -> _Fragment {
        self.0
    }
}

impl From<_Fragment> for Fragment {
    fn from(value: _Fragment) -> Self {
        Self(value)
    }
}

impl AsRef<_Fragment> for Fragment {
    fn as_ref(&self) -> &_Fragment {
        &self.0
    }
}

#[pyclass]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct Interval(_Interval);

impl From<_Interval> for Interval {
    fn from(value: _Interval) -> Self {
        Self(value)
    }
}

impl Into<_Interval> for Interval {
    fn into(self) -> _Interval {
        self.0
    }
}

impl AsRef<_Interval> for Interval {
    fn as_ref(&self) -> &_Interval {
        &self.0
    }
}

#[pymethods]
impl Interval {
    #[new]
    pub fn new(start: usize, end: usize) -> Self {
        _Interval::new(start, end).into()
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

#[pyclass]
#[derive(Debug, Default, Clone, PartialEq)]
#[pyo3(text_signature = "(mass, id, protein_id, start_position, sequence)")]
pub struct Peptide(_Peptide);

impl AsRef<_Peptide> for Peptide {
    fn as_ref(&self) -> &_Peptide {
        &self.0
    }
}

impl From<_Peptide> for Peptide {
    fn from(value: _Peptide) -> Self {
        Self(value)
    }
}

impl Into<_Peptide> for Peptide {
    fn into(self) -> _Peptide {
        self.0
    }
}

#[pymethods]
impl Peptide {
    #[new]
    fn new(mass: f32, id: ParentID, protein_id: ParentID, start_position: u16, sequence: String) -> Self {
        Self(_Peptide::new(mass, id, protein_id, start_position, sequence))
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
            let inner = _Peptide::new(mass, id, protein_id, start_position, proforma.to_string());
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
    fn mass(&self) -> f32 {
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

#[pyclass]
#[derive(Debug)]
pub struct PeptideFragmentIndex(SearchIndex<_Fragment, _Peptide>);

#[pymethods]
impl PeptideFragmentIndex {
    #[new]
    fn new(bins_per_dalton: u32, max_fragment_size: f32) -> PyResult<Self> {
        Ok(Self(SearchIndex::empty(bins_per_dalton, max_fragment_size)))
    }

    fn add_parent(&mut self, parent: Peptide) {
        self.0.add_parent(parent.into())
    }

    fn add(&mut self, fragment: Fragment) {
        self.0.add(fragment.into());
    }

    fn parents_for(&self, mass: f32, error_tolerance: f32) -> Interval {
        self.0.parents_for(mass, error_tolerance).into()
    }

    fn parents_for_range(&self, low: f32, high: f32, error_tolerance: f32) -> Interval {
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

    fn search(&self, query: f32, error_tolerance: f32, parent_interval: Option<Interval>) -> Vec<Fragment> {
        let parent_interval: Option<_Interval> = match parent_interval {
            Some(iv) => {
                Some(iv.as_ref().clone())
            },
            None => None
        };

        let searcher = self.0.search(query, error_tolerance, parent_interval);
        searcher.into_iter().map(|f| f.into()).collect()
    }
}

impl AsRef<SearchIndex<_Fragment, _Peptide>> for PeptideFragmentIndex {
    fn as_ref(&self) -> &SearchIndex<_Fragment, _Peptide> {
        &self.0
    }
}


/// A Python module implemented in Rust.
#[pymodule]
fn pyfragment_index(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Fragment>()?;
    m.add_class::<Interval>()?;
    m.add_class::<Peptide>()?;
    m.add_class::<PeptideFragmentIndex>()?;
    Ok(())
}
