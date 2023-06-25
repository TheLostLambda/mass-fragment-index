use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use fragment_index_rs::fragment::{Fragment, FragmentSeries};
use fragment_index_rs::interval::Interval;


#[pyclass]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PyFragment(Fragment);

#[pymethods]
impl PyFragment {

    #[new]
    fn init(mass: f32, parent_id: usize, ordinal: u16, series: &str) -> PyResult<Self> {
        let series = match series {
            "b" => FragmentSeries::b,
            "y" => FragmentSeries::y,
            "c" => FragmentSeries::c,
            "z" => FragmentSeries::z,
            "a" => FragmentSeries::a,
            "x" => FragmentSeries::x,
            _ => {
                return Err(PyValueError::new_err("Unknown fragment series"))
            }
        };
        Ok(
            Self(
                Fragment::new(mass, parent_id, series, ordinal)
            )
        )
    }

    #[getter]
    fn mass(&self) -> f32 {
        self.0.mass
    }

    #[getter]
    fn parent_id(&self) -> usize {
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
            _ => {
                return Err(PyValueError::new_err("Unknown fragment series"))
            }
        }.to_string();
        Ok(series_label)
    }

    fn __repr__(&self) -> String {
        format!("Py{:?}", self.0).to_string()
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


#[pyclass]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
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

    fn __repr__(&self) -> String {
        format!("Py{:?}", self.0).to_string()
    }
}


/// A Python module implemented in Rust.
#[pymodule]
fn pyfragment_index(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyFragment>()?;
    m.add_class::<PyInterval>()?;
    Ok(())
}