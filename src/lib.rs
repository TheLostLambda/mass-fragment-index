
pub mod sort;
pub mod interval;
pub mod index;
pub mod fragment;
pub mod parent;
pub mod peak;
pub mod r#match;

pub use crate::sort::{IndexSortable, MassType, Tolerance, ToleranceParsingError};
pub use crate::interval::Interval;
pub use crate::index::{SearchIndex};
pub use crate::fragment::{Fragment, FragmentSeriesParsingError};
pub use crate::parent::{ParentMolecule, Peptide, Spectrum};
pub use crate::peak::{DeconvolutedPeak, MZPeak};

pub type PeptideFragmentIndex = SearchIndex<Fragment, Peptide>;
pub type SpectrumIndex = SearchIndex<MZPeak, Spectrum>;
pub type DeconvolutedSpectrumIndex = SearchIndex<DeconvolutedPeak, Spectrum>;