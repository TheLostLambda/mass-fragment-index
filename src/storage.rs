#![cfg(feature = "binary_storage")]

mod peak_parquet;
mod fragment_parquet;

pub use peak_parquet::{read_peak_index, write_peak_index};
pub use fragment_parquet::{read_fragment_index, write_fragment_index};
