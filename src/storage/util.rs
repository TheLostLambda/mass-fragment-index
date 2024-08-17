use std::{io, path::Path};

use arrow::{array::RecordBatch, datatypes::Schema};
use parquet::basic::Compression;


pub trait IndexBinaryStorage<T, P> {
    fn write<D: AsRef<Path>>(
        &self,
        directory: &D,
        compression_level: Option<Compression>
    );

    fn read<D: AsRef<Path>>(directory: &D) -> io::Result<Self> where Self: Sized;

    fn entries_name() -> String;

    fn parents_name() -> String;

    fn make_parent_schema() -> Schema;

    fn make_entry_schema() -> Schema;

    fn make_metadata_schema() -> Schema;

    fn parents_to_arrow(&self, entries: &[P], schema: &Schema) -> RecordBatch;

    fn entries_to_arrow(&self, entries: &[T], schema: &Schema, segment_id: u64) -> RecordBatch;

    fn make_metadata(&self) -> RecordBatch;
}
