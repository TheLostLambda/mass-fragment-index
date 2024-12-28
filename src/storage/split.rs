use std::{
    collections::HashMap,
    fs, io,
    path::{Path, PathBuf},
    sync::Arc,
};

use arrow::{
    array::{
        ArrayRef, Float32Array, Float32Builder, RecordBatch, StructArray,
        UInt32Array, UInt32Builder,
    },
    datatypes::{DataType, Field, Schema, SchemaBuilder},
    json::{LineDelimitedWriter, ReaderBuilder as JSONReaderBuilder},
};
use itertools::izip;
use parquet::{
    arrow::{arrow_reader::ArrowReaderBuilder, ArrowWriter},
    basic::{Compression, ZstdLevel},
};

use crate::{sort::ParentID, IndexSortable, Interval, MassType};

use super::{
    util::{afield, as_array_ref},
    ArrowStorage, IndexBinaryStorage,
};

#[allow(unused)]
pub trait SplitArrowStorage: ArrowStorage {
    fn split_archive_name_prefix() -> String {
        let arch_base = Self::archive_name();
        if let Some((prefix, _suffix)) = arch_base.rsplit_once('.') {
            prefix.into()
        } else {
            arch_base
        }
    }

    fn split_archive_name_for(segment: u32) -> String {
        let arch_base = Self::archive_name();
        if let Some((prefix, suffix)) = arch_base.rsplit_once('.') {
            format!("{prefix}_{segment}.{suffix}")
        } else {
            format!("{arch_base}_{segment}")
        }
    }

    fn collect_segments_from_path(path: impl AsRef<Path>) -> io::Result<Vec<PathBuf>> {
        let path = path.as_ref();
        let content = path.read_dir()?;
        let prefix = Self::split_archive_name_prefix();

        let mut segment_files = Vec::new();
        for ch in content {
            let ch = ch?;
            if ch.file_name().to_string_lossy().starts_with(&prefix) {
                segment_files.push(ch.path());
            }
        }

        Ok(segment_files)
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct SplitBand {
    pub band_id: ParentID,
    pub start_id: ParentID,
    pub end_id: ParentID,
    pub start_mass: MassType,
    pub end_mass: MassType,
}

impl SplitBand {
    pub fn new(
        band_id: ParentID,
        start_id: ParentID,
        end_id: ParentID,
        start_mass: MassType,
        end_mass: MassType,
    ) -> Self {
        Self {
            band_id,
            start_id,
            end_id,
            start_mass,
            end_mass,
        }
    }
}

impl ArrowStorage for SplitBand {
    fn schema() -> arrow::datatypes::SchemaRef {
        let band_id = afield!("band_id", DataType::UInt32);
        let start_id = afield!("start_id", DataType::UInt32);
        let end_id = afield!("end_id", DataType::UInt32);
        let start_mass = afield!("start_mass", DataType::Float32);
        let end_mass = afield!("end_mass", DataType::Float32);
        Arc::new(Schema::new(vec![
            band_id, start_id, end_id, start_mass, end_mass,
        ]))
    }

    fn from_batch<'a>(
        batch: &'a RecordBatch,
        _schema: arrow::datatypes::SchemaRef,
    ) -> impl Iterator<Item = (Self, u64)> + 'a {
        let band_ids: &UInt32Array = batch.column(0).as_any().downcast_ref().unwrap();
        let start_ids: &UInt32Array = batch.column(1).as_any().downcast_ref().unwrap();
        let end_ids: &UInt32Array = batch.column(2).as_any().downcast_ref().unwrap();
        let start_masses: &Float32Array = batch.column(3).as_any().downcast_ref().unwrap();
        let end_masses: &Float32Array = batch.column(4).as_any().downcast_ref().unwrap();

        izip!(band_ids, start_ids, end_ids, start_masses, end_masses).map(
            |(band_id, start_id, end_id, start_mass, end_mass)| {
                (
                    SplitBand::new(
                        band_id.unwrap(),
                        start_id.unwrap(),
                        end_id.unwrap(),
                        start_mass.unwrap(),
                        end_mass.unwrap(),
                    ),
                    0,
                )
            },
        )
    }

    fn to_batch(
        batch: &[Self],
        schema: arrow::datatypes::SchemaRef,
        _segment_id: u64,
    ) -> Result<RecordBatch, arrow::error::ArrowError> {
        let mut band_ids = UInt32Builder::new();
        let mut start_ids = UInt32Builder::new();
        let mut end_ids = UInt32Builder::new();
        let mut start_masses = Float32Builder::new();
        let mut end_masses = Float32Builder::new();

        for item in batch {
            band_ids.append_value(item.band_id);
            start_ids.append_value(item.start_id);
            end_ids.append_value(item.end_id);
            start_masses.append_value(item.start_mass);
            end_masses.append_value(item.end_mass);
        }

        let columns = vec![
            as_array_ref!(band_ids),
            as_array_ref!(start_ids),
            as_array_ref!(end_ids),
            as_array_ref!(start_masses),
            as_array_ref!(end_masses),
        ];

        let batch = RecordBatch::try_new(schema, columns);
        batch
    }

    fn archive_name() -> String {
        "split_bands.jsonl".into()
    }

    fn writer_properties() -> parquet::file::properties::WriterPropertiesBuilder {
        parquet::file::properties::WriterProperties::builder()
    }
}

pub trait SplitIndexBinaryStorage<
    'a,
    T: ArrowStorage + 'a + IndexSortable + Clone,
    P: ArrowStorage + IndexSortable,
    M: ArrowStorage,
>: IndexBinaryStorage<'a, T, P, M>
{
    fn make_item_schema() -> Arc<Schema> {
        let base_schema = T::schema();
        let mut fields = SchemaBuilder::from(base_schema.fields.clone());
        fields.push(afield!("band_id", DataType::UInt32));
        *fields.metadata_mut() = base_schema.metadata.clone();
        Arc::new(fields.finish())
    }

    fn compute_parent_bands(&self, bin_width: MassType) -> Vec<SplitBand> {
        let parents = self.parents();
        let mut min_value = parents.first().map(|p| p.mass()).unwrap_or_default();
        let mut max_value = min_value + bin_width;
        let mut min_index = 0usize;
        let mut bands = Vec::new();
        for (i, val) in parents.iter().enumerate() {
            if val.mass() > max_value {
                bands.push(SplitBand::new(
                    bands.len() as ParentID,
                    min_index as ParentID,
                    i.saturating_sub(1) as ParentID,
                    min_value,
                    max_value,
                ));
                min_index = i;
                min_value = val.mass();
                max_value = min_value + bin_width;
            }
        }
        bands.push(SplitBand::new(
            bands.len() as ParentID,
            min_index as u32,
            parents.len() as u32,
            min_value,
            max_value,
        ));
        bands
    }

    fn write_entries_split(
        &'a self,
        directory: &Path,
        bands: &[SplitBand],
        compression_level: &Compression,
    ) -> io::Result<()> {
        let entries_path = directory.join(T::archive_name());
        let entries_schema = T::schema();
        let props = T::writer_properties()
            .set_compression(compression_level.clone())
            .build();
        let ext_schema = Self::make_item_schema();
        let mut writer = ArrowWriter::try_new(
            fs::File::create(entries_path)?,
            ext_schema.clone(),
            Some(props),
        )?;
        for band in bands.iter() {
            let interval = Interval::new(band.start_id as usize, band.end_id as usize + 1);
            for (i, bin) in self.iter_entries().enumerate() {
                let entries_of: Vec<T> = bin
                    .iter()
                    .filter(|b| interval.contains(b.parent_id() as usize))
                    .cloned()
                    .collect();

                let batch = T::to_batch(&entries_of, entries_schema.clone(), i as u64).unwrap();
                let (_fields, mut arrays, _null_buffer) = StructArray::from(batch).into_parts();
                let band_id_col = vec![band.band_id; entries_of.len()];
                let band_id_col = Arc::new(UInt32Array::from(band_id_col));
                arrays.push(band_id_col);
                let batch = RecordBatch::try_new(ext_schema.clone(), arrays).unwrap();

                writer.write(&batch)?;
            }
        }
        writer.close()?;
        Ok(())
    }

    fn write_split<D: AsRef<Path>>(
        &'a self,
        directory: &D,
        bin_width: MassType,
        compression_level: Option<Compression>,
    ) -> io::Result<()> {
        let directory = directory.as_ref();

        let compression_level =
            compression_level.unwrap_or_else(|| Compression::ZSTD(ZstdLevel::try_new(9).unwrap()));
        let bands = self.compute_parent_bands(bin_width);

        self.write_metadata(directory)?;
        self.write_parents(directory, &compression_level)?;
        self.write_entries_split(directory, &bands, &compression_level)?;
        self.write_split_log(directory, &bands)?;
        Ok(())
    }

    fn band_log_name() -> String {
        let prefix = "search_bands.json".to_string();
        prefix
    }

    fn write_split_log(&self, directory: &Path, split_log: &[SplitBand]) -> io::Result<()> {
        let split_log_path = directory.join(Self::band_log_name());

        let split_log_fh = io::BufWriter::new(fs::File::create(split_log_path)?);
        let mut writer = LineDelimitedWriter::new(split_log_fh);
        let split_log = SplitBand::to_batch(split_log, SplitBand::schema(), 0).unwrap();

        writer.write(&split_log).unwrap();
        writer.finish().unwrap();
        Ok(())
    }

    fn read_split<D: AsRef<Path>>(directory: &D) -> io::Result<Self>
    where
        Self: Sized,
    {
        let root = directory.as_ref();

        let metadata = Self::read_metadata(root)?;
        let parents = Self::read_parents(root)?;
        // let split_log = Self::read_split_log(root)?;

        let entries = {
            let mut bin_collector: HashMap<u64, Vec<T>> = HashMap::default();
            let entries_fh = fs::File::open(root.join(T::archive_name()))?;
            let reader = ArrowReaderBuilder::try_new(entries_fh)?.build()?;
            let entry_schema = T::schema();

            for batch in reader {
                for (entry, segment_id) in T::from_batch(&batch.unwrap(), entry_schema.clone()) {
                    bin_collector.entry(segment_id).or_default().push(entry);
                }
            }

            bin_collector
        };

        let this = Self::from_components(metadata, parents, entries);
        Ok(this)
    }

    fn read_split_log(directory: &Path) -> io::Result<Vec<SplitBand>> {
        let split_log_path = directory.join(Self::band_log_name());
        let split_log_fh = io::BufReader::new(fs::File::open(split_log_path)?);

        let split_log = JSONReaderBuilder::new(SplitBand::schema())
            .build(split_log_fh)
            .unwrap()
            .map(|batch| {
                let batch: Vec<_> = SplitBand::from_batch(&batch.unwrap(), SplitBand::schema())
                    .map(|(x, _)| x)
                    .collect();
                batch
            })
            .flatten()
            .collect();
        Ok(split_log)
    }

    fn read_parents(directory: &Path) -> io::Result<Vec<P>> {
        let parents_path = directory.join(P::archive_name());
        let parent_schema = P::schema();
        let parents_fh = fs::File::open(parents_path)?;

        let reader = ArrowReaderBuilder::try_new(parents_fh)?.build()?;
        let mut parents = Vec::new();
        for batch in reader {
            parents.extend(P::from_batch(&batch.unwrap(), parent_schema.clone()).map(|(p, _)| p));
        }

        Ok(parents)
    }

    fn read_metadata(directory: &Path) -> io::Result<M> {
        let meta_path = directory.join(M::archive_name());

        let metadata = {
            let meta_schema = M::schema();
            let meta_fh = io::BufReader::new(fs::File::open(meta_path)?);
            let meta_rec = JSONReaderBuilder::new(meta_schema.clone())
                .build(meta_fh)
                .unwrap()
                .next()
                .unwrap()
                .unwrap();

            let (metadata, _) = M::from_batch(&meta_rec, meta_schema.clone())
                .next()
                .unwrap();
            metadata
        };
        Ok(metadata)
    }
}
