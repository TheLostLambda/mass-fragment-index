use std::{
    fs, io::{self, prelude::*}, path::{Path, PathBuf}, sync::Arc
};

use arrow::{array::{ArrayRef, AsArray, RecordBatch, StringBuilder, UInt64Builder}, datatypes::{DataType, Field, Schema, UInt64Type}, json::{LineDelimitedWriter, ReaderBuilder as JSONReaderBuilder}};
use parquet::{arrow::{arrow_reader::ArrowReaderBuilder, ArrowWriter}, basic::{Compression, ZstdLevel}};

use super::{ArrowStorage, IndexBinaryStorage, util::{afield, as_array_ref}};

pub fn split_file_log_schema() -> Arc<Schema> {
    let bin_id = afield!("bin_id", DataType::UInt64);
    let filename = afield!("filename", DataType::Utf8);
    Arc::new(Schema::new(vec![bin_id, filename]))
}

pub fn split_file_log_to_arrow(split_log: &[(u64, String)]) -> RecordBatch {
    let mut bin_id_builder = UInt64Builder::new();
    let mut filename_builder = StringBuilder::new();

    for (i, fname) in split_log {
        bin_id_builder.append_value(*i);
        filename_builder.append_value(fname.to_string());
    }

    let columns = vec![
        as_array_ref!(bin_id_builder),
        as_array_ref!(filename_builder),
    ];

    let batch = RecordBatch::try_new(split_file_log_schema(), columns);
    batch.unwrap()
}


pub fn split_file_log_from_arrow(batch: RecordBatch) -> Vec<(u64, String)> {
    let bin_ids = batch.column_by_name("bin_id").unwrap().as_primitive::<UInt64Type>();
    let filenames = batch.column_by_name("filename").unwrap().as_string::<i32>();

    let mut split_log = Vec::new();

    for (bin_id, fname) in bin_ids.into_iter().zip(filenames.into_iter()) {
        split_log.push((bin_id.unwrap(), fname.unwrap().to_string()))
    }

    split_log
}


pub trait SplitArrowStorage : ArrowStorage {

    fn split_archive_name_prefix() -> String {
        let arch_base = Self::archive_name();
        if let Some((prefix, _suffix)) = arch_base.rsplit_once('.') {
            prefix.into()
        } else {
            arch_base
        }
    }

    fn split_archive_name_for(segment: u64) -> String {
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

pub trait SplitIndexBinaryStorage<'a, T: ArrowStorage + 'a + SplitArrowStorage, P: ArrowStorage, M: ArrowStorage> : IndexBinaryStorage<'a, T, P, M> {
    fn write_split<D: AsRef<Path>>(
        &'a self,
        directory: &D,
        compression_level: Option<Compression>,
    ) -> io::Result<()> {
        let directory = directory.as_ref();

        let compression_level =
            compression_level.unwrap_or_else(|| Compression::ZSTD(ZstdLevel::try_new(9).unwrap()));

        self.write_meta(directory)?;

        {
            let parent_path = directory.join(P::archive_name());
            let parent_schema = P::schema();
            let props = P::writer_properties()
                .set_compression(compression_level.clone())
                .build();
            let mut writer = ArrowWriter::try_new(
                fs::File::create(parent_path)?,
                parent_schema.clone(),
                Some(props),
            )?;
            let batch = P::to_batch(self.parents(), parent_schema.clone(), 0).unwrap();
            writer.write(&batch)?;
            writer.close()?;
        }

        {
            let entries_path = directory.join(T::archive_name());
            let entries_schema = T::schema();
            let props = T::writer_properties()
                .set_compression(compression_level.clone())
                .build();
            let mut writer = ArrowWriter::try_new(
                fs::File::create(entries_path)?,
                entries_schema.clone(),
                Some(props),
            )?;
            for (i, bin) in self.iter_entries().enumerate() {
                let batch = T::to_batch(bin, entries_schema.clone(), i as u64).unwrap();
                writer.write(&batch)?;
            }
            writer.close()?;
        }

        Ok(())
    }

    fn write_meta(&self, directory: &Path) -> io::Result<()> {
        let metadata = self.to_metadata();
        let meta_path = directory.join(M::archive_name());
        let meta_schema = M::schema();
        let meta_fh = io::BufWriter::new(fs::File::create(meta_path)?);
        let mut writer = LineDelimitedWriter::new(meta_fh);
        let metadata = M::to_batch(&[metadata], meta_schema, 0).unwrap();

        writer.write(&metadata).unwrap();
        writer.finish().unwrap();
        Ok(())
    }

    fn write_parents(&self, directory: &Path, compression_level: &Compression) -> io::Result<()> {
        let parent_path = directory.join(P::archive_name());
        let parent_schema = P::schema();
        let props = P::writer_properties()
            .set_compression(compression_level.clone())
            .build();
        let mut writer = ArrowWriter::try_new(
            fs::File::create(parent_path)?,
            parent_schema.clone(),
            Some(props),
        )?;
        let batch = P::to_batch(self.parents(), parent_schema.clone(), 0).unwrap();
        writer.write(&batch)?;
        writer.close()?;
        Ok(())
    }

    fn write_entries(&'a self, directory: &Path, compression_level: &Compression) -> io::Result<Vec<(u64, String)>> {
        let entries_schema = T::schema();
        let props = T::writer_properties()
            .set_compression(compression_level.clone())
            .build();
        let mut names = Vec::new();
        for (i, bin) in self.iter_entries().enumerate() {
            let batch = T::to_batch(bin, entries_schema.clone(), i as u64).unwrap();
            let fname = T::split_archive_name_for(i as u64);
            names.push((i as u64, fname.clone()));
            let path_of = directory.join(fname);
            let mut writer = ArrowWriter::try_new(
                fs::File::create(path_of)?,
                entries_schema.clone(),
                Some(props.clone()),
            )?;
            writer.write(&batch)?;
            writer.close()?;
        }

        Ok(names)
    }

    fn split_log_name() -> String {
        let mut prefix = T::split_archive_name_prefix();
        prefix.extend("_splits.json".chars());
        prefix
    }

    fn write_split_log(&self, directory: &Path, split_log: &[(u64, String)]) -> io::Result<()> {
        let split_log_path = directory.join(Self::split_log_name());

        let split_log_fh = io::BufWriter::new(fs::File::create(split_log_path)?);
        let mut writer = LineDelimitedWriter::new(split_log_fh);
        let split_log = split_file_log_to_arrow(split_log);

        writer.write(&split_log).unwrap();
        writer.finish().unwrap();
        Ok(())
    }

    fn read_split<D: AsRef<Path>>(directory: &D) -> io::Result<()> {
        let root = directory.as_ref();
        let meta = Self::read_metadata(root)?;
        let parents = Self::read_parents(root)?;

        Ok(())
    }

    fn read_parents(directory: &Path) -> io::Result<Vec<P>> {
        let parents_path = directory.join(P::archive_name());
        let parent_schema = P::schema();
        let parents_fh = fs::File::open(parents_path)?;

        let reader = ArrowReaderBuilder::try_new(parents_fh)?.build()?;
        let mut parents = Vec::new();
        for batch in reader {
            parents
                .extend(P::from_batch(&batch.unwrap(), parent_schema.clone()).map(|(p, _)| p));
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