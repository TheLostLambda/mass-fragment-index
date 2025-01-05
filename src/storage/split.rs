use std::{
    collections::{HashMap, HashSet},
    fs, io,
    marker::PhantomData,
    ops::Range,
    path::{Path, PathBuf},
    sync::Arc,
};

use arrow::{
    array::{
        ArrayRef, Float32Array, Float32Builder, RecordBatch, StringArray, StringBuilder,
        StructArray, UInt32Array, UInt32Builder,
    },
    datatypes::{DataType, Field, Schema, SchemaBuilder},
    json::{LineDelimitedWriter, ReaderBuilder as JSONReaderBuilder},
};
use itertools::{izip, Itertools};
use parquet::{
    arrow::{
        arrow_reader::{
            statistics::StatisticsConverter, ArrowReaderBuilder, ArrowReaderOptions,
            ParquetRecordBatchReader, ParquetRecordBatchReaderBuilder, RowSelector,
        },
        ArrowWriter,
    },
    basic::{Compression, ZstdLevel},
};

use crate::{
    sort::{IndexBin, ParentID},
    IndexSortable, Interval, MassType, Tolerance,
};

use super::{
    util::{afield, as_array_ref},
    ArrowStorage, IndexBinaryStorage, IndexMetadata,
};

#[allow(unused)]
pub trait SplitArrowStorage: ArrowStorage {
    fn split_schema() -> Arc<Schema> {
        let base_schema = Self::schema();
        let mut fields = SchemaBuilder::from(base_schema.fields.clone());
        fields.push(afield!("band_id", DataType::UInt32));
        *fields.metadata_mut() = base_schema.metadata.clone();
        Arc::new(fields.finish())
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct SplitBand {
    pub band_id: ParentID,
    pub start_id: ParentID,
    pub end_id: ParentID,
    pub start_mass: MassType,
    pub end_mass: MassType,
    pub file_name: Option<String>,
}

impl SplitBand {
    pub fn new(
        band_id: ParentID,
        start_id: ParentID,
        end_id: ParentID,
        start_mass: MassType,
        end_mass: MassType,
        file_name: Option<String>,
    ) -> Self {
        Self {
            band_id,
            start_id,
            end_id,
            start_mass,
            end_mass,
            file_name,
        }
    }

    pub fn as_id_interval(&self) -> Interval {
        Interval::new(self.start_id as usize, self.end_id as usize + 1)
    }
}

impl ArrowStorage for SplitBand {
    fn schema() -> arrow::datatypes::SchemaRef {
        let band_id = afield!("band_id", DataType::UInt32);
        let start_id = afield!("start_id", DataType::UInt32);
        let end_id = afield!("end_id", DataType::UInt32);
        let start_mass = afield!("start_mass", DataType::Float32);
        let end_mass = afield!("end_mass", DataType::Float32);
        let file_name = afield!("file_name", DataType::Utf8);
        Arc::new(Schema::new(vec![
            band_id, start_id, end_id, start_mass, end_mass, file_name,
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
        let file_names: &StringArray = batch.column(5).as_any().downcast_ref().unwrap();

        izip!(
            band_ids,
            start_ids,
            end_ids,
            start_masses,
            end_masses,
            file_names
        )
        .map(
            |(band_id, start_id, end_id, start_mass, end_mass, file_name)| {
                (
                    SplitBand::new(
                        band_id.unwrap(),
                        start_id.unwrap(),
                        end_id.unwrap(),
                        start_mass.unwrap(),
                        end_mass.unwrap(),
                        file_name.map(|s| s.to_string()),
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
        let mut file_names = StringBuilder::new();

        for item in batch {
            band_ids.append_value(item.band_id);
            start_ids.append_value(item.start_id);
            end_ids.append_value(item.end_id);
            start_masses.append_value(item.start_mass);
            end_masses.append_value(item.end_mass);
            file_names.append_option(item.file_name.clone());
        }

        let columns = vec![
            as_array_ref!(band_ids),
            as_array_ref!(start_ids),
            as_array_ref!(end_ids),
            as_array_ref!(start_masses),
            as_array_ref!(end_masses),
            as_array_ref!(file_names),
        ];

        let batch = RecordBatch::try_new(schema, columns);
        batch
    }

    fn archive_name() -> String {
        "search_bands.jsonl".into()
    }

    fn writer_properties() -> parquet::file::properties::WriterPropertiesBuilder {
        parquet::file::properties::WriterProperties::builder()
    }
}

#[derive(Debug, Default, Clone)]
pub enum BinStorageStrategy {
    #[default]
    SingleFile,
    FilePerBin,
    NBinsPerFile(usize),
    NEntriesPerFile(usize),
}

impl BinStorageStrategy {
    pub fn make_file_name<T: ArrowStorage>(&self, bands: &[SplitBand], i: usize) -> String {
        let base = T::archive_name();
        let delim_idx = base.rfind(".").unwrap_or(base.len());
        let (prefix, suffix) = base.split_at(delim_idx);
        match self {
            BinStorageStrategy::SingleFile => base,
            BinStorageStrategy::FilePerBin => format!("{prefix}_{i}.{suffix}"),
            BinStorageStrategy::NBinsPerFile(n) => {
                let (lo, hi) = bands
                    .iter()
                    .fold((ParentID::MAX, ParentID::MIN), |(min, max), band| {
                        (min.min(band.start_id), max.max(band.end_id))
                    });
                format!("{prefix}_{}_b{}_{lo}_{hi}.{suffix}", i, i * *n)
            }
            BinStorageStrategy::NEntriesPerFile(n) => {
                format!("{prefix}_{}_e{}.{suffix}", i, i * *n)
            }
        }
    }
}

#[derive(Debug)]
pub struct SplitStorageOptions {
    pub bin_width: MassType,
    pub bin_storage_strategy: BinStorageStrategy,
}

impl Default for SplitStorageOptions {
    fn default() -> Self {
        Self {
            bin_width: 50.0,
            bin_storage_strategy: Default::default(),
        }
    }
}

pub struct SplitIndexBinaryStorageWriter<
    'a,
    T: ArrowStorage + 'a + IndexSortable + Clone + SplitArrowStorage,
    P: ArrowStorage + IndexSortable,
    M: ArrowStorage,
    I: SplitIndexBinaryStorage<'a, T, P, M>,
> {
    config: SplitStorageOptions,
    index: &'a I,
    bands: Vec<SplitBand>,
    _t: PhantomData<T>,
    _p: PhantomData<P>,
    _m: PhantomData<M>,
}

impl<
        'a,
        T: ArrowStorage + 'a + IndexSortable + Clone + SplitArrowStorage,
        P: ArrowStorage + IndexSortable,
        M: ArrowStorage,
        I: SplitIndexBinaryStorage<'a, T, P, M>,
    > SplitIndexBinaryStorageWriter<'a, T, P, M, I>
{
    pub fn new(config: SplitStorageOptions, index: &'a I, bands: Vec<SplitBand>) -> Self {
        Self {
            config,
            index,
            bands,
            _t: PhantomData,
            _p: PhantomData,
            _m: PhantomData,
        }
    }

    pub fn write(&mut self, directory: &Path, compression_level: &Compression) -> io::Result<()> {
        match self.config.bin_storage_strategy {
            BinStorageStrategy::SingleFile => self.write_single_file(directory, compression_level),
            BinStorageStrategy::FilePerBin => self.write_bin_per_file(directory, compression_level),
            BinStorageStrategy::NBinsPerFile(_) => {
                self.write_n_bands_per_file(directory, compression_level)
            }
            BinStorageStrategy::NEntriesPerFile(_) => {
                self.write_n_entries_per_file(directory, compression_level)
            }
        }
    }

    pub fn write_single_file(
        &mut self,
        directory: &Path,
        compression_level: &Compression,
    ) -> io::Result<()> {
        self.write_bands(0..self.bands.len(), 0, directory, compression_level)?;
        Ok(())
    }

    pub fn write_bin_per_file(
        &mut self,
        directory: &Path,
        compression_level: &Compression,
    ) -> io::Result<()> {
        let idxes = 0..self.bands.len();
        for band_i in idxes {
            self.write_bands(band_i..(band_i + 1), band_i, directory, compression_level)?;
        }
        Ok(())
    }

    pub fn write_n_bands_per_file(
        &mut self,
        directory: &Path,
        compression_level: &Compression,
    ) -> io::Result<()> {
        if let BinStorageStrategy::NBinsPerFile(n_per_file) =
            self.config.bin_storage_strategy.clone()
        {
            let mut start = 0usize;
            let n = self.bands.len();
            let mut bands_i = 0;
            while start < n {
                let indices = start..(start + n_per_file).min(n);
                self.write_bands(indices, bands_i, directory, compression_level)?;
                bands_i += 1;
                start += n_per_file;
            }
        }
        Ok(())
    }

    pub fn write_n_entries_per_file(
        &mut self,
        directory: &Path,
        compression_level: &Compression,
    ) -> io::Result<()> {
        if let BinStorageStrategy::NEntriesPerFile(n_per_file) =
            self.config.bin_storage_strategy.clone()
        {
            let mut band_start = 0usize;
            let n_bands = self.bands.len();
            let mut bands_i = 0;
            while band_start < n_bands {
                let mut n_entries = 0;
                let mut end_band = n_bands;
                for (band_i, band) in self.bands.iter().enumerate().skip(band_start) {
                    let interval = Interval::new(band.start_id as usize, band.end_id as usize + 1);
                    for (_, bin) in self.index.iter_entries().enumerate() {
                        n_entries += bin
                            .iter()
                            .filter(|b| interval.contains(b.parent_id() as usize))
                            .count();
                        if n_entries > n_per_file {
                            end_band = band_i;
                            break;
                        }
                    }
                }
                self.write_bands(
                    band_start..(end_band + 1).min(n_bands),
                    bands_i,
                    directory,
                    compression_level,
                )?;
                bands_i += 1;
                band_start = end_band;
            }
        }

        Ok(())
    }

    pub fn write_bands(
        &mut self,
        band_indices: Range<usize>,
        bands_i: usize,
        directory: &Path,
        compression_level: &Compression,
    ) -> io::Result<()> {
        let archive_name = self
            .config
            .bin_storage_strategy
            .make_file_name::<T>(&self.bands, bands_i);

        let entries_path = directory.join(archive_name.clone());
        let entries_schema = T::schema();
        let props = T::writer_properties()
            .set_compression(compression_level.clone())
            .set_column_encoding("band_id".into(), parquet::basic::Encoding::RLE)
            .set_writer_version(parquet::file::properties::WriterVersion::PARQUET_2_0)
            .set_statistics_enabled(parquet::file::properties::EnabledStatistics::Page)
            .build();

        log::debug!(
            "Opening {archive_name} for bands {}..{}",
            band_indices.start,
            band_indices.end
        );
        let ext_schema = I::make_item_schema();
        let mut writer = ArrowWriter::try_new(
            fs::File::create(entries_path)?,
            ext_schema.clone(),
            Some(props),
        )?;

        for j in band_indices {
            let band = self.bands.get_mut(j).unwrap();
            log::debug!(
                "Writing band {j}: {:0.2}-{:0.2}",
                band.start_mass,
                band.end_mass
            );
            let interval = Interval::new(band.start_id as usize, band.end_id as usize + 1);
            band.file_name = Some(archive_name.clone());
            for (i, bin) in self.index.iter_entries().enumerate() {
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
}

pub trait SplitIndexBinaryStorage<
    'a,
    T: ArrowStorage + 'a + IndexSortable + Clone + SplitArrowStorage,
    P: ArrowStorage + IndexSortable,
    M: ArrowStorage,
>: IndexBinaryStorage<'a, T, P, M> + Sized
{
    fn make_item_schema() -> Arc<Schema> {
        T::split_schema()
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
                    None,
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
            None,
        ));
        bands
    }

    fn write_entries_split(
        &'a self,
        directory: &Path,
        mut bands: Vec<SplitBand>,
        compression_level: &Compression,
    ) -> io::Result<Vec<SplitBand>> {
        let opts = SplitStorageOptions::default();
        let mut writer = SplitIndexBinaryStorageWriter::new(opts, self, bands);
        writer.write(directory, compression_level)?;
        bands = writer.bands;
        Ok(bands)
    }

    fn write_split<D: AsRef<Path>>(
        &'a self,
        directory: &D,
        bin_options: SplitStorageOptions,
        compression_level: Option<Compression>,
    ) -> io::Result<()> {
        let directory = directory.as_ref();

        let compression_level =
            compression_level.unwrap_or_else(|| Compression::ZSTD(ZstdLevel::try_new(9).unwrap()));
        let mut bands = self.compute_parent_bands(bin_options.bin_width);

        self.write_metadata(directory)?;
        self.write_parents(directory, &compression_level)?;
        bands = self.write_entries_split(directory, bands, &compression_level)?;
        self.write_split_log(directory, &bands)?;
        Ok(())
    }

    fn band_log_name() -> String {
        SplitBand::archive_name()
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
        let split_log = Self::read_split_log(root)?;

        let file_names: HashSet<&str> = split_log
            .iter()
            .flat_map(|s| s.file_name.as_deref())
            .collect();
        let mut file_names: Vec<_> = file_names.into_iter().collect();
        file_names.sort();

        let entries = {
            let mut bin_collector: HashMap<u64, Vec<T>> = HashMap::default();
            for archive_name in file_names.iter() {
                let entries_fh = fs::File::open(root.join(archive_name))?;
                let reader = ArrowReaderBuilder::try_new(entries_fh)?.build()?;
                let entry_schema = T::schema();

                for batch in reader {
                    for (entry, segment_id) in T::from_batch(&batch.unwrap(), entry_schema.clone())
                    {
                        bin_collector.entry(segment_id).or_default().push(entry);
                    }
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

#[derive(Debug)]
pub struct SearchIndexOnDisk<
    T: ArrowStorage + IndexSortable + Default + SplitArrowStorage,
    P: ArrowStorage + IndexSortable + Default,
    M: ArrowStorage + Default = IndexMetadata,
> {
    root: PathBuf,
    pub metadata: M,
    parents: IndexBin<P>,
    split_bands: Vec<SplitBand>,
    _t: PhantomData<T>,
    _p: PhantomData<P>,
}

impl<
        T: ArrowStorage + IndexSortable + Default + SplitArrowStorage,
        P: ArrowStorage + IndexSortable + Default,
        M: ArrowStorage + Default,
    > SearchIndexOnDisk<T, P, M>
{
    pub fn new(path: PathBuf) -> io::Result<Self> {
        if !path.exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("Index root {} not found", path.display()),
            ));
        }
        if !path.join(M::archive_name()).exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!(
                    "Index metadata {} not found",
                    path.join(M::archive_name()).display()
                ),
            ));
        }
        if !path.join(P::archive_name()).exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!(
                    "Index parent file {} not found",
                    path.join(P::archive_name()).display()
                ),
            ));
        }
        if !path.join(T::archive_name()).exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!(
                    "Index search target file {} not found",
                    path.join(T::archive_name()).display()
                ),
            ));
        }
        if !path.join(SplitBand::archive_name()).exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!(
                    "Index split bands file {} not found",
                    path.join(SplitBand::archive_name()).display()
                ),
            ));
        }
        let mut this = Self {
            root: path,
            metadata: M::default(),
            split_bands: Vec::new(),
            parents: IndexBin::default(),
            _t: PhantomData,
            _p: PhantomData,
        };
        this.metadata = this.read_metadata()?;
        this.split_bands = this.read_split_bands()?;
        this.parents = this.read_parents()?;
        Ok(this)
    }

    fn read_parents(&self) -> io::Result<IndexBin<P>> {
        let parents_path = self.root.join(P::archive_name());
        let parent_schema = P::schema();
        let parents_fh = fs::File::open(parents_path)?;

        let reader = ArrowReaderBuilder::try_new(parents_fh)?.build()?;
        let mut parents = Vec::new();
        for batch in reader {
            parents.extend(P::from_batch(&batch.unwrap(), parent_schema.clone()).map(|(p, _)| p));
        }
        let parents = IndexBin::from(parents);
        Ok(parents)
    }

    fn read_split_bands(&self) -> io::Result<Vec<SplitBand>> {
        let arch = self.root.join(SplitBand::archive_name());
        let handle = io::BufReader::new(fs::File::open(arch)?);
        let reader = JSONReaderBuilder::new(SplitBand::schema())
            .build(handle)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        let mut buffer = Vec::new();
        let schema = SplitBand::schema();
        for batch in reader.flatten() {
            buffer.extend(SplitBand::from_batch(&batch, schema.clone()).map(|(b, _)| b));
        }
        Ok(buffer)
    }

    fn read_metadata(&self) -> io::Result<M> {
        let arch = self.root.join(M::archive_name());
        let handle = io::BufReader::new(fs::File::open(arch)?);
        let mut reader = JSONReaderBuilder::new(M::schema())
            .build(handle)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        let (meta, _) = M::from_batch(
            &reader
                .next()
                .expect("No metadata record batch found")
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?,
            M::schema(),
        )
        .next()
        .unwrap();
        Ok(meta)
    }

    pub fn parents_for(&self, mass: MassType, error_tolerance: Tolerance) -> Interval {
        let iv = self.parents.search_mass(mass, error_tolerance);
        iv
    }

    pub fn parents_for_range(
        &self,
        low: MassType,
        high: MassType,
        error_tolerance: Tolerance,
    ) -> Interval {
        let mut out = Interval::default();
        out.start = self.parents_for(low, error_tolerance).start;
        out.end = self.parents_for(high, error_tolerance).end;
        out
    }

    fn build_parquet_reader(
        &self,
        handle: fs::File,
        interval: &Interval,
        band_ids: &HashSet<u32>,
    ) -> io::Result<ParquetRecordBatchReader> {
        let builder = ParquetRecordBatchReaderBuilder::try_new_with_options(
            handle,
            ArrowReaderOptions::new()
                .with_page_index(true)
                .with_schema(T::split_schema()),
        )?;

        let metadata = builder.metadata();
        let column_indices = metadata.column_index().unwrap();
        let column_offsets = metadata.offset_index().unwrap();
        let row_groups = metadata.row_groups();

        let schema = builder.schema();
        let parquet_schema = builder.parquet_schema();

        let i = T::parent_id_column().unwrap();
        let sort_col = schema.field(i);
        let (_, band_col) = schema.column_with_name("band_id").unwrap();

        let parent_id_statistics =
            StatisticsConverter::try_new(sort_col.name(), &schema, parquet_schema)?;

        let band_id_statistics =
            StatisticsConverter::try_new(band_col.name(), &schema, parquet_schema)?;

        let row_group_indices: Vec<_> = (0..row_groups.len()).collect();

        let parent_id_maxes = parent_id_statistics.data_page_maxes(
            column_indices,
            column_offsets,
            &row_group_indices,
        )?;
        let parent_id_mins = parent_id_statistics.data_page_mins(
            column_indices,
            column_offsets,
            &row_group_indices,
        )?;

        let band_id_maxes = band_id_statistics.data_page_maxes(
            column_indices,
            column_offsets,
            &row_group_indices,
        )?;
        let band_id_mins = band_id_statistics.data_page_mins(
            column_indices,
            column_offsets,
            &row_group_indices,
        )?;

        let parent_id_maxes: &UInt32Array = parent_id_maxes.as_any().downcast_ref().unwrap();
        let parent_id_mins: &UInt32Array = parent_id_mins.as_any().downcast_ref().unwrap();

        let band_id_maxes: &UInt32Array = band_id_maxes.as_any().downcast_ref().unwrap();
        let band_id_mins: &UInt32Array = band_id_mins.as_any().downcast_ref().unwrap();

        let counts = parent_id_statistics
            .data_page_row_counts(column_offsets, row_groups, &row_group_indices)?
            .unwrap();

        let mut selectors = Vec::new();
        for (parent_id_min, parent_id_max, count, band_id_min, band_id_max) in
            itertools::multizip((
                parent_id_mins.iter(),
                parent_id_maxes.iter(),
                counts.iter(),
                band_id_mins.iter(),
                band_id_maxes.iter(),
            ))
        {
            let parent_id_min = parent_id_min.unwrap();
            let parent_id_max = parent_id_max.unwrap();
            let band_id_min = band_id_min.unwrap();
            let band_id_max = band_id_max.unwrap();

            let band_overlaps = (band_id_min..=band_id_max).any(|i| band_ids.contains(&i));
            let parent_overlaps = Interval::new(parent_id_min as usize, parent_id_max as usize + 1)
                .overlaps(&interval);

            let count = count.unwrap();
            if parent_overlaps && band_overlaps {
                selectors.push(RowSelector::select(count as usize));
            } else {
                selectors.push(RowSelector::skip(count as usize));
            }
        }
        Ok(builder.with_row_selection(selectors.into()).build()?)
    }

    pub fn load_for_parents(&self, interval: Interval) -> io::Result<Vec<T>> {
        let bands: Vec<_> = self
            .split_bands
            .iter()
            .filter(|band| interval.overlaps(&band.as_id_interval()))
            .collect();

        let band_ids: HashSet<_> = bands.iter().map(|band| band.band_id).collect();

        let file_names = bands
            .iter()
            .flat_map(|band| band.file_name.as_deref())
            .unique();

        let schema = T::split_schema();

        let mut items = Vec::new();

        for file_name in file_names {
            let archive_name = self.root.join(file_name);
            let handle = fs::File::open(archive_name)?;

            let reader = self.build_parquet_reader(handle, &interval, &band_ids)?;

            for batch in reader {
                let batch = batch.unwrap();
                items.extend(
                    T::from_batch(&batch, schema.clone())
                        .filter(|(x, _)| interval.contains(x.parent_id() as usize))
                        .map(|(x, _)| x),
                );
            }
        }
        Ok(items)
    }

    pub fn root(&self) -> &PathBuf {
        &self.root
    }

    pub fn split_bands(&self) -> &[SplitBand] {
        &self.split_bands
    }

    pub fn parents(&self) -> &IndexBin<P> {
        &self.parents
    }

    pub fn read_groups(&self, interval: Interval) -> io::Result<HashMap<u32, Vec<T>>> {
        let entries = self.load_for_parents(interval)?;

        let mut index: HashMap<u32, Vec<_>> = HashMap::new();
        for e in entries.into_iter() {
            index.entry(e.parent_id()).or_default().push(e);
        }

        Ok(index)
    }
}
