use std::{env, io};

use mzdata::{io::DetailLevel, prelude::*, MzMLReader};

use mass_fragment_index::{
    sort::SortType, storage::{read_peak_index, write_peak_index}, DeconvolutedPeak, DeconvolutedSpectrumIndex, Spectrum
};

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let mzml_path = args
        .next()
        .unwrap_or_else(|| panic!("Please provide a path to an mzML file"));

    let storage_dir = args
        .next()
        .unwrap_or_else(|| panic!("Please provide a path to write index to"));

    let mut reader = MzMLReader::open_path(mzml_path)?;
    reader.detail_level = DetailLevel::MetadataOnly;
    eprintln!("Loading spectra");
    let mut spectra: Vec<_> = reader
        .iter()
        .filter(|s| s.ms_level() > 1)
        .map(|s| {
            let prec = s.precursor().unwrap();
            Spectrum::new(
                prec.neutral_mass() as f32,
                prec.charge().unwrap(),
                0,
                s.index() as u32,
                u32::MAX,
            )
        })
        .collect();
    eprintln!("Sorting spectra");
    spectra.sort_by(|a, b| {
        a.precursor_mass.total_cmp(&b.precursor_mass).then_with(|| {
            a.precursor_charge.cmp(&b.precursor_charge).then_with(|| {
                a.scan_number
                    .cmp(&b.scan_number)
                    .then_with(|| a.source_file_id.cmp(&b.source_file_id))
            })
        })
    });

    spectra
        .iter_mut()
        .enumerate()
        .for_each(|(i, s)| s.sort_id = i as u32);

    let mut index = DeconvolutedSpectrumIndex::empty(100, 3000.0);

    reader.detail_level = DetailLevel::Full;
    eprintln!("Loading peaks");
    spectra.into_iter().for_each(|s| {
        let mut source = reader
            .get_spectrum_by_index(s.scan_number as usize)
            .unwrap();
        let scan_ref = s.sort_id;
        source.try_build_deconvoluted_centroids().unwrap();
        index.add_parent(s);
        for p in source.deconvoluted_peaks.unwrap().iter() {
            index.add(DeconvolutedPeak::new(
                p.neutral_mass as f32,
                p.charge as i16,
                p.intensity,
                scan_ref,
            ));
        }
    });

    eprintln!("Sorting index");
    index.sort(SortType::ByParentId);

    eprintln!("Writing index");
    write_peak_index(&index, &storage_dir)?;

    let iv = index.parents_for_range(1200.0, 2300.0, mass_fragment_index::Tolerance::Da(0.2));

    let it = index.search(1200.0, mass_fragment_index::Tolerance::Da(0.2), Some(iv));
    eprintln!("Parent Interval: {iv:?}");
    let searched_peaks: Vec<_> = it.collect();
    eprintln!("Found peaks: {}", searched_peaks.len());

    eprintln!("Loading index from disk: {}", storage_dir);
    let duplicate = read_peak_index(&storage_dir)?;
    assert_eq!(duplicate.parents.entries, index.parents.entries);

    for (i, (abin, bbin)) in index.bins.iter().zip(duplicate.bins.iter()).enumerate() {
        assert_eq!(abin.len(), bbin.len(), "Bin {i} not equal: {} != {}", abin.len(), bbin.len());
    }

    let iv = duplicate.parents_for_range(1200.0, 2300.0, mass_fragment_index::Tolerance::Da(0.2));

    let it = duplicate.search(1200.0, mass_fragment_index::Tolerance::Da(0.2), Some(iv));
    eprintln!("Parent Interval: {iv:?}");
    let searched_peaks: Vec<_> = it.collect();
    eprintln!("Found peaks: {}", searched_peaks.len());

    Ok(())
}
