#[cfg(test)]

use std::{mem, io, fs};

use csv;
use rand::Rng;

use fragment_index_rs::index::SearchIndex;
use fragment_index_rs::fragment::{Fragment, FragmentName};
use fragment_index_rs::parent::Peptide;


fn parse_csv<R: io::BufRead>(reader: R) -> io::Result<Vec<(Peptide, Vec<Fragment>)>> {
    let mut csv_reader = csv::Reader::from_reader(reader);
    let mut parent_i: i32 = -1;

    let mut accumulator = Vec::new();
    let mut peprec: Option<Peptide> = None;
    let mut fragments = Vec::new();

    for result in csv_reader.records() {
        let record = result?;
        match record.get(0).expect("Missing record type") {
            "PEPTIDE" => {
                parent_i += 1;
                if let Some(ref peptide) = peprec {
                    accumulator.push((peptide.clone(), mem::take(&mut fragments)));
                }

                peprec = Some(Peptide::new(
                    record.get(2).unwrap().parse::<f32>().unwrap(),
                    parent_i as usize,
                    0,
                    0,
                    record.get(1).unwrap().to_string(),
                ));
            }
            "FRAGMENT" => {
                let name: FragmentName = record.get(1).unwrap().parse().unwrap();
                let frag = Fragment::new(
                    record.get(2).unwrap().parse::<f32>().unwrap(),
                    parent_i as usize,
                    name.0,
                    name.1,
                );
                fragments.push(frag);
            }
            field => {
                panic!("Unknown record type {}", field)
            }
        }
    }
    if let Some(peptide) = peprec {
        accumulator.push((peptide, mem::take(&mut fragments)));
    }
    Ok(accumulator)
}

fn build_index<R: io::BufRead>(reader: R) -> io::Result<SearchIndex<Fragment, Peptide>> {
    let pepfrags = parse_csv(reader)?;
    let mut search_index: SearchIndex<Fragment, Peptide> = SearchIndex::empty(100, 10000.0);
    pepfrags.into_iter().for_each(|(pep, frags)| {
        search_index.add_parent(pep);
        frags.into_iter().for_each(|frag| {
            search_index.add(frag);
        });
    });

    search_index.sort(fragment_index_rs::sort::SortType::ByParentId);
    Ok(search_index)
}


fn test_exact(search_index: &SearchIndex<Fragment, Peptide>, entries: &Vec<(Peptide, Vec<Fragment>)>, precursor_error_tolerance: f32, product_error_tolerance: f32) -> bool {
    for (batch_i, (pept, frags)) in entries.iter().enumerate() {
        let parent_interval = search_index.parents_for(pept.mass, precursor_error_tolerance);
        for i in parent_interval {
            assert!(parent_interval.contains(i));
            let parent = &search_index.parents[i];
            assert!(
                ((parent.mass - pept.mass) / pept.mass).abs() <= precursor_error_tolerance,
                "Parent molecule {}/{:?} {:?} did not pass precursor error tolerance {} for {:?}",
                i,
                parent_interval,
                parent,
                ((parent.mass - pept.mass) / pept.mass).abs(),
                pept
            );
        }

        let n_expected_parents = entries
            .iter()
            .filter(|(alt_pept, _)| {
                ((alt_pept.mass - pept.mass).abs() / pept.mass) < precursor_error_tolerance
            })
            .count();

        let n_parents = parent_interval.clone().into_iter().count();

        assert_eq!(n_expected_parents, n_parents);

        for expected_frag in frags.iter() {
            let search_iter = search_index.search(
                expected_frag.mass,
                product_error_tolerance,
                Some(parent_interval),
            );

            let n_frags =
                search_iter
                    .enumerate()
                    .map(|(i, frag)| {
                        assert!(
                    parent_interval.contains(frag.parent_id),
                    "{}th hit {:?} did not fall within parent interval {:?} {:?} of batch {}",
                    i, frag, parent_interval, expected_frag, batch_i
                );
                        assert!(
                            (frag.mass - expected_frag.mass).abs() / frag.mass
                                < product_error_tolerance
                        );
                        frag
                    })
                    .count();

            let n_expected_frags = entries
                .iter()
                .filter(|(alt_pept, _)| {
                    ((alt_pept.mass - pept.mass).abs() / pept.mass) < precursor_error_tolerance
                })
                .map(|(_, frags)| {
                    frags.iter().filter(|frag| {
                        (frag.mass - expected_frag.mass).abs() / frag.mass < product_error_tolerance
                    })
                })
                .flatten()
                .count();
            assert_eq!(n_expected_frags, n_frags);
        }
    }
    true
}


fn test_permuted(search_index: &SearchIndex<Fragment, Peptide>, entries: &Vec<(Peptide, Vec<Fragment>)>, precursor_error_tolerance: f32, product_error_tolerance: f32) -> bool {
    let mut rng = rand::thread_rng();
    for (batch_i, (pept, frags)) in entries.iter().enumerate() {
        let mass_error_scale = (pept.mass * precursor_error_tolerance) * 0.99;

        // When the random offset is very close to the outer limit of the tolerance window, things
        // seem to fall apart. Try using a larger dataset?
        let precursor_err = rng.gen_range(-mass_error_scale..mass_error_scale);

        let parent_interval = search_index.parents_for(pept.mass + precursor_err, precursor_error_tolerance);
        for i in parent_interval {
            assert!(parent_interval.contains(i));
            let parent = &search_index.parents[i];

            assert!(
                ((parent.mass - (pept.mass + precursor_err)) / (pept.mass + precursor_err)).abs() <= precursor_error_tolerance,
                "Parent molecule {}/{:?} {:?} did not pass precursor error tolerance {} for {:?}  with random offset {} ({}% of {})",
                i,
                parent_interval,
                parent,
                ((parent.mass - (pept.mass + precursor_err)) / (pept.mass  + precursor_err)).abs(),
                pept,
                precursor_err,
                precursor_err * 100.0 / mass_error_scale,
                mass_error_scale
            );
        }

        let n_expected_parents = entries
            .iter()
            .filter(|(alt_pept, _)| {
                ((alt_pept.mass - (pept.mass + precursor_err)).abs() / (pept.mass + precursor_err)) < precursor_error_tolerance
            })
            .count();

        let n_parents = parent_interval.clone().into_iter().count();

        assert_eq!(
            n_expected_parents, n_parents,
            "Expected parent count differed for {} with random offset {} ({}% of {})",
            pept.mass,
            precursor_err,
            precursor_err * 100.0 / mass_error_scale,
            mass_error_scale
        );

        for expected_frag in frags.iter() {
            let search_iter = search_index.search(
                expected_frag.mass,
                product_error_tolerance,
                Some(parent_interval),
            );

            let n_frags =
                search_iter
                    .enumerate()
                    .map(|(i, frag)| {
                        assert!(
                    parent_interval.contains(frag.parent_id),
                    "{}th hit {:?} did not fall within parent interval {:?} {:?} of batch {}",
                    i, frag, parent_interval, expected_frag, batch_i
                );
                        assert!(
                            (frag.mass - expected_frag.mass).abs() / frag.mass
                                < product_error_tolerance
                        );
                        frag
                    })
                    .count();

            let n_expected_frags = entries
                .iter()
                .filter(|(alt_pept, _)| {
                    ((alt_pept.mass - (pept.mass + precursor_err)).abs() / (pept.mass + precursor_err)) < precursor_error_tolerance
                })
                .map(|(_, frags)| {
                    frags.iter().filter(|frag| {
                        (frag.mass - expected_frag.mass).abs() / frag.mass < product_error_tolerance
                    })
                })
                .flatten()
                .count();
            assert!(n_expected_frags == n_frags);
        }
    }
    true
}


#[test]
fn test_index_build_traversal() -> io::Result<()> {
    let reader = io::BufReader::new(fs::File::open("tests/data/test_data.csv")?);
    let search_index = build_index(reader)?;
    let reader = io::BufReader::new(fs::File::open("tests/data/test_data.csv")?);
    let entries = parse_csv(reader)?;

    assert!(1013803 == search_index.bins.len());
    assert!(search_index.parents.len() == entries.len());
    assert_eq!(
        search_index.bins.iter().map(|b| b.len()).sum::<usize>(),
        entries.iter().map(|(_, frags)| frags.len()).sum::<usize>()
    );

    let precursor_error_tolerance = 5e-6;
    let product_error_tolerance = 1e-5;

    assert!(test_exact(&search_index, &entries, precursor_error_tolerance, product_error_tolerance));

    for _ in 0..10 {
        assert!(test_permuted(&search_index, &entries, precursor_error_tolerance, product_error_tolerance));
    }

    let parent_interval = search_index.parents_for_range(1158.0, 1163.0, 1e-5);
    assert!(parent_interval.start == 49);
    assert!(parent_interval.end == 51);
    let iter = search_index.search(
        174.111, 2e-5, Some(parent_interval));
    let matched_fragments: Vec<Fragment> = iter.collect();
    assert!(matched_fragments.len() == 2);
    assert!(matched_fragments.iter().map(|f| ((f.mass - 174.111).abs() / f.mass <= 2e-5) as u8).sum::<u8>() == 2);
    Ok(())
}