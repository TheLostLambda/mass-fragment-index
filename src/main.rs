use std::{io, fs};

use csv;

use fragment_index_rs::index::SearchIndex;
use fragment_index_rs::fragment::{Fragment, FragmentName};
use fragment_index_rs::parent::{Peptide};

fn main() -> io::Result<()> {
    let reader = io::BufReader::new(fs::File::open("tests/data/test_data.csv")?);
    let mut csv_reader = csv::Reader::from_reader(reader);

    let mut search_index: SearchIndex<Fragment, Peptide> = SearchIndex::empty(100, 3000.0);
    let mut parent_i = -1;
    for result in csv_reader.records() {
        let record = result?;
        match record.get(0).expect("Missing record type") {
            "PEPTIDE" => {
                parent_i += 1;
                let peprec = Peptide::new(
                    record.get(2).unwrap().parse::<f32>().unwrap(),
                    parent_i as usize,
                    0,
                    0,
                    record.get(1).unwrap().to_string()
                );
                search_index.add_parent(peprec);
            },
            "FRAGMENT" => {
                let name: FragmentName = record.get(1).unwrap().parse().unwrap();
                let frag = Fragment::new(
                    record.get(2).unwrap().parse::<f32>().unwrap(),
                    parent_i as usize,
                    name.0,
                    name.1,
                );
                search_index.add(frag);
            },
            field => {
                panic!("Unknown record type {}", field)
            }
        }
    }
    search_index.sort(fragment_index_rs::sort::SortType::ByParentId);
    println!("Num Bins: {}", search_index.bins.len());
    let parent_interval = search_index.parents_for_range(1158.0, 1163.0, 1e-5);
    println!("Parent Interval: {}-{}", parent_interval.start, parent_interval.end);
    for i in parent_interval.start..parent_interval.end {
        let parent = &search_index.parents.entries[i];
        println!("Peptide: {}:{:0.3} @ {}", parent.sequence, parent.mass, parent.id);
    }
    let iter = search_index.search(174.111, 4e-5, Some(parent_interval));

    println!("Fragments:");
    for frag in iter {
        println!("\tFragment {:?}", frag);
    }
    Ok(())
}
