
use std::ops::Index;

use crate::interval::Interval;


#[derive(Debug, Clone, Copy)]
pub enum SortType {
    ByMass,
    ByParentId,
    Unsorted,
}


impl Default for SortType {
    fn default() -> Self {
        Self::Unsorted
    }
}

pub trait IndexSortable {
    fn mass(&self) -> f32;
    fn parent_id(&self) -> usize;
}



#[derive(Debug, Clone, Default)]
pub struct IndexBin<T: IndexSortable> {
    pub entries: Vec<T>,
    pub sort_type: SortType,
    pub min_mass: f32,
    pub max_mass: f32,
}

impl<T: IndexSortable> IndexBin<T> {
    pub fn new(
        entries: Vec<T>,
        sort_type: SortType,
        min_mass: f32,
        max_mass: f32,
    ) -> Self {
        Self {
            entries,
            sort_type,
            min_mass,
            max_mass,
        }
    }

    pub fn append(&mut self, entry: T) {
        self.entries.push(entry);
        self.sort_type = SortType::Unsorted;
    }

    fn find_min_max_masses(&self) -> (f32, f32) {
        let mut min_mass = f32::INFINITY;
        let mut max_mass = 0f32;

        for f in self.entries.iter() {
            if f.mass() < min_mass {
                min_mass = f.mass();
            }
            if f.mass() > max_mass {
                max_mass = f.mass();
            }
        }
        return (min_mass, max_mass)
    }

    pub fn sort(&mut self, ordering: SortType) {
        match ordering {
            SortType::ByMass => {
                self.entries.sort_by(|a,b| {
                    a.mass().partial_cmp(&b.mass()).unwrap()
                });
                if let Some(f) = self.entries.first() {
                    self.min_mass = f.mass()
                }
                if let Some(f) = self.entries.last() {
                    self.max_mass = f.mass()
                }
            },
            SortType::ByParentId => {
                self.entries.sort_by(|a,b| {
                    a.parent_id().cmp(&b.parent_id())
                });
                (self.min_mass, self.max_mass) = self.find_min_max_masses();
            },
            SortType::Unsorted => {
                (self.min_mass, self.max_mass) = self.find_min_max_masses();
            }
        }
        self.sort_type = ordering;
    }

    pub fn len(&self) -> usize {
        return self.entries.len()
    }

    pub fn iter(&self) -> std::slice::Iter<T> {
        self.entries.iter()
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<T> {
        self.entries.iter_mut()
    }

    pub fn search_mass(&self, query: f32, error_tolerance: f32, hint_low: Option<usize>, hint_high: Option<usize>) -> Interval {
        let n = self.entries.len();
        let mut lo = hint_low.unwrap_or(0);
        let mut hi = hint_high.unwrap_or(n);

        let mut out = Interval::new(lo, hi);

        while hi != lo {
            let mid = (hi + lo) / 2 as usize;
            let x = self.entries[mid].mass();
            let err = (x - query) / query;
            if lo == hi - 1 || err.abs() <= error_tolerance {
                let mut i = mid;
                while i != usize::MAX {
                    let x = self.entries[i].mass();
                    if (x - query).abs() / query > error_tolerance {
                        i += 1;
                        break;
                    }
                    if i == 0 {
                        break;
                    } else {
                        i -= 1;
                    }
                }
                out.start = i;
                i = mid;
                while i < n {
                    let x = self.entries[i].mass();
                    if (x - query).abs() / query > error_tolerance {
                        break;
                    }
                    i += 1;
                }
                i += 1;
                out.end = (i).min(n);
                return out;
            }
            else if err > 0f32 {
                hi = mid;
            } else if err < 0f32 {
                lo = mid;
            }
        }
        return out
    }

    pub fn search_parent_id(&self, parent_id_range: Interval) -> Interval {
        let mut result = Interval::default();

        let mut index = match self.entries.binary_search_by(|e| e.parent_id().cmp(&parent_id_range.start)) {
            Ok(found) =>  found ,
            Err(location) => location
        };

        while index >= 1 && index < self.len() {
            if parent_id_range.contains(self.entries[index - 1].parent_id()) {
                index -= 1;
            }
            break;
        }
        result.start = index;

        let end = if parent_id_range.end > 0 { parent_id_range.end - 1 } else { 0 };
        index = match self.entries.binary_search_by(|e| e.parent_id().cmp(&end)) {
            Ok(found) =>  found ,
            Err(location) => location
        };

        while index < self.len() {
            if parent_id_range.contains(self.entries[index + 1].parent_id()) {
                index += 1;
            }
            break;
        }
        result.end = index;

        result
    }
}


impl<T: IndexSortable + Default> Index<usize> for IndexBin<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.entries[index]
    }
}



#[cfg(test)]
mod test {
    use super::*;
    use crate::parent::Spectrum;

    #[test]
    fn test_build() {
        let spectra = vec![
            Spectrum::new(2300.0, 2, 0, 0),
            Spectrum::new(2301.0, 2, 0, 1),
            Spectrum::new(2401.0, 2, 1, 0),
            Spectrum::new(4100.0, 4, 0, 2),
        ];
        let mut parent_list = IndexBin::new(spectra, SortType::Unsorted, 0f32, 0f32);
        parent_list.sort(SortType::ByMass);
        assert!(parent_list.len() == 4);

        let search_out = parent_list.search_mass(2300.0, 5e-6, None, None);
        assert!(search_out.start == 0);
        assert!(search_out.end == 1);

        let search_out = parent_list.search_mass(2401.0, 5e-6, None, None);
        assert!(search_out.start == 2);
        assert!(search_out.end == 3);
    }
}