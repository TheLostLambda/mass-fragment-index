use std::fmt::Debug;

use crate::interval::Interval;
use crate::sort::{IndexBin, IndexSortable, SortType};

#[derive(Debug, Default)]
pub struct SearchIndex<T: IndexSortable + Default, P: IndexSortable + Default> {
    pub bins: Vec<IndexBin<T>>,
    pub parents: IndexBin<P>,
    pub bins_per_dalton: u32,
    pub max_item_mass: f32,
    pub sort_type: SortType,
}

impl<T: IndexSortable + Default, P: IndexSortable + Default> SearchIndex<T, P> {
    pub fn empty(bins_per_dalton: u32, max_fragment_size: f32) -> Self {
        let mut inst = Self {
            bins_per_dalton,
            max_item_mass: max_fragment_size,
            ..Self::default()
        };
        inst.initialize_bins();
        inst
    }

    fn initialize_bins(&mut self) {
        let mut mass_step: f32 = 0.0;
        self.bins = Vec::new();
        while mass_step < self.max_item_mass {
            let bin: IndexBin<T> = IndexBin::default();
            self.bins.push(bin);
            mass_step += 1 as f32 / self.bins_per_dalton as f32;
        }
        let bin: IndexBin<T> = IndexBin::default();
        self.bins.push(bin);
    }

    pub fn describe_bins(&self) {
        for (i, bin, bin_mass) in self.bins.iter().enumerate().map(
                        |(i, bin)| (i, bin, (i as f32) * 1.0 / self.bins_per_dalton as f32)) {
            if bin.len() > 0 {
                println!("Bin {} @ {:0.3} has {} elements", i, bin_mass, bin.len());
            }
        }
    }

    pub fn new(
        bins: Vec<IndexBin<T>>,
        parents: IndexBin<P>,
        bins_per_dalton: u32,
        max_item_mass: f32,
        sort_type: SortType,
    ) -> Self {
        Self {
            bins,
            parents: parents,
            bins_per_dalton,
            max_item_mass,
            sort_type,
        }
    }

    pub fn total_bins_for_mass(&self) -> u32 {
        self.bins_per_dalton * (self.max_item_mass.round() as u32)
    }

    pub fn bin_for_mass(&self, mass: f32) -> usize {
        let i = (mass * self.bins_per_dalton as f32).round() as usize;
        let bin_index = if i > self.bins.len() {
            self.bins.len() - 1
        } else {
            i
        };
        return bin_index;
    }

    pub fn sort(&mut self, ordering: SortType) {
        for bin in self.bins.iter_mut() {
            bin.sort(ordering)
        }
        self.sort_type = ordering
    }

    pub fn add_parent(&mut self, parent_molecule: P) {
        self.parents.append(parent_molecule)
    }

    pub fn add(&mut self, entry: T) -> usize {
        let mass = entry.mass();
        let bin_index = self.bin_for_mass(mass);
        self.bins[bin_index].append(entry);
        bin_index
    }

    pub fn parents_for(&self, mass: f32, error_tolerance: f32) -> Interval {
        let iv = self.parents
            .search_mass(mass, error_tolerance);
        iv
    }

    pub fn parents_for_range(&self, low: f32, high: f32, error_tolerance: f32) -> Interval {
        let mut out = Interval::default();
        out.start = self.parents_for(low, error_tolerance).start;
        out.end = self.parents_for(high, error_tolerance).end;
        out
    }

    pub fn search(
        &self,
        query: f32,
        error_tolerance: f32,
        parent_interval: Option<Interval>,
    ) -> SearchIndexSearcher<T, P> {
        let low = query - (query * error_tolerance);
        let high = query + (query * error_tolerance);

        let mut low_bin = self.bin_for_mass(low);
        let mut high_bin = self.bin_for_mass(high);

        if low_bin != 0 {
            if self.bins[low_bin].max_mass < low {
                low_bin -= 1;
            }
        }

        if high_bin < self.bins.len() - 2 {
            high_bin += 1;
            if self.bins[high_bin].min_mass > high {
                high_bin -= 1;
            }
        }

        let parent_interval_used = match parent_interval {
            Some(iv) => iv,
            None => Interval::new(0, self.parents.len()),
        };

        SearchIndexSearcher::new(
            self,
            query,
            error_tolerance,
            low_bin,
            high_bin,
            low_bin,
            0,
            parent_interval_used,
        )
    }
}

#[derive(Debug, Clone)]
pub struct SearchIndexSearcher<'a, T: IndexSortable + Default, P: IndexSortable + Default> {
    index: &'a SearchIndex<T, P>,
    pub query: f32,
    pub error_tolerance: f32,
    pub low_bin: usize,
    pub high_bin: usize,
    pub current_bin: usize,
    pub bin_position_range: Interval,
    pub bin_position: usize,
    pub parent_id_range: Interval,
}

impl<'a, T: IndexSortable + Default, P: IndexSortable + Default> SearchIndexSearcher<'a, T, P> {
    pub fn new(
        index: &'a SearchIndex<T, P>,
        query: f32,
        error_tolerance: f32,
        low_bin: usize,
        high_bin: usize,
        current_bin: usize,
        position: usize,
        parent_id_range: Interval,
    ) -> Self {
        let mut inst = Self {
            index,
            query,
            error_tolerance,
            low_bin,
            high_bin,
            current_bin,
            bin_position_range: Interval { start: 0, end: 0 },
            bin_position: position,
            parent_id_range,
        };

        // println!("Searching for {} in {}-{}", inst.query, inst.low_bin, inst.high_bin);
        inst.initialize_position_range();
        // println!("{:?} is valid? {}", inst.peek(), inst.is_peek_valid());
        inst
    }

    pub fn peek(&self) -> &T {
        let bin = self.get_current_bin();
        &bin[self.bin_position]
    }

    fn advance(&mut self) -> bool {
        let current_bin = &self.index.bins[self.current_bin];
        let mut hit = false;
        match self.index.sort_type {
            SortType::ByMass => {
                self.bin_position += 1;
                while self.current_bin_has_more() {
                    if self.is_peek_valid() {
                        hit = true;
                        break;
                    }
                    self.bin_position += 1;
                }
                hit
            }
            SortType::ByParentId => {
                self.bin_position += 1;
                for i in self.bin_position..self.bin_position_range.end {
                    self.bin_position = i;
                    if self._test_entry_parent_sort(&current_bin[i]) {
                        hit = true;
                        break;
                    }
                }
                hit
            }
            SortType::Unsorted => {true}
        }
    }

    fn _update_position_range_from_mass_sort(&mut self, current_bin: &IndexBin<T>) -> bool {
        self.bin_position_range = current_bin.search_mass(self.query, self.error_tolerance);
        self.bin_position = self.bin_position_range.start;
        let mut hit = false;
        for i in self.bin_position..current_bin.len() {
            if !self.parent_id_range.contains(current_bin[i].parent_id()) {
                continue;
            } else {
                self.bin_position_range.start = i;
                hit = true;
                break;
            }
        }
        if !hit {
            self.bin_position_range.start = current_bin.len();
        }
        self.bin_position = self.bin_position_range.start;
        hit
    }

    fn _update_position_range_from_parent_id_sort(&mut self, current_bin: &IndexBin<T>) -> bool {
        // Starting from the start of the bin, walk along it sequentially until we find a
        // valid entry.
        self.bin_position_range.start = 0;
        self.bin_position_range.end = current_bin.len();
        let mut hit = false;
        let guess_range = current_bin.search_parent_id(self.parent_id_range);
        let starting_guess = if guess_range.start > 0 {guess_range.start - 1} else { 0 };
        self.bin_position_range.end = (guess_range.end + 1).min(current_bin.len());
        self.bin_position_range.start = starting_guess;
        for i in starting_guess..current_bin.len() {
            if self._test_entry_parent_sort(&current_bin[i]) {
                self.bin_position_range.start = i;
                hit = true;
                break;
            }
        }
        if !hit {
            self.bin_position_range.start = current_bin.len();
        }
        hit
    }

    fn is_peek_valid(&self) -> bool {
        if self.current_bin_has_more() {
            match self.index.sort_type {
                SortType::ByParentId => self._test_entry_parent_sort(self.peek()),
                SortType::ByMass => self._test_entry_mass_sort(self.peek()),
                SortType::Unsorted => true
            }
        } else {
            false
        }
    }

    fn _test_entry_parent_sort(&self, entry: &T) -> bool {
        let entry_parent_id = entry.parent_id();
        if !self.parent_id_range.contains(entry_parent_id) {
            return false;
        }
        let entry_mass = entry.mass();
        (entry_mass - self.query).abs() / self.query < self.error_tolerance
    }

    fn _test_entry_mass_sort(&self, entry: &T) -> bool {
        let entry_parent_id = entry.parent_id();
        if !self.parent_id_range.contains(entry_parent_id) {
            return false;
        }
        let entry_mass = entry.mass();
        (entry_mass - self.query).abs() / self.query < self.error_tolerance
    }

    fn update_bin(&mut self) -> bool {
        let mut hit = false;
        while (self.current_bin <= self.high_bin) && ((self.current_bin + 1) < self.index.bins.len()) {
            self.current_bin += 1;
            let current_bin = &self.index.bins[self.current_bin];
            hit = match self.index.sort_type {
                SortType::ByMass => {
                    self._update_position_range_from_mass_sort(current_bin)
                }
                SortType::ByParentId => {
                    self._update_position_range_from_parent_id_sort(current_bin)
                }
                SortType::Unsorted => {
                    true
                }
            };

            self.bin_position = self.bin_position_range.start;
            if self.current_bin_has_more() {
                break;
            }
        }
        hit
    }

    fn next_entry(&mut self) -> Option<&T> {
        let mut entry = None;
        let mut hit = false;
        if self.bin_position < self.bin_position_range.end {
            entry = Some(&(self.index.bins[self.current_bin])[self.bin_position]);
            hit = self.advance();
        }

        if self.bin_position >= self.bin_position_range.end
            || self.bin_position >= self.get_current_bin().len()
            || !hit
        {
            self.update_bin();
        }

        entry
    }

    fn get_current_bin(&self) -> &IndexBin<T> {
        &self.index.bins[self.current_bin]
    }

    fn current_bin_has_more(&self) -> bool {
        let bin = self.get_current_bin();
        self.bin_position < bin.len()
    }

    fn initialize_position_range(&mut self) -> bool {
        let current_bin = &self.index.bins[self.current_bin];
        self.bin_position_range.start = 0;
        self.bin_position_range.end = current_bin.len();
        self.bin_position = 0;
        match self.index.sort_type {
            SortType::ByMass => {
                let result = current_bin.search_mass(self.query, self.error_tolerance);
                self.bin_position_range = result;
                if !self.advance() {
                    self.bin_position_range.start = self.bin_position_range.end;
                } else {
                    self.bin_position_range.start = self.bin_position;
                }
            }
            SortType::ByParentId => {
                let mut hit = false;
                for i in self.bin_position..current_bin.len() {
                    if self._test_entry_parent_sort(&current_bin[i]) {
                        hit = true;
                        self.bin_position_range.start = i;
                        break;
                    }
                }

                if !hit {
                    self.bin_position_range.start = self.bin_position_range.end;
                }
            }
            SortType::Unsorted => {}
        }
        self.bin_position = self.bin_position_range.start;
        if self.bin_position == self.bin_position_range.end {
            if self.current_bin <= self.high_bin {
                self.current_bin += 1;
                return self.initialize_position_range();
            }
            return true;
        }
        // println!("Starting search in bin {} in range {:?}", self.current_bin, self.bin_position_range);
        return false;
    }
}

impl<'a, T: IndexSortable + Default + Clone + Debug, P: IndexSortable + Default + Clone> Iterator for SearchIndexSearcher<'a, T, P> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(entry) = self.next_entry() {
            Some(entry.clone())
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test{
    use super::*;
    use crate::parent::Spectrum;
    use crate::peak::DeconvolutedPeak;

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

        let peaks = vec![
            DeconvolutedPeak::new(251.5, 1, 0.0, 0),
            DeconvolutedPeak::new(251.5, 1, 0.0, 2),
            DeconvolutedPeak::new(251.6, 1, 0.0, 1),
            DeconvolutedPeak::new(303.7, 1, 0.0, 0),
            DeconvolutedPeak::new(501.2, 1, 0.0, 1)
        ];
        let mut index: SearchIndex<DeconvolutedPeak, Spectrum> = SearchIndex::empty(10, 1000.0);
        for peak in peaks {
            index.add(peak);
        }

        let idx = index.bin_for_mass(251.5);
        assert!(idx == 2515);
    }
}