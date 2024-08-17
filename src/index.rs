#[cfg(feature = "serde")]
use serde::{Serialize, Deserialize};

#[cfg(feature = "parallelism")]
use rayon::prelude::*;

use crate::interval::Interval;
use crate::sort::{IndexBin, IndexSortable, MassType, SortType, Tolerance};

#[derive(Debug, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SearchIndex<T: IndexSortable + Default, P: IndexSortable + Default> {
    pub bins: Vec<IndexBin<T>>,
    pub parents: IndexBin<P>,
    pub(crate) bins_per_dalton: u32,
    pub(crate) max_item_mass: MassType,
    pub(crate) sort_type: SortType,
}

impl<T: IndexSortable + Default, P: IndexSortable + Default> SearchIndex<T, P> {
    pub fn empty(bins_per_dalton: u32, max_fragment_size: MassType) -> Self {
        let mut inst = Self {
            bins_per_dalton,
            max_item_mass: max_fragment_size,
            ..Self::default()
        };
        inst.initialize_bins();
        inst
    }

    fn initialize_bins(&mut self) {
        let mut mass_step: MassType = 0.0;
        self.bins = Vec::new();
        while mass_step < self.max_item_mass {
            let bin: IndexBin<T> = IndexBin::default();
            self.bins.push(bin);
            mass_step += 1 as MassType / self.bins_per_dalton as MassType;
        }
        let bin: IndexBin<T> = IndexBin::default();
        self.bins.push(bin);
    }

    pub fn num_bins(&self) -> usize {
        self.bins.len()
    }

    pub fn num_entries(&self) -> usize {
        self.bins.iter().map(|b| b.len()).sum()
    }

    pub fn iter_bins(&self) -> std::slice::Iter<'_, IndexBin<T>> {
        self.bins.iter()
    }

    pub fn describe_bins(&self) {
        for (i, bin, bin_mass) in self.bins.iter().enumerate().map(|(i, bin)| {
            (
                i,
                bin,
                (i as MassType) * 1.0 / self.bins_per_dalton as MassType,
            )
        }) {
            if bin.len() > 0 {
                println!("Bin {} @ {:0.3} has {} elements", i, bin_mass, bin.len());
            }
        }
    }

    pub fn new(
        bins: Vec<IndexBin<T>>,
        parents: IndexBin<P>,
        bins_per_dalton: u32,
        max_item_mass: MassType,
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

    pub fn bin_for_mass(&self, mass: MassType) -> usize {
        let i = (mass * self.bins_per_dalton as MassType).round() as usize;
        let bin_index = if i >= self.bins.len() {
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

    #[cfg(feature = "parallelism")]
    pub fn par_sort(&mut self, ordering: SortType) where T: Send {
        self.bins.par_iter_mut().for_each(|bin| bin.sort(ordering));
        self.sort_type = ordering;
    }

    pub fn add_parent(&mut self, parent_molecule: P) {
        self.parents.push(parent_molecule)
    }

    pub fn add(&mut self, entry: T) -> usize {
        let mass = entry.mass();
        let bin_index = self.bin_for_mass(mass);
        self.bins[bin_index].push(entry);
        bin_index
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

    pub fn search(
        &self,
        query: MassType,
        error_tolerance: Tolerance,
        parent_interval: Option<Interval>,
    ) -> SearchIndexSearcher<T, P> {
        let (low, high) = error_tolerance.bounds(query);
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

    pub fn bins_per_dalton(&self) -> u32 {
        self.bins_per_dalton
    }

    pub fn max_item_mass(&self) -> f32 {
        self.max_item_mass
    }

    pub fn sort_type(&self) -> SortType {
        self.sort_type
    }
}

#[derive(Debug, Clone)]
pub struct SearchIndexSearcher<'a, T: IndexSortable + Default, P: IndexSortable + Default> {
    index: &'a SearchIndex<T, P>,
    pub query: MassType,
    pub error_tolerance: Tolerance,
    pub low_bin: usize,
    pub high_bin: usize,
    pub current_bin: usize,
    pub bin_position_range: Interval,
    pub bin_position: usize,
    pub parent_id_range: Interval,
    pub query_spans_bin: bool,
}

impl<'a, T: IndexSortable + Default, P: IndexSortable + Default> SearchIndexSearcher<'a, T, P> {
    pub fn new(
        index: &'a SearchIndex<T, P>,
        query: MassType,
        error_tolerance: Tolerance,
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
            query_spans_bin: false,
        };

        inst.initialize_position_range();
        inst
    }

    /// Update the `query_spans_bin` flag to indicate whether we can assume that the mass accuracy
    /// requirement is guaranteed and we can skip checking it.
    fn _update_bin_state(&mut self) {
        if self.current_bin == self.index.bins.len() - 1 {
            self.query_spans_bin = false
        } else {
            let current_bin = self.get_current_bin();
            let (lo, hi) = self.min_max_mass_acceptable();
            if lo <= current_bin.min_mass && hi >= current_bin.max_mass {
                self.query_spans_bin = true
            } else {
                self.query_spans_bin = false
            }
        }
    }

    /// Compute the minimum and maximum acceptable masses
    fn min_max_mass_acceptable(&self) -> (MassType, MassType) {
        self.error_tolerance.bounds(self.query)
    }

    /// Peek at the next item this iterator would yield.
    /// This might be unsafe because it doesn't check if the
    /// next position is valid
    pub fn peek(&self) -> &T {
        let bin = self.get_current_bin();
        &bin[self.bin_position]
    }

    /// Advance the iterator to the next valid bin position
    /// and possibly the next bin, or terminate in an unusable
    /// state
    fn advance(&mut self) -> bool {
        let current_bin = &self.index.bins[self.current_bin];
        let mut hit = false;
        hit = match self.index.sort_type {
            SortType::ByMass => {
                self.bin_position += 1;
                while self.current_bin_has_more() && self.bin_position < self.bin_position_range.end
                {
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
                    if self._test_entry(&current_bin[i]) {
                        hit = true;
                        break;
                    }
                }
                hit
            }
            SortType::Unsorted => true,
        };
        if !hit {
            self.bin_position = self.bin_position_range.end;
        }
        hit
    }

    /// Update the `bin_position_range` and `bin_position` for a new bin when
    /// sorted by product mass
    fn _update_position_range_from_mass_sort(&mut self, current_bin: &IndexBin<T>) -> bool {
        self.bin_position_range = current_bin.search_mass(self.query, self.error_tolerance);
        self.bin_position = self.bin_position_range.start;
        let mut hit = false;
        for i in self.bin_position..current_bin.len() {
            if !self.parent_id_range.contains(current_bin[i].parent_id() as usize) {
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

    /// Update the `bin_position_range` and `bin_position` for a new bin when
    /// sorted by parent id
    fn _update_position_range_from_parent_id_sort(&mut self, current_bin: &IndexBin<T>) -> bool {
        // Starting from the start of the bin, walk along it sequentially until we find a
        // valid entry.
        self.bin_position_range.start = 0;
        self.bin_position_range.end = current_bin.len();
        let mut hit = false;
        let guess_range = current_bin.search_parent_id(self.parent_id_range);
        let starting_guess = if guess_range.start > 0 {
            guess_range.start - 1
        } else {
            0
        };
        self.bin_position_range.end = (guess_range.end + 1).min(current_bin.len());
        self.bin_position_range.start = starting_guess;
        for i in starting_guess..current_bin.len() {
            if self._test_entry(&current_bin[i]) {
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

    /// Peek if the next position is valid
    fn is_peek_valid(&self) -> bool {
        if self.current_bin_has_more() {
            match self.index.sort_type {
                SortType::ByParentId => self._test_entry(self.peek()),
                SortType::ByMass => self._test_entry(self.peek()),
                SortType::Unsorted => true,
            }
        } else {
            false
        }
    }

    /// Test if the entry is valid under parent id sorting
    #[inline(always)]
    fn _test_entry(&self, entry: &T) -> bool {
        let entry_parent_id = entry.parent_id();
        if !self.parent_id_range.contains(entry_parent_id as usize) {
            return false;
        }
        let entry_mass = entry.mass();
        self.error_tolerance.test(self.query, entry_mass)
        // if self.query_spans_bin {
        //     true
        // } else {
        //     let entry_mass = entry.mass();
        //     self.error_tolerance.test(self.query, entry_mass)
        // }
    }

    fn update_bin(&mut self) -> bool {
        let mut hit = false;
        self.query_spans_bin = false;
        while (self.current_bin <= self.high_bin)
            && ((self.current_bin + 1) < self.index.bins.len())
        {
            self.current_bin += 1;
            let current_bin = &self.index.bins[self.current_bin];
            hit = match self.index.sort_type {
                SortType::ByMass => self._update_position_range_from_mass_sort(current_bin),
                SortType::ByParentId => {
                    self._update_position_range_from_parent_id_sort(current_bin)
                }
                SortType::Unsorted => true,
            };

            self.bin_position = self.bin_position_range.start;
            if self.current_bin_has_more() {
                self._update_bin_state();
                break;
            }
        }
        hit
    }

    fn next_entry(&mut self) -> Option<&'a T> {
        let mut entry = None;
        let mut hit = false;
        if self.bin_position < self.bin_position_range.end {
            assert!(self._test_entry(&(self.index.bins[self.current_bin])[self.bin_position]));
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
        self.query_spans_bin = false;
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
                let guess_range = current_bin.search_parent_id(self.parent_id_range);
                let starting_guess = if guess_range.start > 0 {
                    guess_range.start - 1
                } else {
                    0
                };
                self.bin_position_range.end = (guess_range.end + 1).min(current_bin.len());
                self.bin_position_range.start = starting_guess;

                for i in starting_guess..current_bin.len() {
                    if self._test_entry(&current_bin[i]) {
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
        self._update_bin_state();
        return false;
    }
}

impl<'a, T: IndexSortable + Default, P: IndexSortable + Default> Iterator
    for SearchIndexSearcher<'a, T, P>
{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(entry) = self.next_entry() {
            Some(entry)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::parent::Spectrum;
    use crate::peak::DeconvolutedPeak;

    #[test]
    fn test_build() {
        let spectra = vec![
            Spectrum::new(2300.0, 2, 0, 0, 0),
            Spectrum::new(2301.0, 2, 0, 1, 1),
            Spectrum::new(2401.0, 2, 1, 0, 2),
            Spectrum::new(4100.0, 4, 0, 2, 3),
        ];
        let mut parent_list = IndexBin::new(spectra, SortType::Unsorted, 0.0, 0.0);
        parent_list.sort(SortType::ByMass);

        let peaks = vec![
            DeconvolutedPeak::new(251.5, 1, 0.0, 0),
            DeconvolutedPeak::new(251.5, 1, 0.0, 2),
            DeconvolutedPeak::new(251.6, 1, 0.0, 1),
            DeconvolutedPeak::new(303.7, 1, 0.0, 0),
            DeconvolutedPeak::new(501.2, 1, 0.0, 1),
        ];
        let mut index: SearchIndex<DeconvolutedPeak, Spectrum> = SearchIndex::empty(10, 1000.0);
        for peak in peaks {
            index.add(peak);
        }

        let idx = index.bin_for_mass(251.5);
        assert!(idx == 2515);
    }
}
