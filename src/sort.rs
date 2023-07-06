
use std::{ops::{Index, Mul}, str::FromStr, error::Error, fmt::Display};

use crate::interval::Interval;


pub type ParentID = u32;
pub type MassType = f32;


pub fn _isclose(x: MassType, y: MassType, rtol: MassType, atol: MassType) -> bool
{
    (x - y).abs() <= (atol + rtol * y.abs())
}

pub fn isclose(x: MassType, y: MassType) -> bool
{
    _isclose(x, y, 1e-5, 1e-8)
}

pub fn aboutzero(x: MassType) -> bool
{
    isclose(x, 0.0)
}


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


#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum Tolerance {
    PPM(MassType),
    Da(MassType)
}

#[derive(Debug, PartialEq, Eq)]
pub enum ToleranceParsingError {
    UnknownUnit,
    InvalidMagnitude
}

impl Display for ToleranceParsingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{:?}", self))
    }
}

impl Error for ToleranceParsingError {}

impl FromStr for Tolerance {
    type Err = ToleranceParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let n = s.len();
        if n <= 2 {
            return Err(ToleranceParsingError::InvalidMagnitude)
        }
        let s = s.to_lowercase();
        if s.ends_with("da") {
            if let Ok(magnitude) = s[0..n-2].parse::<MassType>() {
                Ok(Self::Da(magnitude))
            } else {
                Err(ToleranceParsingError::InvalidMagnitude)
            }
        } else if s.ends_with("ppm") {
            if let Ok(magnitude) = s[0..n-3].parse::<MassType>() {
                Ok(Self::PPM(magnitude))
            } else {
                Err(ToleranceParsingError::InvalidMagnitude)
            }
        } else {
            Err(ToleranceParsingError::UnknownUnit)
        }
    }
}

impl Tolerance {
    pub fn bounds(&self, query: MassType) -> (MassType, MassType) {
        match self {
            Tolerance::PPM(tol) => {
                let width = query * *tol / 1e6;
                (query - width, query + width)
            }
            Tolerance::Da(tol) => {
                (query - *tol, query + *tol)
            }
        }
    }

    pub fn test(&self, query: MassType, reference: MassType) -> bool {
        let (lower_bound, upper_bound) = self.bounds(reference);
        query >= lower_bound && query <= upper_bound
    }

    pub fn format_error(&self, query: MassType, reference: MassType) -> String {
        match self {
            Self::PPM(_tol) => {
                let magnitude = (query - reference) / reference * 1e6;
                format!("{}PPM", magnitude).to_string()
            },
            Self::Da(_tol) => {
                let magnitude = query - reference;
                format!("{}Da", magnitude).to_string()
            }
        }
    }
}

impl Mul<MassType> for Tolerance {
    type Output = Tolerance;

    fn mul(self, rhs: MassType) -> Self::Output {
        match self {
            Self::Da(val) => Self::Da(rhs * val),
            Self::PPM(val) => Self::PPM(rhs * val)
        }
    }
}

pub trait IndexSortable {
    fn mass(&self) -> MassType;
    fn parent_id(&self) -> ParentID;
}



#[derive(Debug, Clone, Default)]
pub struct IndexBin<T: IndexSortable> {
    pub entries: Vec<T>,
    pub sort_type: SortType,
    pub min_mass: MassType,
    pub max_mass: MassType,
}

impl<T: IndexSortable> IndexBin<T> {
    pub fn new(
        entries: Vec<T>,
        sort_type: SortType,
        min_mass: MassType,
        max_mass: MassType,
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

    fn find_min_max_masses(&self) -> (MassType, MassType) {
        let mut min_mass = MassType::INFINITY;
        let mut max_mass = 0.0 as MassType;

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

    pub fn search_mass(&self, query: MassType, error_tolerance: Tolerance) -> Interval {
        let (lower_bound, upper_bound) = error_tolerance.bounds(query);

        let mut lower_i = self.entries.partition_point(|entry| entry.mass() <= lower_bound);
        let mut upper_i = self.entries[lower_i..self.len()].partition_point(|entry| entry.mass() <= upper_bound) + lower_i;

        while lower_i > 0 {
            if error_tolerance.test(query, self.entries[lower_i - 1].mass()) {
                lower_i -= 1;
            } else {
                break;
            }
        }

        while upper_i + 1 < self.len() {
            if error_tolerance.test(query, self.entries[upper_i + 1].mass()) {
                upper_i += 1;
            } else {
                break;
            }
        }

        return Interval::new(lower_i, upper_i);
    }

    pub fn search_parent_id(&self, parent_id_range: Interval) -> Interval {
        let mut result = Interval::default();

        let start = parent_id_range.start as ParentID;
        let mut index = match self.entries.binary_search_by(|e| e.parent_id().cmp(&start)) {
            Ok(found) =>  found ,
            Err(location) => location
        };

        while index >= 1 && index < self.len() {
            if parent_id_range.contains(self.entries[index - 1].parent_id() as usize) {
                index -= 1;
            } else {
                break;
            }
        }
        result.start = index;

        let end = if parent_id_range.end > 0 { parent_id_range.end - 1 } else { 0 } as ParentID;
        index = match self.entries.binary_search_by(|e| e.parent_id().cmp(&end)) {
            Ok(found) =>  found ,
            Err(location) => location
        };

        while index + 1 < self.len() {
            if parent_id_range.contains(self.entries[index + 1].parent_id() as usize) {
                index += 1;
            } else {
                break;
            }
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
        let mut parent_list = IndexBin::new(spectra, SortType::Unsorted, 0.0, 0.0);
        parent_list.sort(SortType::ByMass);
        assert!(parent_list.len() == 4);

        let search_out = parent_list.search_mass(2300.01, "5ppm".parse().unwrap());
        assert!(search_out.start == 0);
        assert!(search_out.end == 1);

        let search_out = parent_list.search_mass(2401.0, "5ppm".parse().unwrap());
        assert!(search_out.start == 2);
        assert!(search_out.end == 3);
    }
}