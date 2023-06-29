use crate::sort::{IndexSortable, ParentID};


#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct DeconvolutedPeak {
    pub mass: f32,
    pub charge: i16,
    pub intensity: f32,
    pub scan_id: ParentID,
}

impl IndexSortable for DeconvolutedPeak {
    fn mass(&self) -> f32 {
        self.mass
    }

    fn parent_id(&self) -> ParentID {
        self.scan_id
    }
}

impl PartialOrd for DeconvolutedPeak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.mass.partial_cmp(&other.mass)
    }
}

impl DeconvolutedPeak {
    pub fn new(mass: f32, charge: i16, intensity: f32, scan_id: ParentID) -> Self {
        Self {
            mass,
            charge,
            intensity,
            scan_id: scan_id,
        }
    }
}


#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MZPeak {
    pub mz: f32,
    pub intensity: f32,
    pub scan_id: ParentID,
}

impl IndexSortable for MZPeak {
    fn mass(&self) -> f32 {
        self.mz
    }

    fn parent_id(&self) -> ParentID {
        self.scan_id
    }
}

impl PartialOrd for MZPeak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.mz.partial_cmp(&other.mz)
    }
}

impl MZPeak {
    pub fn new(mz: f32 , intensity: f32, scan_id: ParentID) -> Self {
        Self {
            mz,
            intensity,
            scan_id: scan_id,
        }
    }
}





#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_creation() {
        let peak = DeconvolutedPeak::new(256.03, 1, 0.0, 300);
        assert!(peak.mass == 256.03);
        assert!(peak.mass() == 256.03);
        assert!(peak.scan_id == 300);
        assert!(peak.parent_id() == 300);
    }
}