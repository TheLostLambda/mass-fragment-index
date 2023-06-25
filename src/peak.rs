use crate::sort::IndexSortable;


#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Peak {
    pub mass: f32,
    pub charge: i16,
    pub intensity: f32,
    pub scan_id: usize,
}

impl IndexSortable for Peak {
    fn mass(&self) -> f32 {
        self.mass
    }

    fn parent_id(&self) -> usize {
        self.scan_id
    }
}

impl PartialOrd for Peak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.mass.partial_cmp(&other.mass)
    }
}

impl Peak {
    pub fn new(mass: f32, charge: i16, intensity: f32, scan_id: usize) -> Self {
        Self {
            mass,
            charge,
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
        let peak = Peak::new(256.03, 1, 0.0, 300);
        assert!(peak.mass == 256.03);
        assert!(peak.mass() == 256.03);
        assert!(peak.scan_id == 300);
        assert!(peak.parent_id() == 300);
    }
}