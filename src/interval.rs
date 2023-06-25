
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    pub start: usize,
    pub end: usize,
}

impl Interval {
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    pub fn contains(&self, point: usize) -> bool {
        self.start <= point && self.end > point
    }

    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_contains() {
        let iv = Interval::new(5, 10);
        assert!(iv.contains(5));
        assert!(iv.contains(9));
        assert!(!iv.contains(10));
    }
}