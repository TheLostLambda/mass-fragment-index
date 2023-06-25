use crate::sort::IndexSortable;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ParentMolecule {
    pub mass: f32,
    pub id: usize,
    pub source_id: usize,
    pub start_position: u16,
    pub size: u16,
}

impl ParentMolecule {
    pub fn new(mass: f32, id: usize, parent_id: usize, start_position: u16, size: u16) -> Self {
        Self {
            mass,
            id,
            source_id: parent_id,
            start_position,
            size,
        }
    }
}

impl IndexSortable for ParentMolecule {
    fn mass(&self) -> f32 {
        self.mass
    }

    fn parent_id(&self) -> usize {
        self.source_id
    }
}


#[derive(Debug, Clone, Default, PartialEq)]
pub struct Peptide {
    pub mass: f32,
    pub id: usize,
    pub protein_id: usize,
    pub start_position: u16,
    pub sequence: String,
}

impl Peptide {
    pub fn new(mass: f32, id: usize, protein_id: usize, start_position: u16, sequence: String) -> Self { Self { mass, id, protein_id, start_position, sequence } }
}


impl IndexSortable for Peptide {
    fn mass(&self) -> f32 {
        self.mass
    }

    fn parent_id(&self) -> usize {
        self.protein_id
    }
}


#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Spectrum {
    pub precursor_mass: f32,
    pub precursor_charge: i32,
    pub source_file_id: usize,
    pub scan_number: usize,
}

impl Spectrum {
    pub fn new(
        precursor_mass: f32,
        precursor_charge: i32,
        source_file_id: usize,
        scan_number: usize,
    ) -> Self {
        Self {
            precursor_mass,
            precursor_charge,
            source_file_id,
            scan_number,
        }
    }
}

impl IndexSortable for Spectrum {
    fn mass(&self) -> f32 {
        self.precursor_mass
    }

    fn parent_id(&self) -> usize {
        self.source_file_id
    }
}
