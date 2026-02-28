#[derive(Debug, Clone)]
pub struct RawMatrix {
    pub genes: Vec<String>,
    pub cells: Vec<String>,
    pub triplets: Vec<(u32, u32, u32)>,
}

pub mod h5ad;
pub mod mtx;
