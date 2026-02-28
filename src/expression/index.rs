use crate::input::error::InputError;

#[derive(Debug, Clone)]
pub struct IndexMapping {
    pub sorted_names: Vec<String>,
    pub old_to_new: Vec<u32>,
}

pub fn build_index(names: Vec<String>, is_gene: bool) -> Result<IndexMapping, InputError> {
    let mut pairs: Vec<(String, usize)> =
        names.into_iter().enumerate().map(|(i, n)| (n, i)).collect();
    pairs.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

    let mut sorted_names = Vec::with_capacity(pairs.len());
    let mut old_to_new = vec![0u32; pairs.len()];

    for (new_idx, (name, old_idx)) in pairs.into_iter().enumerate() {
        if new_idx > u32::MAX as usize {
            return if is_gene {
                Err(InputError::GeneIndexOverflow)
            } else {
                Err(InputError::CellIndexOverflow)
            };
        }
        sorted_names.push(name);
        old_to_new[old_idx] = new_idx as u32;
    }

    Ok(IndexMapping {
        sorted_names,
        old_to_new,
    })
}
