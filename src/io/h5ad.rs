use kira_scio::detect::DetectedFormat;

use crate::input::H5ADInput;
use crate::input::error::InputError;
use crate::io::{RawMatrix, read_via_scio};

pub fn read_h5ad(input: &H5ADInput) -> Result<RawMatrix, InputError> {
    read_via_scio(&input.path, DetectedFormat::H5ad)
}
