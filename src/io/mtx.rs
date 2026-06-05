use kira_scio::detect::DetectedFormat;

use crate::input::TenXInput;
use crate::input::error::InputError;
use crate::io::{RawMatrix, read_via_scio};

pub fn read_tenx(input: &TenXInput) -> Result<RawMatrix, InputError> {
    read_via_scio(&input.root, DetectedFormat::Mtx10x)
}
