use std::cmp::Ordering;

pub fn median(values: &[f32]) -> f32 {
    let mut data: Vec<f32> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if data.is_empty() {
        return f32::NAN;
    }
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let mid = data.len() / 2;
    if data.len() % 2 == 1 {
        data[mid]
    } else {
        (data[mid - 1] + data[mid]) * 0.5
    }
}

pub fn mad(values: &[f32], med: f32) -> f32 {
    let mut devs: Vec<f32> = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .map(|v| (v - med).abs())
        .collect();
    if devs.is_empty() {
        return f32::NAN;
    }
    devs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let mid = devs.len() / 2;
    if devs.len() % 2 == 1 {
        devs[mid]
    } else {
        (devs[mid - 1] + devs[mid]) * 0.5
    }
}
