#[cfg(all(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature = "avx2"
))]
const BACKEND: &str = "avx2";

#[cfg(all(
    target_arch = "aarch64",
    target_feature = "neon",
    not(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_feature = "avx2"
    ))
))]
const BACKEND: &str = "neon";

#[cfg(not(any(
    all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_feature = "avx2"
    ),
    all(target_arch = "aarch64", target_feature = "neon")
)))]
const BACKEND: &str = "scalar";

pub const fn backend() -> &'static str {
    BACKEND
}
