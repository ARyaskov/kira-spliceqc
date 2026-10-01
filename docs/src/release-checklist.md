# Release checklist

1. `CHANGELOG.md`: move "Unreleased" into a `## [x.y.z] - YYYY-MM-DD` section;
   check "Changed (breaking)" against the [compatibility policy](compatibility.md).
2. Bump `version` in `Cargo.toml`, `python/pyproject.toml`,
   `python/kira_spliceqc/__init__.py`, `packaging/bioconda/meta.yaml`,
   `packaging/nf-core/modules/kira/spliceqc/{main.nf,environment.yml}`,
   `CITATION.cff` (with `date-released`); `cargo update -p kira-spliceqc`
   refreshes `Cargo.lock`.
3. `cargo fmt --check`, `cargo clippy --lib --bins -- -D warnings`,
   `cargo test --no-fail-fast`, `./benchmarks/run_tier1.sh /tmp/bench`.
4. Commit `chore(release): vx.y.z`, tag `vx.y.z`, push the tag. The release
   workflow builds the binaries, drafts the GitHub release from the changelog
   section and publishes to crates.io (`CARGO_REGISTRY_TOKEN`).
5. Zenodo: the GitHub integration archives the release; copy the DOI into
   `README.md`, `CITATION.cff` (`identifiers`) and `packaging/bioconda/meta.yaml`.
6. bioconda: update `sha256` in `packaging/bioconda/meta.yaml` from the tag
   tarball and open a pull request to bioconda-recipes (first release: add the
   recipe). BioContainers images follow automatically.
7. PyPI: `cd python && python -m build && twine upload dist/*`.
8. nf-core: update the container tag in the module and open a pull request to
   nf-core/modules (first release: add the module with tests).
9. Documentation: the docs workflow publishes `docs/` to GitHub Pages on push
   to `main`; check the rendered metric cards and tutorials.
10. Announce: GitHub release notes, scverse Discourse, Biostars, Bluesky; link
    the MultiQC and nf-core integrations.
