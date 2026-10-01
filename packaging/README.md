# Packaging

Distribution channels and the files that drive them. Versions in these
files are set at release time (see `docs/src/release-checklist.md`).

| channel | files | status |
| --- | --- | --- |
| GitHub releases (static binaries) | `.github/workflows/release.yml` | automatic on `v*` tags |
| crates.io | `Cargo.toml`, release workflow | automatic on `v*` tags when `CARGO_REGISTRY_TOKEN` is set |
| Docker / BioContainers | `packaging/docker/Dockerfile` | build locally; BioContainers builds from the bioconda recipe |
| bioconda | `packaging/bioconda/meta.yaml`, `build.sh` | recipe ready; submit to bioconda-recipes after the first tagged release |
| nf-core module | `packaging/nf-core/modules/kira/spliceqc/` | module ready; submit to nf-core/modules after the bioconda package exists |
| MultiQC | `kira_spliceqc_mqc.json` (custom content, pipeline mode) | shipped |
| Python wrapper (PyPI) | `python/` | package ready; publish after the first tagged release |
| Zenodo | `.zenodo.json`, `CITATION.cff` | metadata ready; archive created by the GitHub-Zenodo integration on release |

The binary embeds the default geneset catalog, so the copied `resources/`
directory is informational: `--catalog` is only needed for a custom catalog.
