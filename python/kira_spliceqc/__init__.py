"""Python wrapper for the kira-spliceqc binary.

The functions here only build command lines, run the binary and read its
outputs back; every number comes from the Rust implementation.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable

import pandas as pd

__all__ = ["run", "simulate", "validate", "reference_build", "binary", "read_results", "__version__"]
__version__ = "0.5.0"

PathLike = str | os.PathLike[str]


def binary() -> str:
    """Path of the kira-spliceqc binary: `KIRA_SPLICEQC_BIN`, else `PATH`."""
    env = os.environ.get("KIRA_SPLICEQC_BIN")
    if env:
        if not Path(env).is_file():
            raise FileNotFoundError(f"KIRA_SPLICEQC_BIN={env} is not a file")
        return env
    found = shutil.which("kira-spliceqc")
    if found is None:
        raise FileNotFoundError(
            "kira-spliceqc binary not found: install it (conda install -c bioconda kira-spliceqc, "
            "cargo install kira-spliceqc, or a GitHub release) or set KIRA_SPLICEQC_BIN"
        )
    return found


def _call(args: Iterable[str]) -> None:
    cmd = [binary(), *args]
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed ({result.returncode}):\n{result.stderr}")


def _opt(flag: str, value) -> list[str]:
    return [] if value is None else [flag, str(value)]


def read_results(out: PathLike, *, pipeline: bool = True) -> tuple[pd.DataFrame, dict]:
    """`cells.tsv` as a DataFrame indexed by cell name and `summary.json` as a dict."""
    base = Path(out)
    if pipeline:
        base = base / "kira-spliceqc"
    cells = pd.read_csv(base / "cells.tsv", sep="\t", index_col="cell_name", low_memory=False)
    summary_path = base / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.is_file() else {}
    return cells, summary


def run(
    input: PathLike | "anndata.AnnData",
    out: PathLike | None = None,
    *,
    layers: PathLike | None = None,
    junctions: PathLike | None = None,
    metadata: PathLike | None = None,
    stratify_by: str | None = None,
    reference: PathLike | None = None,
    catalog: PathLike | None = None,
    min_counts: int | None = None,
    min_genes: int | None = None,
    extended: bool = False,
    experimental_signatures: bool = False,
    threads: int | None = None,
    extra_args: Iterable[str] = (),
):
    """Run `kira-spliceqc run --run-mode pipeline` and return an AnnData.

    `input` is a 10x directory, a STARsolo `Gene/<subset>` directory, an
    `.h5ad` file or an AnnData. `out` defaults to a temporary directory
    (the files are read and discarded). The per-cell columns of `cells.tsv`
    are added to `obs`, `summary.json` to `uns["kira_spliceqc"]`.
    """
    import anndata as ad

    tmp_in: tempfile.TemporaryDirectory | None = None
    tmp_out: tempfile.TemporaryDirectory | None = None
    try:
        if isinstance(input, ad.AnnData):
            tmp_in = tempfile.TemporaryDirectory(prefix="kira-spliceqc-in-")
            input_path = Path(tmp_in.name) / "input.h5ad"
            input.write_h5ad(input_path)
            adata = input.copy()
        else:
            input_path = Path(input)
            adata = None
        if out is None:
            tmp_out = tempfile.TemporaryDirectory(prefix="kira-spliceqc-out-")
            out = tmp_out.name
        args = [
            "run",
            "--input", str(input_path),
            "--out", str(out),
            "--run-mode", "pipeline",
            "--tsv",
            *_opt("--layers", layers),
            *_opt("--junctions", junctions),
            *_opt("--metadata", metadata),
            *_opt("--stratify-by", stratify_by),
            *_opt("--reference", reference),
            *_opt("--catalog", catalog),
            *_opt("--min-counts", min_counts),
            *_opt("--min-genes", min_genes),
            *_opt("--threads", threads),
        ]
        if extended:
            args.append("--extended")
        if experimental_signatures:
            args.append("--experimental-signatures")
        args.extend(extra_args)
        _call(args)
        cells, summary = read_results(out)
        if adata is None:
            adata = _read_input(input_path)
        return _attach(adata, cells, summary)
    finally:
        for tmp in (tmp_in, tmp_out):
            if tmp is not None:
                tmp.cleanup()


def _read_input(path: Path):
    import anndata as ad

    if path.suffix == ".h5ad":
        return ad.read_h5ad(path)
    try:
        import scanpy as sc
    except ImportError as e:  # pragma: no cover
        raise ImportError("reading a 10x directory needs scanpy (pip install 'kira-spliceqc[tenx]')") from e
    return sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=True)


def _attach(adata, cells: pd.DataFrame, summary: dict):
    cells = cells.reindex(adata.obs_names)
    for column in cells.columns:
        adata.obs[column] = cells[column].values
    adata.uns["kira_spliceqc"] = summary
    return adata


def simulate(out: PathLike, *, n_cells: int | None = None, seed: int | None = None, extra_args: Iterable[str] = ()) -> Path:
    """`kira-spliceqc simulate`; returns the dataset directory (truth in `truth.tsv`)."""
    _call(["simulate", "--out", str(out), *_opt("--n-cells", n_cells), *_opt("--seed", seed), *extra_args])
    return Path(out)


def validate(run_dir: PathLike, truth: PathLike, out: PathLike, *, pairs: Iterable[str] = ()) -> dict:
    """`kira-spliceqc validate`; returns the parsed validation report."""
    args = ["validate", "--run", str(run_dir), "--truth", str(truth), "--out", str(out)]
    for pair in pairs:
        args.extend(["--pair", pair])
    _call(args)
    return json.loads(Path(out).read_text())


def reference_build(input: PathLike, out: PathLike, *, stratify_by: str | None = None, extra_args: Iterable[str] = ()) -> Path:
    """`kira-spliceqc reference build`; returns the reference file path."""
    _call(["reference", "build", "--input", str(input), "--out", str(out), *_opt("--stratify-by", stratify_by), *extra_args])
    return Path(out)
