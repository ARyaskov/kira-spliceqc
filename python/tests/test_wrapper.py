"""End-to-end test of the wrapper; skipped when the binary is not installed."""

import shutil

import pytest

import kira_spliceqc as ksq

pytestmark = pytest.mark.skipif(
    shutil.which("kira-spliceqc") is None and "KIRA_SPLICEQC_BIN" not in __import__("os").environ,
    reason="kira-spliceqc binary not available",
)


def test_simulate_run_validate(tmp_path):
    sim = ksq.simulate(tmp_path / "sim", n_cells=600, seed=3)
    adata = ksq.run(sim, tmp_path / "run", junctions=sim / "sj")
    assert "unspliced_fraction" in adata.obs.columns
    assert adata.uns["kira_spliceqc"]["input"]["n_cells"] == adata.n_obs
    report = ksq.validate(tmp_path / "run", sim / "truth.tsv", tmp_path / "validation.json")
    assert any(r["truth"] == "truth_cryptic" for r in report["results"])
