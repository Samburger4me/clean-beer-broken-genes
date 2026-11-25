"""Cross-checks on results/results.json (run after `python run_all.py`): numbers the paper quotes must agree
with each other, so a silent change in one step can't leave the manuscript inconsistent."""
import json, pathlib

R = json.loads((pathlib.Path(__file__).resolve().parents[1] / "results/results.json").read_text())


def test_genotype_phenotype_confusion_adds_up():
    c = R["a3_confusion"]
    assert sum(c.values()) == R["a3_n_compared"]
    assert abs(R["a3_sensitivity"][0] - c["TP"] / (c["TP"] + c["FN"])) < 1e-12


def test_cmh_strata_are_informative_clades():
    assert R["a1_cmh_n_strata"] == len(R["a1_cmh_strata"]) >= 2
    lo, hi = R["a1_cmh_or_ci"]
    assert lo < R["a1_cmh_or"] < hi


def test_origin_count_matches_allele_table():
    called = [a for a in R["a2_alleles"] if a["called_LOF"]]
    assert len(called) == R["a2_n_distinct_lof"]
    assert all(a["parsimony_changes"] >= 1 for a in R["a2_alleles"])


def test_strain_level_fisher_table_matches_group_counts():
    f = R["core_fisher_domestic"]
    g = R["core_pofneg_by_group"]
    ref = sum(g[k]["pof_neg"] for k in ["Wild", "Wine", "Bioethanol"])
    assert f["b_POFneg"] == ref
