"""Assert-based checks for the helper cell of src/yeastflavor/04_phylogeny_clustering_cnv.py.
Run from the repo root: python -m pytest tests/  (or python tests/test_phylo.py)"""
import random, pathlib
src = pathlib.Path(__file__).resolve().parents[1] / "src/yeastflavor/04_phylogeny_clustering_cnv.py"
cells = src.read_text(encoding="utf-8").split("# %%")
ns = {}
exec(next(c for c in cells if "import PairwiseAligner" in c).replace("os.makedirs", "(lambda *a, **k: None)"), ns)  # imports cell, no dirs created
exec(next(c for c in cells if "# helpers" in c), ns)                  # pure helper cell

random.seed(0)
ref = "".join(random.choice("ACGT") for _ in range(600))
snp = ref[:100] + ("A" if ref[100] != "A" else "C") + ref[101:]
ins = ref[:300] + "T" + ref[300:]

assert ns["pct_identity"](ref, ref) == 100
assert abs(ns["pct_identity"](ref, snp) - 100 * 599 / 600) < 1e-9
assert abs(ns["pct_identity"](ref, ins) - 100 * 600 / 601) < 1e-9

seqs = {"ref": ref, "snp": snp, "ins": ins, "dup": ref}
msa = ns["reference_guided_msa"](ref, seqs)
D = ns["p_distance"](ns["msa_array"](msa))
assert (D == D.T).all() and D[0, 3] == 0
assert abs(100 * (1 - D[0, 1]) - ns["pct_identity"](ref, snp)) < 1e-9   # MSA identity = pairwise identity
assert {len(r.seq) for r in msa} == {601}                       # insertion adds one gap column
assert str(msa[0].seq).replace("-", "") == ref and str(msa[2].seq).replace("-", "") == ins

dm, nj, up = ns["build_trees"](msa)
for t in (nj, up):
    assert sorted(c.name for c in t.get_terminals()) == sorted(seqs)
assert dm["ref", "dup"] == 0 and dm["ref", "snp"] > 0
X = ns["one_hot_variable"](msa)
assert X.shape[0] == 4 and X.shape[1] == 4                      # 2 variable columns x 2 states
print("test_phylo: all checks passed")


def test_module_checks_passed():
    """The asserts above run at import; reaching here means they all passed (lets pytest report it)."""
    assert True
