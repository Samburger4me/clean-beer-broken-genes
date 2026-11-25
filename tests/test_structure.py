"""Self-check for src/yeastflavor/05_structural_analysis.py.
Run from the repo root: python -m pytest tests/  (or python tests/test_structure.py)"""
import json, os, subprocess, sys
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
r = subprocess.run([sys.executable, "src/yeastflavor/05_structural_analysis.py"], capture_output=True, text=True)
assert r.returncode == 0, r.stderr[-2000:]

feat = json.load(open("data/structure/uniprot_features.json"))
p = PDBParser(QUIET=True)
for path, acc in [("data/structure/4ZAC_A_CA_het.pdb", "Q03034"), ("data/structure/AF-P33751-F1_CA.pdb", "P33751")]:
    ch = p.get_structure("x", path)[0]["A"]
    res = [x for x in ch if x.id[0] == " "]
    seq = feat[acc]["sequence"]
    assert len(seq) - 5 <= len(res) <= len(seq), (path, len(res), len(seq))  # allow a few unmodelled residues
    assert all(seq1(x.get_resname()) == seq[x.id[1] - 1] for x in res), f"{path}: residue names differ from UniProt"
ch = p.get_structure("x", "data/structure/4ZAC_A_CA_het.pdb")[0]["A"]
assert sum(1 for a in ch.get_atoms() if a.get_parent().get_resname() == "4LU") >= 30, "prFMN missing"
# known geometry: binding residue R175 CA must be close to prFMN, distal C-terminus far
lu = [a for a in ch.get_atoms() if a.get_parent().get_resname() == "4LU"]
assert min(ch[175]["CA"] - a for a in lu) < 8 and min(ch[503]["CA"] - a for a in lu) > 20

m = pd.read_csv("data/results/structure_missense_context.csv")
d = m.filter(like="d_").to_numpy(dtype=float)
assert np.all(d[~np.isnan(d)] >= 0)
l = pd.read_csv("data/results/structure_lof_context.csv")
assert ((l.fraction_retained > 0) & (l.fraction_retained < 1)).all()
assert (l.n_binding_residues_lost <= l.n_binding_residues_total).all()
q = l[l.allele.str.contains(r"Q154\*")]
assert q.empty or (q.native_residues == 153).all() and (q.n_binding_residues_lost == q.n_binding_residues_total).all()
for f in ["05a_feature_tracks", "05b_fdc1_ca_pca_truncations", "05c_fdc1_contact_map", "05d_missense_distances"]:
    assert os.path.getsize(f"figures/{f}.png") > 10_000, f
print("test_structure: all checks passed")


def test_module_checks_passed():
    """The asserts above run at import; reaching here means they all passed (lets pytest report it)."""
    assert True
