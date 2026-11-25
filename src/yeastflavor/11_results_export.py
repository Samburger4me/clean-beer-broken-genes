# Technique source: pandas I/O conventions from class notebooks — notebook_notes/lecture_05.md

# %%
# 11. Collect the core numbers into results/results.json
import pandas as pd
from Bio import Phylo

RES = "data/results"

# Step 1: dataset size (strains with a genotype row, records per gene, callable records).
geno = pd.read_csv(f"{RES}/pof_genotype_by_strain.csv")
alle = pd.read_csv(f"{RES}/allele_summary.csv")
RESULTS.update({
    "core_n_strains": int(len(geno)),
    "core_n_records": {g: int((alle.gene == g).sum()) for g in ["PAD1", "FDC1"]},
    "core_n_callable": {g: int(((alle.gene == g) & alle.callable).sum()) for g in ["PAD1", "FDC1"]},
    "core_pof_calls": {k: int(v) for k, v in geno.predicted_POF.value_counts().items()},
    "core_genotype_class": {k: int(v) for k, v in geno.genotype_class.value_counts().items()},
})

# Step 2: the strain-level enrichment tests from section 3.6 (kept, but quoted with the relatedness caveat).
enr = pd.read_csv(f"{RES}/lof_enrichment_tests.csv")
for key, row in zip(["core_fisher_domestic", "core_fisher_beer", "core_chi2_all"], enr.itertuples()):
    RESULTS[key] = {k: (None if pd.isna(v) else float(v)) for k, v in row._asdict().items() if k not in ("Index", "test")}
    RESULTS[key]["test"] = row.test

# Step 3: POF- fraction per ecological group (the numbers behind the enrichment bar chart).
called = geno[geno.predicted_POF.isin(["POF-", "POF+"])]
RESULTS["core_pofneg_by_group"] = {g: {"n": int(len(d)), "pof_neg": int((d.predicted_POF == "POF-").sum())}
                                   for g, d in called.groupby("group")}

# Step 4: recurrent alleles and the Mertens comparison (share of POF- Gallone beers carrying FDC1 Q154*).
rec = pd.read_csv(f"{RES}/recurrent_lof_alleles.csv")
RESULTS["core_recurrent_alleles"] = rec[["gene", "hgvs_c", "protein_change", "n_strains", "n_called_LOF", "groups"]].to_dict("records")
vc = pd.read_csv(f"{RES}/variant_calls.csv")
bneg = set(geno[(geno.group == "Beer (Gallone)") & (geno.predicted_POF == "POF-")].strain)
q154 = set(vc[(vc.gene == "FDC1") & (vc.hgvs_c == "c.460C>T")].strain)
RESULTS["core_q154_in_pofneg_beer"] = {"k": len(bneg & q154), "n": len(bneg)}
assert RESULTS["core_q154_in_pofneg_beer"]["n"] > 0, "no POF- Gallone beers found - check step 03 output"

# Step 5: well-supported clades on the haplotype NJ tree (>= 70% bootstrap is the usual cut-off) and what
# strains they hold, so the paper can quote "91% bootstrap" with its contents.
hap = pd.read_csv(f"{RES}/PAD1_FDC1_haplotypes.csv")
tree = Phylo.read(f"{RES}/PAD1_FDC1_haplotype_nj_bootstrap.nwk", "newick")
BOOT_MIN = 70
n_tips = len(tree.get_terminals())
strong = []
for c in tree.get_nonterminals():
    if c.confidence is not None and c.confidence >= BOOT_MIN and len(c.get_terminals()) < n_tips:
        sub = hap[hap.haplotype.isin([int(t.name[1:]) for t in c.get_terminals()])]
        strong.append({"support": float(c.confidence), "n_haplotypes": len(c.get_terminals()), "n_strains": int(len(sub)),
                       "groups": {k: int(v) for k, v in sub.group.value_counts().items()},
                       "POF": {k: int(v) for k, v in sub.POF.value_counts().items()},
                       "LOF_alleles": sorted({x for x in sub.LOF_allele.dropna()})})
RESULTS["core_boot_clades"] = sorted(strong, key=lambda d: -d["support"])
RESULTS["core_n_haplotypes"] = int(hap.haplotype.nunique())

# Step 6: structural context of each truncation (how much protein is left, how many cofactor/metal/active-site
# residues are lost).
st = pd.read_csv(f"{RES}/structure_lof_context.csv")
RESULTS["core_structure"] = st[["gene", "allele", "native_residues", "fraction_retained", "LOF_called",
                                "n_binding_residues_lost", "n_binding_residues_total"]].to_dict("records")

# Step 7: the original assembly-based MAL/MPH copy-number comparison (superseded by step 09, kept for Methods).
mal = pd.read_csv(f"{RES}/mal_cnv_comparison.csv")
RESULTS["core_mal_assembly"] = mal.to_dict("records")
print(f"RESULTS now holds {len(RESULTS)} keys")
sorted(RESULTS)
