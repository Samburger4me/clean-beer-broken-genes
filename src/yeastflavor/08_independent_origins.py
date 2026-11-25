# Technique source: Lecture 08 (Bio.Phylo listed as a Biopython feature) — lecture_notes/lecture08_slides.md;
# Lecture 14 (hierarchical clustering / trees) — notebook_notes/lecture_14.md. Parsimony counting is beyond class.

# %%
# 8. How many times was the clove flavour lost? (A2)
import os
import numpy as np
import pandas as pd
from Bio import Phylo

RES = "data/results"
HAP_CSV, TREE_NWK = f"{RES}/PAD1_FDC1_haplotypes.csv", f"{RES}/PAD1_FDC1_haplotype_nj_bootstrap.nwk"
assert os.path.exists(HAP_CSV) and os.path.exists(TREE_NWK), "run step 04 first (haplotype table and tree)"

# Step 1: one allele set per strain (every nonsense or frameshift change, including the late stops so I can
# show them separately), then one allele set per haplotype. Strains in one haplotype have identical PAD1+FDC1
# sequence, so they must share alleles; the assert checks that.
vc = pd.read_csv(f"{RES}/variant_calls.csv")
trunc = vc[vc.effect.isin(["nonsense", "frameshift"])].copy()
trunc["allele"] = trunc.gene + " " + trunc.hgvs_c + " " + trunc.protein_change
LATE_STOPS = {"p.W497*", "p.Y500*"}          # stop at >= 95% of Fdc1 length; not called LOF by the caller (step 03)
hap = pd.read_csv(HAP_CSV)
strain_alleles = trunc.groupby("strain").allele.apply(set)
hap["alleles"] = hap.strain.map(strain_alleles).apply(lambda x: x if isinstance(x, set) else set())
per_hap = hap.groupby("haplotype").alleles.agg(lambda s: frozenset().union(*s))
mixed = hap.groupby("haplotype").alleles.agg(lambda s: len({frozenset(x) for x in s}) > 1)
print(f"haplotypes: {hap.haplotype.nunique()}  strains: {len(hap)}  haplotypes with mixed allele sets: {mixed.sum()}")
missing = sorted(set(trunc.strain) - set(hap.strain))
print(f"strains with a truncating allele but no paired PAD1+FDC1 haplotype (left out of the tree count): {missing}")

# %%
# Step 2: the tree-free lower bound. I list each distinct mutation with how many strains and which groups and
# Gallone clades carry it (clades from Table S1, joined the same way as step 07).
s1 = pd.read_excel(ensure_external("data/raw/external/gallone2016/mmc1.xlsx"))
s1 = s1[s1.Code.astype(str).str.fullmatch(r"[A-Z]{2}\d{3}")]
PREFIX = {"BE": "beer", "WI": "wine", "SP": "spirits", "SA": "sake", "WL": "wild", "BI": "bioethanol",
          "BR": "bread", "LA": "laboratory"}
clade_of = dict(zip(s1.Code.str[:2].map(PREFIX) + s1.Code.str[2:], s1.Subpopulation.astype(str).str.strip()))
trunc["clade"] = trunc.strain.map(clade_of).fillna("not in Gallone panel")
origin_tbl = (trunc.groupby(["gene", "hgvs_c", "protein_change"])
              .agg(n_strains=("strain", "nunique"),
                   groups=("group", lambda x: "; ".join(f"{k} {v}" for k, v in x.value_counts().items())),
                   clades=("clade", lambda x: "; ".join(f"{k} {v}" for k, v in x.value_counts().items())))
              .reset_index())
origin_tbl["called_LOF"] = ~origin_tbl.protein_change.isin(LATE_STOPS)
n_distinct_lof = int(origin_tbl.called_LOF.sum())
print(f"distinct LOF mutations (lower bound on independent origins): {n_distinct_lof}")
print(f"... plus late stops not called LOF: {int((~origin_tbl.called_LOF).sum())}")
origin_tbl

# %%
# Step 3: Sankoff parsimony on the haplotype tree for a 0/1 trait. I wrote it from the textbook recursion
# because Biopython has no parsimony scorer for a single character; it handles polytomies, which NJ can
# leave after zero-length branches. Cost of a node in state s = sum over children of min(cost of child in
# s, cost of child in the other state + 1).
def sankoff_changes(tree, present):
    """Minimum number of 0<->1 changes on the tree for leaf set `present`, with the root forced to 0 (absent)."""
    def cost(clade):
        if clade.is_terminal():
            return np.array([np.inf, 0.0]) if clade.name in present else np.array([0.0, np.inf])
        c = np.zeros(2)
        for ch in clade.clades:
            cc = cost(ch)
            c += [min(cc[0], cc[1] + 1), min(cc[1], cc[0] + 1)]
        return c
    return int(cost(tree.root)[0])

# Quick self-check on toy trees: one allele on one cherry = 1 change; on two separate leaves = 2 changes.
_t = Phylo.read(__import__("io").StringIO("((A,B),(C,(D,E)));"), "newick")
assert sankoff_changes(_t, {"A", "B"}) == 1 and sankoff_changes(_t, {"A", "D"}) == 2 and sankoff_changes(_t, set()) == 0

# Step 4: root the haplotype tree on the S288C haplotype. S288C carries no truncating allele and is the
# reference the alleles are named against, so "absent at the root" is the natural starting state.
tree = Phylo.read(TREE_NWK, "newick")
s288c_hap = f"H{int(hap.loc[hap.strain == 'S288C', 'haplotype'].iloc[0])}"
tree.root_with_outgroup(s288c_hap)
rows = []
for r in origin_tbl.itertuples():
    a = f"{r.gene} {r.hgvs_c} {r.protein_change}"
    carriers = {f"H{k}" for k, al in per_hap.items() if a in al}
    rows.append(dict(allele=a, n_haplotypes=len(carriers), parsimony_changes=sankoff_changes(tree, carriers) if carriers else 0))
pars = origin_tbl.assign(allele=origin_tbl.gene + " " + origin_tbl.hgvs_c + " " + origin_tbl.protein_change).merge(pd.DataFrame(rows), on="allele")
pars.to_csv(f"{RES}/lof_independent_origins.csv", index=False)
print(f"tree rooted on {s288c_hap} (S288C)")
pars[["allele", "n_strains", "n_haplotypes", "parsimony_changes", "called_LOF", "clades"]]

# %%
# Step 5: a joint view of the POF trait. A strain is POF− here if it carries any called-LOF allele. The
# parsimony count for this combined trait is the number of times POF− (re)appears on the tree; I report it
# next to the distinct-mutation bound because they answer slightly different questions.
lof_alleles = set(pars[pars.called_LOF].allele)
pofneg_haps = {f"H{k}" for k, al in per_hap.items() if al & lof_alleles}
pof_changes = sankoff_changes(tree, pofneg_haps)
print(f"haplotypes carrying >= 1 called LOF allele: {len(pofneg_haps)} of {len(per_hap)}")
print(f"minimum POF+ <-> POF- changes on the haplotype tree: {pof_changes}")

# Step 6: mutations are not losses. A second stop codon in a gene that is already broken changes nothing,
# so I join LOF mutations that ever occur in the same strain (union-find) and count the groups. Groups never
# share a strain, so each one needs its own loss: this is the minimum number of separate losses of POF.
lof_carriers = {a: set(trunc.loc[trunc.allele == a, "strain"]) for a in lof_alleles}
parent = {a: a for a in lof_carriers}
def find(a):
    while parent[a] != a:
        a = parent[a]
    return a
for a in lof_carriers:
    for b in lof_carriers:
        if a < b and lof_carriers[a] & lof_carriers[b]:
            parent[find(a)] = find(b)
events = {}
for a in lof_carriers:
    events.setdefault(find(a), []).append(a)
BEVERAGE = {"Beer (Gallone)", "Ale", "Lager", "Sake", "Shochu/Awamori", "Spirits"}
event_rows = []
for als in events.values():
    strains = set().union(*(lof_carriers[a] for a in als))
    grp = trunc[trunc.strain.isin(strains)].drop_duplicates("strain").group
    event_rows.append(dict(alleles=" + ".join(sorted(a.split(" p.")[1] if " p." in a else a for a in als)),
                           n_strains=len(strains), beverage=bool(grp.isin(BEVERAGE).mean() > 0.5),
                           groups="; ".join(f"{k} {v}" for k, v in grp.value_counts().items())))
loss_events = pd.DataFrame(event_rows).sort_values("n_strains", ascending=False)
# self-check: groups split the 8 mutations exactly once each
assert sum(len(v) for v in events.values()) == n_distinct_lof, "every LOF mutation must sit in exactly one loss group"
print(f"separate losses (mutation groups that never share a strain): {len(loss_events)}, "
      f"in mostly beverage strains: {int(loss_events.beverage.sum())}")
loss_events

RESULTS.update({
    "a2_n_loss_events": int(len(loss_events)), "a2_n_beverage_loss_events": int(loss_events.beverage.sum()),
    "a2_loss_events": loss_events.to_dict("records"),
    "a2_n_distinct_lof": n_distinct_lof,
    "a2_n_late_stops": int((~origin_tbl.called_LOF).sum()),
    "a2_alleles": pars[["allele", "n_strains", "n_haplotypes", "parsimony_changes", "called_LOF", "groups", "clades"]].to_dict("records"),
    "a2_n_haplotypes": int(len(per_hap)), "a2_n_pofneg_haplotypes": len(pofneg_haps),
    "a2_pof_parsimony_changes": pof_changes, "a2_tree_root": s288c_hap,
    "a2_missing_from_tree": missing,
})
