# %%
# 4. Phylogeny, clustering and MAL/MPH copy number
import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import seaborn as sns
from Bio import SeqIO, Phylo
from Bio.Align import PairwiseAligner, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.Phylo.Consensus import get_support
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score

os.makedirs("data/results", exist_ok=True); os.makedirs("figures", exist_ok=True)
GROUPS = ["Beer (Gallone)", "Ale", "Lager", "Wheat/Hefeweizen", "Wine", "Sake", "Shochu/Awamori", "Spirits", "Baking",
          "Bioethanol", "Laboratory", "Wild", "Clinical", "Other/Unknown"]
PALETTE = dict(zip(GROUPS, sns.color_palette("tab10").as_hex() + ["#000000", "#556b2f", "#b8860b", "#00008b"]))
POF_COL = {"POF-": "#000000", "POF+": "#dddddd", "?": "#ff69b4"}
N_BOOT = int(os.environ.get("N_BOOT", 100))
rng = np.random.default_rng(0)

# %%
# helpers (pure functions; tests/test_phylo.py runs this cell on its own)
aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                          open_gap_score=-5, extend_gap_score=-0.5)

def pct_identity(a, b):
    """% identical columns in the best global alignment (gaps count as differences)."""
    aln = aligner.align(a, b)[0]
    return 100 * aln.counts().identities / aln.length

def reference_guided_msa(ref, seqs):
    """Align every sequence to the reference and project it onto reference coordinates.
    Insertions relative to the reference get their own gap-padded columns, so the result is a
    true MSA (all rows same length). shortcut: no realignment between non-reference sequences;
    fine for near-identical alleles, use Clustal Omega/MUSCLE for divergent sets."""
    L = len(ref)
    rows = {}
    for name, s in seqs.items():
        aln = aligner.align(ref, s)[0]
        main, ins, k = ["-"] * L, [""] * (L + 1), 0   # ins[k] = bases inserted before ref position k
        for r, q in zip(aln[0], aln[1]):
            if r == "-": ins[k] += q
            else: main[k] = q; k += 1
        rows[name] = (main, ins)
    width = [max(len(rows[n][1][k]) for n in rows) for k in range(L + 1)]
    out = {}
    for n, (main, ins) in rows.items():
        parts = []
        for k in range(L + 1):
            parts.append(ins[k].ljust(width[k], "-"))
            if k < L: parts.append(main[k])
        out[n] = "".join(parts)
    return MultipleSeqAlignment([SeqRecord(Seq(v), id=n) for n, v in out.items()])

def msa_array(msa):
    return np.array([np.frombuffer(str(r.seq).encode(), dtype=np.uint8) for r in msa])

def p_distance(A):
    """Pairwise p-distance between MSA rows: differing columns / columns not gapped in both
    (= 1 - alignment-induced identity). Same idea as DistanceCalculator('identity'), vectorised with numpy."""
    gap = A == ord("-")
    D = np.zeros((len(A), len(A)))
    for i in range(len(A)):
        both = gap[i] & gap
        D[i] = ((A[i] != A) & ~both).sum(1) / (~both).sum(1)
    return D

def to_dm(names, D):
    return DistanceMatrix(list(names), [list(D[i, :i + 1]) for i in range(len(names))])

def one_hot_variable(msa):
    """One-hot encode only the MSA columns that are not constant."""
    A = np.array([list(str(r.seq)) for r in msa])
    var = [j for j in range(A.shape[1]) if len(set(A[:, j])) > 1]
    df = pd.DataFrame(A[:, var], index=[r.id for r in msa], columns=var)
    return pd.get_dummies(df.astype(str)).astype(int)

def build_trees(msa):
    names = [r.id for r in msa]
    dm = to_dm(names, p_distance(msa_array(msa)))
    c = DistanceTreeConstructor()
    nj, up = c.nj(dm), c.upgma(dm)
    nj.root_at_midpoint()
    for t in (nj, up):
        for cl in t.get_nonterminals(): cl.name = None   # drop 'Inner1' labels so they don't print
    return dm, nj, up

# %%
# 4a. Load alleles, strain names, groups and LOF / POF calls
strains = pd.read_csv("data/metadata/strains.csv").set_index("sacc")
summ = pd.read_csv("data/results/allele_summary.csv")
callable_ = {(g, s): bool(c) for g, s, c in zip(summ.gene, summ.sacc, summ.callable)}
lof_of = {(g, s): bool(l) for g, s, l in zip(summ.gene, summ.sacc, summ.LOF) if str(l) in ("True", "False")}
lofallele = {(g, s): a for g, s, a in zip(summ.gene, summ.sacc, summ.LOF_allele) if isinstance(a, str)}
pof = pd.read_csv("data/results/pof_genotype_by_strain.csv").set_index("strain").predicted_POF.to_dict()

def strain_group(sid):
    if sid == "S288C_REF": return "S288C", "Laboratory"
    key = sid if sid in strains.index else re.sub(r"_\d+$", "", sid)
    if key not in strains.index: return sid, "Other/Unknown"
    st, g = strains.loc[key, ["strain", "group"]]
    return str(st), (g if g in PALETTE else "Other/Unknown")

alleles, meta = {}, []
for gene in ["PAD1", "FDC1"]:
    recs = [r for r in SeqIO.parse(f"data/seqs/{gene}_alleles.fasta", "fasta")
            if not re.search("PARTIAL|CONTAINS_N", r.description) and callable_.get((gene, r.id), True)]
    alleles[gene] = {r.id: str(r.seq).upper() for r in recs}
    for r in recs:
        st, gr = strain_group(r.id)
        meta.append(dict(gene=gene, sacc=r.id, strain=st, group=gr, LOF=lof_of.get((gene, r.id), False),
                         LOF_allele=lofallele.get((gene, r.id), ""), POF=pof.get(st, "POF+" if r.id == "S288C_REF" else "?")))
meta = pd.DataFrame(meta)

# Concatenated PAD1+FDC1 for strains with both (first accession per gene; strain names from strains.csv)
gm = meta.drop_duplicates("sacc").set_index("sacc")
p = pd.read_csv("data/seqs/paired_strains.csv")
pairs = [("S288C", "S288C_REF", "S288C_REF")] + [
    (gm.loc[a, "strain"], a, b) for a, b in zip(p.pad1_sacc.str.split(";").str[0], p.fdc1_sacc.str.split(";").str[0])
    if a in alleles["PAD1"] and b in alleles["FDC1"]]
pairs = list({st: (st, a, b) for st, a, b in pairs}.values())
cid = lambda st: re.sub(r"[^A-Za-z0-9.]+", "_", st)   # safe leaf names for Newick
alleles["PAD1_FDC1"] = {cid(st): alleles["PAD1"][a] + alleles["FDC1"][b] for st, a, b in pairs}
fm = meta[meta.gene == "FDC1"].set_index("sacc")
cmeta = pd.DataFrame([dict(gene="PAD1_FDC1", sacc=cid(st), strain=st, group=gm.loc[a, "group"],
                           LOF=bool(meta[(meta.gene == "PAD1") & (meta.sacc == a)].LOF.any() or fm.loc[b, "LOF"]),
                           LOF_allele=" + ".join(x for x in [meta[(meta.gene == "PAD1") & (meta.sacc == a)].LOF_allele.iloc[0],
                                                              fm.loc[b, "LOF_allele"]] if x),
                           POF=gm.loc[a, "POF"]) for st, a, b in pairs])
meta = pd.concat([meta, cmeta], ignore_index=True)
meta.to_csv("data/results/phylo_allele_groups.csv", index=False)
info = {(r.gene, r.sacc): r for r in meta.itertuples()}
dup = meta.groupby("gene").strain.apply(lambda v: set(v[v.duplicated()]))   # strains with >1 allele (e.g. RIB6002)
label = lambda gene, x: info[gene, x].strain + (f" ({x})" if info[gene, x].strain in dup[gene] else "")
print(meta.groupby(["group", "gene"]).size().unstack().fillna(0).astype(int))
print("POF status (concatenated set):", cmeta.POF.value_counts().to_dict())

# %%
# 4b. Reference-guided MSA and alignment similarity heatmap
msas, ident, outliers = {}, {}, {}
for gene, seqs in alleles.items():
    ref = seqs["S288C_REF"] if "S288C_REF" in seqs else seqs[cid("S288C")]
    msa = reference_guided_msa(ref, seqs)
    names = [r.id for r in msa]
    M = pd.DataFrame(100 * (1 - p_distance(msa_array(msa))), index=names, columns=names)
    # divergent outliers (median identity < 97%, e.g. a non-S. cerevisiae allele) would flatten every tree; report and drop
    out = [x for x in names if M.loc[x].drop(x).median() < 97]
    for x in out:
        print(f"{gene}: excluded divergent allele {label(gene, x)} [{info[gene, x].group}], median identity {M.loc[x].drop(x).median():.2f}%")
        outliers[gene, x] = M.loc[x].drop(x)
    keep = [x for x in names if x not in out]
    M = M.loc[keep, keep]
    msa = MultipleSeqAlignment([r for r in msa if r.id in keep]); names = keep
    alleles[gene] = {x: seqs[x] for x in keep}
    msas[gene] = msa
    with open(f"data/results/{gene}_msa.fasta", "w") as fh: SeqIO.write(msa, fh, "fasta")
    ident[gene] = M
    M.to_csv(f"data/results/{gene}_pct_identity.csv")
    chk = [tuple(rng.choice(names, 2, replace=False)) for _ in range(10)]
    err = max(abs(M.loc[a, b] - pct_identity(seqs[a], seqs[b])) for a, b in chk)
    off = M.values[np.triu_indices(len(M), 1)]
    print(f"{gene}: n={len(M)}, MSA {len(msa)} x {msa.get_alignment_length()}; identity min {off.min():.2f}%, "
          f"median {np.median(off):.2f}%; max |MSA - PairwiseAligner| on 10 pairs = {err:.3f} points")

    stripes = pd.DataFrame({"group": [PALETTE[info[gene, x].group] for x in names],
                            ("POF" if gene == "PAD1_FDC1" else f"{gene} LOF"):
                                [POF_COL[info[gene, x].POF] if gene == "PAD1_FDC1" else ("#000000" if info[gene, x].LOF else "#dddddd")
                                 for x in names]}, index=names)
    Z = linkage(squareform(100 - M.values, checks=False), "average")
    cg = sns.clustermap(M, row_linkage=Z, col_linkage=Z, row_colors=stripes, col_colors=stripes, cmap="viridis",
                        vmin=max(np.percentile(off, 1), 90), vmax=100, figsize=(12, 11), xticklabels=False,
                        yticklabels=len(M) <= 40, cbar_kws={"label": "% identity"})
    cg.fig.suptitle(f"Alignment similarity heatmap: {gene.replace('_', '+')} (n={len(M)})", y=1.02)
    present = [g for g in GROUPS if any(info[gene, x].group == g for x in names)]
    hs = [plt.Rectangle((0, 0), 1, 1, color=PALETTE[g]) for g in present]
    ls = present
    if gene == "PAD1_FDC1": hs += [plt.Rectangle((0, 0), 1, 1, color=c) for c in POF_COL.values()]; ls += list(POF_COL)
    else: hs += [plt.Rectangle((0, 0), 1, 1, color="#000000")]; ls += [f"{gene} LOF"]
    cg.ax_heatmap.legend(hs, ls, loc="upper left", bbox_to_anchor=(1.22, 1.25), fontsize=8)
    cg.savefig(f"figures/heatmap_identity_{gene}.png", dpi=150, bbox_inches="tight")
    plt.show()

# %%
# 4c. Haplotypes
haps = {}
for gene, M in ident.items():
    names = list(M.index)
    hid = pd.Series(M.round(8).apply(tuple, axis=1).factorize()[0], index=names)
    order = hid.value_counts().index           # H0 = most common haplotype
    hid = hid.map({h: i for i, h in enumerate(order)})
    h = pd.DataFrame({"sacc": names, "haplotype": hid.values, "strain": [label(gene, x) for x in names],
                      "group": [info[gene, x].group for x in names], "LOF": [info[gene, x].LOF for x in names],
                      "LOF_allele": [info[gene, x].LOF_allele for x in names], "POF": [info[gene, x].POF for x in names]})
    haps[gene] = h
    h.to_csv(f"data/results/{gene}_haplotypes.csv", index=False)
    tab = pd.crosstab(h.haplotype, h.group)
    tab["n"] = tab.sum(axis=1); tab["LOF"] = h.groupby("haplotype").LOF.sum()
    tab["LOF_allele"] = h.groupby("haplotype").LOF_allele.agg(lambda v: "; ".join(sorted(set(x for x in v if x)))[:60])
    print(f"\n{gene}: {h.haplotype.nunique()} haplotypes among {len(h)} strains/alleles; "
          f"haplotypes per group: {h.groupby('group').haplotype.nunique().to_dict()}")
    print(tab[tab.n >= 3].to_string())

# %%
# 4d. Phylogeny: Neighbor-Joining and UPGMA
def draw(tree, gene, ax, title, fs=5):
    names = {x: label(gene, x) for x in alleles[gene]}
    cols = {names[x]: PALETTE[info[gene, x].group] for x in alleles[gene]}
    Phylo.draw(tree, label_func=lambda c: names.get(c.name, c.name) if c.is_terminal() else "",
               label_colors=lambda s: cols.get(s, "grey"), axes=ax, do_show=False)
    for t in ax.texts: t.set_fontsize(fs)
    tr = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    for y, cl in enumerate(tree.get_terminals(), start=1):   # Phylo.draw puts terminal i at y = i + 1
        r = info[gene, cl.name]
        dark = (r.POF == "POF-") if gene == "PAD1_FDC1" else r.LOF
        ax.scatter([-0.03, -0.015], [y, y], c=[PALETTE[r.group], "#000000" if dark else "#dddddd"],
                   marker="s", s=12, transform=tr, clip_on=False)
    ax.set_title(title); ax.set_ylabel(""); ax.set_yticks([])

trees = {}
for gene, msa in msas.items():
    dm, nj, up = build_trees(msa)
    trees[gene] = dict(dm=dm, nj=nj, up=up)
    Phylo.write(nj, f"data/results/{gene}_nj.nwk", "newick")
    Phylo.write(up, f"data/results/{gene}_upgma.nwk", "newick")
    n = len(msa)
    fig, axs = plt.subplots(1, 2, figsize=(20, max(8, 0.11 * n)))
    draw(nj, gene, axs[0], f"{gene.replace('_', '+')} Neighbor-Joining (midpoint root), n={n}")
    draw(up, gene, axs[1], f"{gene.replace('_', '+')} UPGMA")
    present = [g for g in GROUPS if any(info[gene, x].group == g for x in alleles[gene])]
    fig.legend(handles=[plt.Line2D([], [], color=PALETTE[g], lw=6, label=g) for g in present]
               + [plt.Line2D([], [], color="k", marker="s", lw=0, label="POF- (strip)" if gene == "PAD1_FDC1" else "LOF allele (strip)")],
               loc="upper center", ncol=8, fontsize=9)
    fig.tight_layout(rect=(0.02, 0, 1, 0.97))
    fig.savefig(f"figures/tree_{gene}.png", dpi=150); plt.show()

# %%
# Haplotype NJ tree with bootstrap support (concatenated PAD1+FDC1)
gene = "PAD1_FDC1"
h = haps[gene]
rep = h.drop_duplicates("haplotype").set_index("haplotype").sacc          # one representative per haplotype
A = msa_array(msas[gene])
idx = {r.id: i for i, r in enumerate(msas[gene])}
Ah = A[[idx[s] for s in rep]]
hn = [f"H{k}" for k in rep.index]
nj_h = DistanceTreeConstructor().nj(to_dm(hn, p_distance(Ah)))
reps = [DistanceTreeConstructor().nj(to_dm(hn, p_distance(Ah[:, rng.integers(0, Ah.shape[1], Ah.shape[1])])))
        for _ in range(N_BOOT)]
nj_h = get_support(nj_h, reps)                     # confidence = % of replicates containing the clade
nj_h.root_at_midpoint()
for cl in nj_h.get_nonterminals(): cl.name = None
Phylo.write(nj_h, f"data/results/{gene}_haplotype_nj_bootstrap.nwk", "newick")

def hap_label(name):
    k = int(name[1:]); sub = h[h.haplotype == k]
    comp = ", ".join(f"{g} {c}" for g, c in sub.group.value_counts().items())
    pofs = sub.POF.value_counts().to_dict()
    return f"{name} (n={len(sub)}: {comp}) {'POF-' if pofs.get('POF-', 0) == len(sub) else ''}"
major = lambda name: PALETTE[h[h.haplotype == int(name[1:])].group.mode()[0]]
fig, ax = plt.subplots(figsize=(16, max(8, 0.22 * len(hn))))
Phylo.draw(nj_h, label_func=lambda c: hap_label(c.name) if c.is_terminal() else "",
           label_colors=lambda s: major(s.split()[0]) if s.startswith("H") else "grey", axes=ax, do_show=False)
for t in ax.texts: t.set_fontsize(7)
ax.set_title(f"PAD1+FDC1 haplotype NJ tree, bootstrap % of {N_BOOT} (label colour = majority group)")
ax.set_yticks([]); ax.set_ylabel("")
fig.tight_layout(); fig.savefig("figures/tree_haplotypes_PAD1_FDC1_bootstrap.png", dpi=150); plt.show()
strong = [c for c in nj_h.get_nonterminals() if c.confidence is not None and c.confidence >= 70]
print(f"{len(hn)} haplotypes; {len(strong)} clades with >= 70% bootstrap support")
for c in strong:
    ks = [int(t.name[1:]) for t in c.get_terminals()]
    if len(ks) < len(hn):
        sub = h[h.haplotype.isin(ks)]
        print(f"  {c.confidence:.0f}%: {len(ks)} haplotypes, {len(sub)} strains: {sub.group.value_counts().to_dict()}; POF {sub.POF.value_counts().to_dict()}")

# %%
# 4e. Where do the 'off-group' strains fall?
gene = "PAD1_FDC1"
M = ident[gene]; L = msas[gene].get_alignment_length()
for q in ["laboratory001", "wine012", "wine014", "sake002", "bioethanol002", "Y7", "RIB6002", "N85"]:
    hits = [x for x in M.index if info[gene, x].strain == q] or [x for x in alleles["FDC1"] if info["FDC1", x].strain == q]
    if not hits:
        o = [(g, x) for g, x in outliers if info[g, x].strain == q]
        for g, x in o:
            d = (100 - outliers[g, x]).sort_values()
            print(f"{q} [{info[g, x].group}; {info[g, x].POF}] ({g}) EXCLUDED as divergent: median identity "
                  f"{outliers[g, x].median():.2f}%, closest {info[g, d.index[0]].strain} at {100 - d.iloc[0]:.2f}%")
        if not o: print(f"{q}: not in the data")
        continue
    x = hits[0]
    if x not in M.index:   # not paired -> use FDC1 alone
        gx, Mx = "FDC1", ident["FDC1"]
    else: gx, Mx = gene, M
    r = info[gx, x]
    d = ((100 - Mx.loc[x].drop(x)) / 100 * msas[gx].get_alignment_length()).round().astype(int).sort_values()
    nn = ", ".join(f"{info[gx, y].strain} [{info[gx, y].group}, {info[gx, y].POF}] {v}" for y, v in d.head(5).items())
    print(f"{q} [{r.group}; {r.POF}; {r.LOF_allele or 'intact'}] ({gx}, haplotype H{haps[gx].set_index('sacc').haplotype[x]}, "
          f"{Mx.loc[x].drop(x).median():.2f}% median identity to others)\n    nearest (diffs): {nn}")

# %%
# 4f. Dimensionality reduction and KMeans
gene = "PAD1_FDC1"
X = one_hot_variable(msas[gene])
names = list(X.index)
grp = [info[gene, x].group for x in names]
pofs = [info[gene, x].POF for x in names]
nc = min(10, X.shape[0] - 1, X.shape[1])
pca = PCA(n_components=nc).fit(X.values)
P = pca.transform(X.values)
Mds = MDS(n_components=2, metric="precomputed", init="classical_mds", random_state=0).fit_transform(
    1 - ident[gene].loc[names, names].values / 100)

ks = range(1, min(15, len(set(map(tuple, X.values)))) + 1)
fits = {k: KMeans(k, n_init=10, random_state=0).fit(P) for k in ks}
inertia = [fits[k].inertia_ for k in ks]
sil = {k: silhouette_score(P, fits[k].labels_) for k in ks if 1 < k < len(names)}
best_k = max(sil, key=sil.get) if sil else 1
km = fits[best_k].labels_
known = [i for i, v in enumerate(pofs) if v != "?"]
ari_g = adjusted_rand_score(grp, km)
ari_p = adjusted_rand_score([pofs[i] for i in known], km[known])

fig, axs = plt.subplots(1, 3, figsize=(21, 7))
axs[0].plot(list(ks), inertia, "o-"); axs[0].set(xlabel="k", ylabel="inertia", title=f"KMeans elbow (chosen k={best_k} by silhouette)")
for ax, Y, t in [(axs[1], P[:, :2], f"PCA of one-hot variable sites ({X.shape[1]} features)"), (axs[2], Mds, "MDS on p-distance")]:
    jit = rng.normal(0, 0.01 * (np.ptp(Y, axis=0) + 1e-9), Y.shape)   # identical haplotypes overlap
    for g in [g for g in GROUPS if g in grp]:
        for pv, mk in [("POF+", "o"), ("POF-", "X"), ("?", "s")]:
            i = [n for n in range(len(names)) if grp[n] == g and pofs[n] == pv]
            if i: ax.scatter(Y[i, 0] + jit[i, 0], Y[i, 1] + jit[i, 1], c=PALETTE[g], marker=mk, s=45, edgecolor="k",
                             linewidth=0.4, label=f"{g} {pv}")
    for n, x in enumerate(names):
        if info[gene, x].strain in ["laboratory001", "wine012", "wine014", "sake002", "bioethanol002", "Y7", "RIB6002", "N85", "S288C"]:
            ax.annotate(info[gene, x].strain, Y[n] + jit[n], fontsize=7)
    ax.set_title(t)
ev = pca.explained_variance_ratio_
axs[1].set(xlabel=f"PC1 ({100 * ev[0]:.0f}%)", ylabel=f"PC2 ({100 * ev[1]:.0f}%)")
axs[2].legend(fontsize=7, title="group / POF (o = POF+, X = POF-)", bbox_to_anchor=(1, 1), ncol=1)
fig.suptitle(f"PAD1+FDC1: strains in sequence space (n={len(names)})"); fig.tight_layout()
fig.savefig(f"figures/pca_mds_kmeans_{gene}.png", dpi=150); plt.show()

clusters = pd.DataFrame({"sacc": names, "strain": [info[gene, x].strain for x in names], "group": grp, "POF": pofs,
                         "LOF_allele": [info[gene, x].LOF_allele for x in names], "kmeans": km, "PC1": P[:, 0], "PC2": P[:, 1]})
clusters.to_csv(f"data/results/{gene}_pca_kmeans.csv", index=False)
print(f"PCA explained variance (first 5): {np.round(100 * ev[:5], 1)} %")
print(f"best k = {best_k}, silhouette = {sil[best_k]:.2f}; ARI(kmeans vs group) = {ari_g:.2f}; "
      f"ARI(kmeans vs POF, n={len(known)}) = {ari_p:.2f}")
# silhouette keeps rising when points are discrete haplotypes, so show agreement for every k, not just the chosen one
print(pd.DataFrame([dict(k=k, silhouette=sil[k], ARI_group=adjusted_rand_score(grp, fits[k].labels_),
                         ARI_POF=adjusted_rand_score([pofs[i] for i in known], fits[k].labels_[known])) for k in sil]).round(2).to_string(index=False))
print(pd.crosstab(clusters.group, clusters.kmeans))
print(pd.crosstab(clusters.POF, clusters.kmeans))
print("\nLOF alleles per cluster (beer = 'Beer (Gallone)' + Ale/Lager/Wheat):")
beer = clusters[clusters.group.isin(["Beer (Gallone)", "Ale", "Lager", "Wheat/Hefeweizen"])]
print(pd.crosstab(beer.LOF_allele.replace("", "intact").str.replace(r" -> \d+ aa", "", regex=True), beer.kmeans).to_string())

# %%
# 4g. MAL / MPH maltose-transporter copy number
cnv_path = "data/cnv/mal_copy_number.csv"
if os.path.exists(cnv_path):
    cnv = pd.read_csv(cnv_path)
    genes = [g for g in ["MAL11", "MAL31", "MPH2_3", "total_loci"] if g in cnv.columns]
    cnv["group"] = cnv.group.where(cnv.group.isin(PALETTE), "Other/Unknown")
    cnv["name"] = cnv.strain + np.where(cnv.strain == "Jean-Talon", "*", "") + " [" + cnv.group + "]"
    cnv = cnv.sort_values(["group", "strain"], key=lambda c: c.map(GROUPS.index) if c.name == "group" else c)
    H = cnv.set_index("name")[genes]
    fig, axs = plt.subplots(1, 2, figsize=(18, max(6, 0.3 * len(cnv))), gridspec_kw={"width_ratios": [1, 2]})
    sns.heatmap(H, annot=True, fmt="g", cmap="YlOrRd", ax=axs[0], cbar_kws={"label": "copies"})
    gcol = cnv.set_index("name").group
    for t in axs[0].get_yticklabels(): t.set_color(PALETTE[gcol[t.get_text()]])
    axs[0].set_title("MAL/MPH copy number per genome (* = inferred label)")
    long = cnv.melt(id_vars=["strain", "group"], value_vars=genes, var_name="gene", value_name="copies")
    order = [g for g in GROUPS if g in set(cnv.group)]
    sns.barplot(data=long, x="gene", y="copies", hue="group", hue_order=order, palette=PALETTE, ax=axs[1], errorbar=None, alpha=0.7)
    sns.stripplot(data=long, x="gene", y="copies", hue="group", hue_order=order, palette=PALETTE, dodge=True, ax=axs[1],
                  size=4, edgecolor="k", linewidth=0.5, legend=False)
    axs[1].set_title("Mean copy number by group (bars) and per genome (points)")
    fig.tight_layout(); fig.savefig("figures/mal_copy_number.png", dpi=150); plt.show()

    summary = cnv.groupby("group")[genes].agg(["mean", "min", "max"]).round(2)
    summary.insert(0, "n", cnv.groupby("group").size())
    print(summary.loc[order].to_string())
    dom = cnv.group.isin(["Beer (Gallone)", "Spirits", "Ale", "Lager", "Wheat/Hefeweizen", "Wine", "Sake", "Shochu/Awamori", "Baking", "Bioethanol"])
    wl = cnv.group.isin(["Wild", "Laboratory"])
    rows = []
    for g in genes:
        a, b = cnv.loc[dom, g], cnv.loc[wl, g]
        p = mannwhitneyu(a, b).pvalue if len(a) and len(b) and (a.nunique() + b.nunique() > 2) else np.nan
        rows.append(dict(gene=g, domesticated_mean=a.mean(), domesticated_median=a.median(), domesticated_n=len(a),
                         wild_lab_mean=b.mean(), wild_lab_median=b.median(), wild_lab_n=len(b), mannwhitney_p=p))
    cnv_test = pd.DataFrame(rows).round(3)
    cnv_test.to_csv("data/results/mal_cnv_comparison.csv", index=False)
    print(cnv_test.to_string())
    print("\nBeer genomes:"); print(cnv[cnv.group.isin(["Ale", "Lager"])][["name"] + genes].to_string(index=False))
else:
    print("data/cnv/mal_copy_number.csv not found - CNV section skipped")
