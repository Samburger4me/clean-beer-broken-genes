# Technique source: matplotlib / seaborn plotting from class notebooks — notebook_notes/lecture_05.md,
# notebook_notes/lecture_14.md (heatmaps); Bio.PDB coordinates as in step 05. Figure layout follows Nature specs.

# %%
# 12. Main figures (Nature format)
import io, json, os, re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from Bio import Phylo
from Bio.PDB import PDBParser
from statsmodels.stats.proportion import proportion_confint

RES, OUT = "data/results", "figures/main"
os.makedirs(OUT, exist_ok=True)
MM = 1 / 25.4                       # matplotlib sizes are in inches; Nature gives widths in mm
W2, W1 = 183 * MM, 89 * MM

# Okabe-Ito colours. I checked this exact set with the dataviz palette validator (all checks pass; the
# light ones need direct labels, which every panel has).
C = {"neg": "#D55E00", "pos": "#0072B2", "green": "#009E73", "orange": "#E69F00", "pink": "#CC79A7",
     "sky": "#56B4E9", "grey": "#9a9a9a", "light": "#e6e6e6", "ink": "#222222"}
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5, "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": 300, "mathtext.default": "regular", "figure.dpi": 150, "lines.linewidth": 0.8,
})

def panel(ax_or_fig, letter, x=-0.12, y=1.04):
    """Put a bold lowercase panel letter at the top-left of an axes (axes coordinates)."""
    ax_or_fig.text(x, y, letter, transform=ax_or_fig.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")

def save(fig, name):
    """Save one figure as vector PDF and 300-dpi PNG in figures/main/."""
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/{name}.{ext}", bbox_inches="tight", pad_inches=0.02)
    print(f"saved {OUT}/{name}.pdf/.png  ({fig.get_size_inches()[0] / MM:.0f} mm wide)")

R = RESULTS  # every number on the figures comes from here
geno = pd.read_csv(f"{RES}/pof_genotype_by_strain.csv")
vc = pd.read_csv(f"{RES}/variant_calls.csv")
uni = json.load(open("data/structure/uniprot_features.json"))

# %%
# Figure 1 — the pathway and the study design (schematic)
def box(ax, x, y, w, h, text, fc, ec="none", fs=5.5, color="#222222", weight="normal"):
    """Rounded text box in data coordinates (0-1 axes)."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012", fc=fc, ec=ec, lw=0.6))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=color, weight=weight, linespacing=1.25)

def arrow(ax, x0, y0, x1, y1, color="#444444", lw=0.8):
    """Simple straight arrow between two points."""
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=6))

fig = plt.figure(figsize=(W2, 72 * MM))
axA = fig.add_axes([0.0, 0.0, 0.40, 0.92]); axB = fig.add_axes([0.43, 0.0, 0.57, 0.92])
for ax in (axA, axB):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

# a: ferulic acid -> 4-vinylguaiacol needs Fdc1 AND the prFMN cofactor made by Pad1
box(axA, 0.02, 0.70, 0.30, 0.14, "Ferulic acid\n(from malt / wheat)", C["light"])
box(axA, 0.62, 0.70, 0.36, 0.14, "4-Vinylguaiacol\nclove, 'phenolic' (POF+)", "#fbe3d3", ec=C["neg"])
arrow(axA, 0.33, 0.77, 0.61, 0.77, lw=1.0)
axA.text(0.47, 0.80, "Fdc1\n(decarboxylase)", ha="center", va="bottom", fontsize=5.5, weight="bold", color=C["pos"])
box(axA, 0.02, 0.40, 0.30, 0.13, "FMN + DMAP", C["light"])
box(axA, 0.40, 0.40, 0.22, 0.13, "prFMN\ncofactor", "#dcecf7", ec=C["pos"])
arrow(axA, 0.33, 0.465, 0.39, 0.465)
axA.text(0.355, 0.54, "Pad1", ha="center", fontsize=5.5, weight="bold", color=C["pos"])
arrow(axA, 0.51, 0.54, 0.47, 0.74, color=C["pos"])
axA.text(0.52, 0.30, "Loss of function in PAD1 or FDC1\n→ no 4-VG → 'clean' (POF−) beer",
         ha="center", va="center", fontsize=5.5, color=C["neg"], weight="bold")
axA.text(0.5, 0.13, "Selected in lager, ale and saké yeasts;\nretained in wheat-beer and wild strains",
         ha="center", va="center", fontsize=5.2, color="#555555")
panel(axA, "a", x=0.0, y=0.98)

# b: pipeline, with sample sizes pulled from RESULTS
src = [("NCBI core_nt BLAST\n+ Gallone 2016 genomes\n" f"{R['core_n_strains']} strains", 0.82),
       ("1,011 Yeast Genomes\nread-depth copy number\n" f"{R['a4_n_strains_joined']} strains", 0.52),
       ("PDB 4ZAC (Fdc1–prFMN)\nAlphaFold Pad1\nlager assemblies", 0.22)]
steps = [("Extract CDS, global alignment\nto S288C, call SNV/indel", 0.82),
         ("LOF rule: stop/frameshift\nbefore 95% of protein", 0.62),
         ("Haplotype NJ tree, bootstrap;\nparsimony origin count", 0.42),
         ("Clade-stratified CMH test;\nvs measured POF phenotype", 0.22),
         ("Map truncations onto\nFdc1 structure; MAL/MPH CN", 0.02)]
for t, y in src:
    box(axB, 0.0, y, 0.30, 0.15, t, C["light"], fs=5.2)
for t, y in steps:
    box(axB, 0.40, y, 0.33, 0.15, t, "#dcecf7", fs=5.2)
outs = [("Which strains\nlost POF?", 0.72), ("How many\nindependent losses?", 0.42), ("Why are the\nalleles null?", 0.12)]
for t, y in outs:
    box(axB, 0.82, y, 0.18, 0.15, t, "#fbe3d3", fs=5.2, weight="bold")
for _, y in src:
    arrow(axB, 0.305, y + 0.075, 0.395, 0.495)
for (_, y0), (_, y1) in zip(steps[:-1], steps[1:]):
    arrow(axB, 0.565, y0, 0.565, y1 + 0.15)
for _, y in outs:
    arrow(axB, 0.735, 0.495, 0.815, y + 0.075)
panel(axB, "b", x=0.0, y=0.98)
save(fig, "Fig1_pathway_design")
plt.show()

# %%
# Figure 2 — where the loss-of-function alleles sit in PAD1 and FDC1 (lollipop map)
# Step 1: one lollipop per distinct truncating allele. The protein position of a nonsense/frameshift is the
# codon that carries it; height = number of strains (log scale, because 1 and 71 must both be readable).
trunc = vc[vc.effect.isin(["nonsense", "frameshift"])].copy()
trunc["aa"] = (trunc.cds_pos - 1) // 3 + 1
lol = trunc.groupby(["gene", "hgvs_c", "protein_change", "aa", "effect"]).strain.nunique().reset_index(name="n")
LATE = {"p.W497*", "p.Y500*"}
DOMAINS = {"FDC1": [(14, 102, "UbiD N-term"), (119, 320, "UbiD prFMN-binding"), (328, 464, "UbiD C-term")],
           "PAD1": [(1, 58, "mito. transit"), (50, 201, "Flavoprotein (UbiX-like)")]}

def sites(acc):
    """Binding-site and active-site residue numbers from the UniProt feature list."""
    f = uni[acc].get("features", [])
    b = sorted({p for x in f if x["type"] == "Binding site" for p in range(x["start"], x["end"] + 1)})
    a = sorted({p for x in f if x["type"] == "Active site" for p in range(x["start"], x["end"] + 1)})
    return b, a

fig, axes = plt.subplots(2, 1, figsize=(W2, 70 * MM), gridspec_kw={"height_ratios": [1, 1.25], "hspace": 0.55})
for ax, gene, acc, letter in [(axes[0], "PAD1", "P33751", "a"), (axes[1], "FDC1", "Q03034", "b")]:
    L = uni[acc]["length"]
    # protein backbone + Pfam/transit boxes
    ax.add_patch(Rectangle((1, -0.18), L, 0.36, fc=C["light"], ec="none"))
    for i, (s, e, name) in enumerate(DOMAINS[gene]):
        fc = "#f3f3f3" if "transit" in name else "#cfe3f2"
        ax.add_patch(Rectangle((s, -0.3), e - s, 0.6, fc=fc, ec=C["pos"] if "transit" not in name else C["grey"], lw=0.5))
        ax.text((s + e) / 2, -0.62, name, ha="center", va="top", fontsize=5, color="#444444")
    b, a = sites(acc)
    ax.vlines(b, -0.3, 0.3, color=C["ink"], lw=0.5)
    ax.vlines(a, -0.42, 0.42, color=C["neg"], lw=1.0)
    for p in a:
        ax.text(p - 3, -0.45, f"E{p}" if gene == "FDC1" else str(p), ha="right", va="top", fontsize=5, color=C["neg"])
    # lollipops
    d = lol[lol.gene == gene].sort_values("aa").reset_index(drop=True)
    CLOSE = 25                       # residues; closer lollipops get stacked heights and left/right labels
    prev_h = -9
    for k, r in d.iterrows():
        h = 0.8 + np.log10(r.n) * 1.3
        if k > 0 and r.aa - d.aa[k - 1] < CLOSE:
            h = max(h, prev_h + 0.8)          # raise the stick above its neighbour
        prev_h = h
        left = k + 1 < len(d) and d.aa[k + 1] - r.aa < CLOSE      # label the earlier one on its left side
        col = C["grey"] if r.protein_change in LATE else C["neg"]
        mk = "D" if r.effect == "frameshift" else "o"
        ax.plot([r.aa, r.aa], [0.3, h], color=col, lw=0.6, zorder=1)
        ax.scatter([r.aa], [h], s=10 + 6 * np.log2(r.n + 1), color=col, marker=mk, edgecolor="white", lw=0.4, zorder=3)
        ax.text(r.aa + (-1 if left else 1) * L * 0.008, h + 0.12, f"{r.protein_change[2:]} ({r.n})", fontsize=5,
                ha="right" if left else "left", va="bottom", color=C["ink"])
    ax.set_xlim(0, L + 30); ax.set_ylim(-1.3, 5.0)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlabel(f"{gene} residue (of {L})", labelpad=1)
    panel(ax, letter, x=-0.04, y=0.98)
# legend drawn once, by hand, so marker shapes (nonsense vs frameshift) are explained without colour alone
lx = axes[1].get_xlim()[1]
axes[0].scatter([], [], color=C["neg"], marker="o", s=12, label="nonsense, called LOF")
axes[0].scatter([], [], color=C["neg"], marker="D", s=12, label="frameshift, called LOF")
axes[0].scatter([], [], color=C["grey"], marker="o", s=12, label="late stop (≥95% length)")
axes[0].plot([], [], color=C["ink"], lw=0.8, label="prFMN/FMN/Mn$^{2+}$ binding residue")
axes[0].plot([], [], color=C["neg"], lw=1.2, label="active-site residue")
axes[0].legend(loc="upper right", frameon=False, ncol=2, handlelength=1.2, columnspacing=1.0)
save(fig, "Fig2_lof_lollipop")
plt.show()

# %%
# Figure 3 — POF− genotypes are enriched in domesticated lineages, inside clades too, and match the measured phenotype
fig = plt.figure(figsize=(W2, 105 * MM))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1], wspace=0.55, hspace=0.55)

# a: POF- fraction by ecological group, Wilson 95% CI, n on the right
g = pd.DataFrame(R["core_pofneg_by_group"]).T
g = g[g.n >= 3].copy()          # groups with < 3 strains (ale, lager) are listed in the legend text instead
g["frac"] = g.pof_neg / g.n
lo_hi = [proportion_confint(r.pof_neg, r.n, method="wilson") for r in g.itertuples()]
g["lo"], g["hi"] = [x[0] for x in lo_hi], [x[1] for x in lo_hi]
DOM = {"Beer (Gallone)", "Sake", "Shochu/Awamori", "Spirits"}   # same groups as the Fisher test (baking is not one)
g["dom"] = g.index.isin(DOM)
g = g.sort_values(["dom", "frac"], ascending=[True, True])
ax = fig.add_subplot(gs[0, 0])
y = np.arange(len(g))
ax.barh(y, g.frac, color=[C["neg"] if d else C["pos"] for d in g.dom], height=0.62)
ax.errorbar(g.frac, y, xerr=np.clip([g.frac - g.lo, g.hi - g.frac], 0, None), fmt="none", ecolor=C["ink"], elinewidth=0.5, capsize=1.2)
ax.set_yticks(y, [f"{i}  (n={int(n)})" for i, n in zip(g.index, g.n)])
ax.set_xlim(0, 1.0); ax.set_xlabel("fraction predicted POF−")
fd = R["core_fisher_domestic"]
ax.set_title(f"Strain level: OR {fd['odds_ratio']:.0f}, Fisher P = {fd['fisher_p']:.1e}", loc="left")
panel(ax, "a", x=-0.55)

# b: by Gallone clade
cv = pd.read_csv(f"{RES}/pof_by_gallone_clade.csv").sort_values("frac_POFneg")
ax = fig.add_subplot(gs[0, 1])
y = np.arange(len(cv))
lo_hi = [proportion_confint(r.n_POFneg, r.n, method="wilson") for r in cv.itertuples()]
ax.barh(y, cv.frac_POFneg, color=C["grey"], height=0.62)   # one colour: clades mix niches, so orange/blue would mean nothing here
ax.errorbar(cv.frac_POFneg, y, xerr=np.clip([cv.frac_POFneg - [a for a, _ in lo_hi], [b for _, b in lo_hi] - cv.frac_POFneg], 0, None),
            fmt="none", ecolor=C["ink"], elinewidth=0.5, capsize=1.2)
ax.set_yticks(y, [f"{c}  (n={n})" for c, n in zip(cv.clade, cv.n)])
ax.set_xlim(0, 1.0); ax.set_xlabel("fraction predicted POF−")
ax.set_title(f"Gallone 2016 clades: {R['a1_n_pofneg_clades']} of {R['a1_n_clades']} mostly POF−", loc="left")
panel(ax, "b", x=-0.45)

# c: forest plot, per-clade OR (Haldane +0.5 for display) and the pooled Mantel-Haenszel OR
st = pd.read_csv(f"{RES}/cmh_strata_by_clade.csv")
st = st[st.informative]
ax = fig.add_subplot(gs[1, 0])
rows = []
for r in st.itertuples():
    a_, b_, c_, d_ = [v + 0.5 for v in (r.dom_POFneg, r.dom_POFpos, r.ref_POFneg, r.ref_POFpos)]
    lor, se = np.log(a_ * d_ / (b_ * c_)), np.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    rows.append((f"{r.clade} ({r.dom_POFneg + r.dom_POFpos} dom. / {r.ref_POFneg + r.ref_POFpos} ref.)", lor, se))
y = np.arange(len(rows))[::-1] + 2
for (lab, lor, se), yy in zip(rows, y):
    ax.plot([np.exp(lor - 1.96 * se), np.exp(lor + 1.96 * se)], [yy, yy], color=C["ink"], lw=0.6)
    ax.scatter([np.exp(lor)], [yy], s=12, color=C["grey"], zorder=3)
mh, (cl, ch) = R["a1_cmh_or"], R["a1_cmh_or_ci"]
ax.plot([cl, ch], [0.5, 0.5], color=C["neg"], lw=1.0)
ax.scatter([mh], [0.5], marker="D", s=22, color=C["neg"], zorder=3)
ax.axvline(1, color=C["grey"], lw=0.5, ls="--")
ax.set_xscale("log"); ax.set_xlim(0.01, 2000)
ax.set_yticks(list(y) + [0.5], [r[0] for r in rows] + ["Mantel–Haenszel pooled"])
ax.set_xlabel("odds ratio, POF− (domesticated vs wine/wild/bioethanol)")
ax.set_title(f"Within clades: CMH P = {R['a1_cmh_p']:.1e}", loc="left")
panel(ax, "c", x=-0.75)

# d: measured phenotype (Gallone Table S5) for each allele signature
ap = pd.read_csv(f"{RES}/pof_phenotype_by_allele_signature.csv").set_index("alleles")
ap = ap[ap.n >= 3].sort_values("n")
ax = fig.add_subplot(gs[1, 1])
y = np.arange(len(ap))
ax.barh(y, ap["POF-"], color=C["neg"], height=0.62, label="measured POF−")
ax.barh(y, ap["POF+"], left=ap["POF-"], color=C["pos"], height=0.62, label="measured POF+")
for yy, (a_, b_) in zip(y, zip(ap["POF-"], ap["POF+"])):
    ax.text(a_ + b_ + 1, yy, f"{a_}/{a_ + b_}", va="center", fontsize=5)
ax.set_yticks(y, [s.replace("FDC1 ", "").replace("PAD1 ", "PAD1 ").replace(" + ", "+") for s in ap.index])
ax.set_xlabel("strains (Gallone 2016 phenotype)")
ax.legend(frameon=False, loc="lower right")
sens = R["a3_sensitivity"][0]; spec = R["a3_specificity"][0]
ax.set_title(f"Genotype vs phenotype: sens. {sens:.0%}, spec. {spec:.0%}", loc="left")
panel(ax, "d", x=-0.62)
save(fig, "Fig3_enrichment_clades_phenotype")
plt.show()

# %%
# Figure 4 — the haplotype tree with the independent loss-of-function origins marked (headline figure)
hap = pd.read_csv(f"{RES}/PAD1_FDC1_haplotypes.csv")
tree = Phylo.read(f"{RES}/PAD1_FDC1_haplotype_nj_bootstrap.nwk", "newick")
tree.root_with_outgroup(R["a2_tree_root"])
tree.ladderize()
trunc["allele"] = trunc.gene + " " + trunc.protein_change.str.replace("p.", "", regex=False)
trunc.loc[trunc.hgvs_c == "c.306G>A", "allele"] = "PAD1 W102* (c.306)"
s_all = trunc.groupby("strain").allele.apply(set)
hap["alleles"] = hap.strain.map(s_all).apply(lambda x: x if isinstance(x, set) else set())
per_hap = hap.groupby("haplotype").alleles.agg(lambda s: set().union(*s))

def sankoff_gains(tree, present):
    """Branches (child clades) where the 0/1 trait switches 0->1 in one most-parsimonious reconstruction."""
    cost = {}
    for cl in tree.find_clades(order="postorder"):
        if cl.is_terminal():
            cost[cl] = np.array([np.inf, 0.0]) if cl.name in present else np.array([0.0, np.inf])
        else:
            cost[cl] = sum(np.array([min(cost[ch][0], cost[ch][1] + 1), min(cost[ch][1], cost[ch][0] + 1)]) for ch in cl.clades)
    gains, stack = [], [(tree.root, 0)]
    while stack:
        cl, s = stack.pop()
        for ch in cl.clades:
            t = int(np.argmin([cost[ch][0] + (s != 0), cost[ch][1] + (s != 1)]))
            if s == 0 and t == 1:
                gains.append(ch)
            stack.append((ch, t))
    return gains

# layout: x = cumulative branch length, y = tip order
tips = tree.get_terminals()
ypos = {t: i for i, t in enumerate(tips)}
for cl in tree.find_clades(order="postorder"):
    if not cl.is_terminal():
        ypos[cl] = np.mean([ypos[c] for c in cl.clades])
xpos = {tree.root: 0.0}
for cl in tree.find_clades(order="preorder"):
    for ch in cl.clades:
        xpos[ch] = xpos[cl] + max(ch.branch_length or 0, 0)
xmax = max(xpos.values())

fig, ax = plt.subplots(figsize=(W2, 150 * MM))
for cl in tree.find_clades():
    for ch in cl.clades:
        ax.plot([xpos[cl], xpos[cl]], [ypos[cl], ypos[ch]], color="#555555", lw=0.5)
        ax.plot([xpos[cl], xpos[ch]], [ypos[ch], ypos[ch]], color="#555555", lw=0.5)
    if not cl.is_terminal() and cl.confidence is not None and cl.confidence >= 70 and cl is not tree.root:
        ax.scatter([xpos[cl]], [ypos[cl]], s=6, color=C["ink"], zorder=4)
        ax.text(xpos[cl] - xmax * 0.004, ypos[cl] - 0.35, f"{cl.confidence:.0f}", fontsize=4.5, ha="right", va="bottom", color="#333333")

# tip labels: haplotype, strain count, groups; coloured square = POF status
for t in tips:
    sub = hap[hap.haplotype == int(t.name[1:])]
    # lower-case niche names so "beer 1" (one beer strain) is not misread as the Beer 1 clade
    grp = ", ".join(f"{k.replace(' (Gallone)', '').lower()} ×{v}" for k, v in sub.group.value_counts().head(3).items())
    neg = (sub.POF == "POF-").mean() > 0.5
    ax.plot([xpos[t], xmax * 1.01], [ypos[t]] * 2, color=C["light"], lw=0.4, ls=":", zorder=0)   # leader line tip -> label
    ax.scatter([xmax * 1.02], [ypos[t]], marker="s", s=9, color=C["neg"] if neg else C["pos"])
    ax.text(xmax * 1.045, ypos[t], f"{t.name} (n={len(sub)}): {grp}", fontsize=4.6, va="center")

# origins: one marker per gain, colour per allele (fixed order), label beside the marker
# allele colours avoid the POF-/POF+ square colours (neg/pos) so a star is never read as a POF status
ALLELE_COL = {"FDC1 Q154*": C["ink"], "FDC1 K167fs": C["orange"], "PAD1 W102*": C["pink"], "FDC1 K54*": C["green"],
              "FDC1 R309*": C["sky"], "FDC1 G288fs": "#8C6D1F", "PAD1 Y98*": "#8C6D1F", "PAD1 W102* (c.306)": C["pink"],
              "FDC1 W497*": C["grey"], "FDC1 Y500*": C["grey"]}
# collect every gain first, so alleles gained on the same branch share one marker and one label
gains = {}
for al, col in ALLELE_COL.items():
    carriers = {f"H{k}" for k, s_ in per_hap.items() if al in s_}
    for ch in (sankoff_gains(tree, carriers) if carriers else []):
        gains.setdefault(id(ch), [ch, []])[1].append(al)
n_gain_marks = 0
for ch, als in gains.values():
    par = next(p_ for p_ in tree.find_clades() if ch in p_.clades)
    xm = (xpos[par] + xpos[ch]) / 2 if xpos[ch] > xpos[par] else xpos[ch]
    col = ALLELE_COL[als[0]]
    late = all("W497" in x or "Y500" in x for x in als)
    ax.scatter([xm], [ypos[ch]], marker="|" if late else "*", s=40, color=col, zorder=5, linewidths=0.8)
    lab = " + ".join(x.replace("FDC1 ", "") for x in als)
    ax.text(xm, ypos[ch] - 0.55, lab, fontsize=4.8, color=col, ha="center", va="bottom", weight="bold")
    n_gain_marks += len(als)
ax.set_ylim(len(tips) - 0.3, -1.2); ax.set_xlim(-xmax * 0.02, xmax * 1.55)
ax.axis("off")
ax.plot([0, 0.001], [-0.6] * 2, color=C["ink"], lw=0.8)
ax.text(0.0005, -0.8, "0.001 subst./site", ha="center", va="bottom", fontsize=5)
ax.scatter([], [], marker="*", s=30, color=C["ink"], label="LOF origin (parsimony gain)")
ax.scatter([], [], marker="s", s=9, color=C["neg"], label="haplotype mostly POF−")
ax.scatter([], [], marker="s", s=9, color=C["pos"], label="haplotype mostly POF+")
ax.scatter([], [], s=6, color=C["ink"], label="bootstrap ≥ 70%")
ax.legend(loc="lower left", frameon=False, fontsize=5)
ax.set_title(f"{R['a2_n_haplotypes']} PAD1+FDC1 haplotypes, NJ, rooted on S288C: {R['a2_n_distinct_lof']} LOF mutations, ≥{R['a2_n_loss_events']} separate losses",
             loc="left", fontsize=6.5)
print(f"gain markers drawn: {n_gain_marks}")
save(fig, "Fig4_haplotype_tree_origins")
plt.show()

# %%
# Figure 5 — every common truncation removes the whole cofactor/metal/active site of Fdc1
s = PDBParser(QUIET=True).get_structure("fdc1", "data/structure/4ZAC_A_CA_het.pdb")
ca = [(r.id[1], r["CA"].coord) for r in s[0]["A"] if "CA" in r and r.id[0] == " "]
resn = np.array([n for n, _ in ca]); X = np.array([c for _, c in ca])
het = [a for a in s[0].get_atoms() if a.get_parent().id[0].startswith("H_")]
fmn = np.array([a.coord for a in het if a.get_parent().get_resname() == "4LU"])
mn = np.array([a.coord for a in het if a.get_parent().get_resname() == "MN"])
mu = X.mean(axis=0); _, _, Vt = np.linalg.svd(X - mu, full_matrices=False)
P, Pf, Pm = (X - mu) @ Vt[:2].T, (fmn - mu) @ Vt[:2].T, (mn - mu) @ Vt[:2].T
e285 = P[resn == 285][0]
stc = pd.DataFrame(R["core_structure"])
show = [("FDC1 K54*", 53, "Saké/shochu"), ("FDC1 Q154*", 153, "Ale (Beer 1)"), ("FDC1 K167fs", 166, "Lager / beer"),
        ("FDC1 R309*", 308, "wine014 (private)")]
fig, axes = plt.subplots(1, 5, figsize=(W2, 50 * MM), gridspec_kw={"width_ratios": [1, 1, 1, 1, 1.1], "wspace": 0.15})
for ax, (lab, keep, who), letter in zip(axes[:4], show, "abcd"):
    kept = resn <= keep
    ax.plot(P[:, 0], P[:, 1], color=C["light"], lw=0.6, zorder=1)
    ax.plot(P[kept, 0], P[kept, 1], color=C["pos"], lw=0.7, zorder=2)
    ax.scatter(Pf[:, 0], Pf[:, 1], s=2, color=C["orange"], zorder=3)
    ax.scatter(Pm[:, 0], Pm[:, 1], s=14, color=C["pink"], edgecolor="k", lw=0.3, zorder=4)
    ax.scatter([e285[0]], [e285[1]], marker="*", s=30, color=C["neg"], zorder=5)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{lab.split()[1]}: keeps 1–{keep}\n{who}", fontsize=5.8)
    panel(ax, letter, x=0.0, y=1.0)
# e: binding residues lost per allele (both proteins), straight from the structure table
stc = stc[~stc.allele.str.contains(";")]       # the Q154*+K167fs haplotype is the same truncation as Q154*
stc["lab"] = stc.gene + " " + stc.allele.str.extract(r"(p\.\S+?)(?:\s|$)")[0].str.replace("p.", "", regex=False)
stc = stc.drop_duplicates("lab").sort_values(["gene", "native_residues"])
ax = axes[4]
y = np.arange(len(stc))
ax.barh(y, stc.n_binding_residues_lost / stc.n_binding_residues_total, color=[C["neg"] if c else C["grey"] for c in stc.LOF_called], height=0.6)
ax.set_yticks(y, stc.lab); ax.set_xlim(0, 1.05)
ax.set_xlabel("fraction of binding /\nactive-site residues lost")
panel(ax, "e", x=-0.75, y=1.0)
ax.set_position([ax.get_position().x0 + 0.08, ax.get_position().y0, ax.get_position().width * 0.85, ax.get_position().height])
axes[0].scatter([], [], s=6, color=C["orange"], label="prFMN"); axes[0].scatter([], [], s=12, color=C["pink"], label="Mn$^{2+}$")
axes[0].scatter([], [], marker="*", s=20, color=C["neg"], label="E285 (proton donor)")
axes[0].plot([], [], color=C["pos"], label="retained"); axes[0].plot([], [], color=C["light"], label="lost")
axes[0].legend(loc="upper left", bbox_to_anchor=(0, -0.02), ncol=5, frameon=False, fontsize=5)
save(fig, "Fig5_fdc1_structure_truncations")
plt.show()

# %%
# Figure 6 — maltose-transporter copy number trades off between genes; lager assemblies lack a S. eubayanus pathway copy
# Step 1: per-strain read-depth copy number from step 09 (1,011 genomes). I show beer vs wild strains per gene,
# with a light jitter so stacked integer values stay visible, and the median as a black bar.
cn = pd.read_csv(f"{RES}/mal_cnv_1011_per_strain.csv")
GENES6 = [("MAL31", "MAL31\n(+paralogue)"), ("MPH2/3", "MPH2/\nMPH3"), ("MAL11", "MAL11"), ("Total", "total")]
rng = np.random.default_rng(0)       # fixed seed so the jitter is identical on every run
fig = plt.figure(figsize=(W2, 62 * MM))
gs = fig.add_gridspec(1, 3, width_ratios=[2.3, 1.25, 1.0], wspace=0.55)
ax = fig.add_subplot(gs[0, 0])
for i, (col, lab) in enumerate(GENES6):
    for j, (grp, colr) in enumerate([("Beer", C["neg"]), ("Wild", C["pos"])]):
        v = cn.loc[cn.group == grp, col].dropna().values
        x = i * 3 + j
        ax.scatter(x + rng.uniform(-0.28, 0.28, len(v)), v, s=2.5, color=colr, alpha=0.55, lw=0)
        ax.plot([x - 0.32, x + 0.32], [np.median(v)] * 2, color=C["ink"], lw=1.0)
    key = col.lower().replace("/", "_")
    pbh = R.get(f"a4_beer_vs_wild_{key}_p_bh", R.get(f"a4_beer_vs_wild_{key}_p"))
    ax.text(i * 3 + 0.5, 8.9, f"P$_{{adj}}$ = {pbh:.1g}", ha="center", fontsize=5)   # one row, above the data
ax.set_xticks([i * 3 + 0.5 for i in range(len(GENES6))], [g[1] for g in GENES6])
ax.set_ylabel("copy number (read depth)")
ax.scatter([], [], color=C["neg"], s=8, label=f"beer (n={R['a4_n_beer']})")
ax.scatter([], [], color=C["pos"], s=8, label=f"wild (n={R['a4_n_wild']})")
ax.set_ylim(-0.3, 10.9)
ax.legend(frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02))
ax.set_title("1,011 Yeast Genomes: beer vs wild isolates", loc="left")
panel(ax, "a", x=-0.12)

# Step 2: total copy number for the three beer clades vs the wild clades (clade definition from step 09).
ax = fig.add_subplot(gs[0, 1])
clades = [("11. Ale beer", "Ale beer"), ("7. Mosaic beer", "Mosaic beer"), ("6. African beer", "African beer")]
sets = [(cn.loc[cn.clade == c, "Total"].dropna().values, lab) for c, lab in clades]
sets.append((cn.loc[cn.clade.isin(R["a4_wild_clades"]), "Total"].dropna().values, "wild clades"))
for j, (v, lab) in enumerate(sets):
    colr = C["pos"] if lab == "wild clades" else C["neg"]
    ax.scatter(j + rng.uniform(-0.28, 0.28, len(v)), v, s=2.5, color=colr, alpha=0.6, lw=0)
    ax.plot([j - 0.32, j + 0.32], [np.median(v)] * 2, color=C["ink"], lw=1.0)
ax.set_xticks(range(len(sets)), [lab.replace(" ", "\n") + f"\n(n={len(v)})" for v, lab in sets], fontsize=5)
ax.set_ylabel("total MAL/MPH copies")
ax.set_title("Beer vs wild clades", loc="left")
panel(ax, "b", x=-0.3)

# Step 3: lager subgenome status grid from step 10 (one cell per strain x gene x parental subgenome).
lg = pd.read_csv(f"{RES}/lager_subgenome_pad1_fdc1.csv")
status_col = lg.columns[-1]
ax = fig.add_subplot(gs[0, 2])
STAT = {"intact": C["pos"], "LOF": C["neg"], "absent": C["light"]}
rows6 = [(s_, g_) for s_ in lg.strain.unique() for g_ in ["PAD1", "FDC1"]]
for yi, (s_, g_) in enumerate(rows6):
    for xi, sub in enumerate(["Sc", "Se"]):
        st_ = lg[(lg.strain == s_) & (lg.gene == g_) & (lg.subgenome == sub)][status_col].iloc[0]
        ax.add_patch(Rectangle((xi, yi), 0.92, 0.88, fc=STAT[st_], ec="white"))
        ax.text(xi + 0.46, yi + 0.44, st_ if st_ != "LOF" else "K167fs", ha="center", va="center", fontsize=5,
                color="white" if st_ != "absent" else "#555555")
ax.set_xlim(0, 1.95); ax.set_ylim(len(rows6), 0)
ax.set_xticks([0.46, 1.46], ["S. cerevisiae\nsubgenome", "S. eubayanus\nsubgenome"], fontsize=5)
ax.set_yticks([i + 0.44 for i in range(len(rows6))], [f"{s_} {g_}" for s_, g_ in rows6], fontsize=5)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.tick_params(length=0)
ax.set_title("Lager (S. pastorianus) assemblies", loc="left")
panel(ax, "c", x=-0.75)
save(fig, "Fig6_mal_cnv_lager")
plt.show()
