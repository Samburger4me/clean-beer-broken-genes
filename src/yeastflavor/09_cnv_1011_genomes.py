# Technique source: Notebook 05 (pandas, plotting) + Notebook 12 (adjusted p-values) — notebook_notes/lecture_05.md, notebook_notes/lecture_12.md; rank tests, BH and Cliff's delta are beyond class

# %%
# 9. MAL/MPH transporter copy number across the 1,011 S. cerevisiae genomes
# Libraries used in this step
import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kruskal       # rank tests (beyond class)
from statsmodels.stats.multitest import multipletests  # Benjamini-Hochberg (beyond class)

# Input files (downloaded earlier, see data/MANIFEST.csv)
cn_path = ensure_external("data/raw/external/peter2018_1002genomes/genesMatrix_CopyNumber.tab.gz")
meta_path = ensure_external("data/raw/external/peter2018_nature_supp/41586_2018_30_MOESM3_ESM.xls")
os.makedirs("data/results", exist_ok=True)
os.makedirs("figures/extended", exist_ok=True)

# The matrix columns for my genes. I checked these names by searching the header for the
# systematic IDs; the '_NumOfGenes_2' suffix marks a column that merges two paralogs.
gene_cols = {
    "MAL11": "X4188.YGR289C",                    # AGT1, the maltotriose transporter
    "MAL31": "X2289.YBR298C_NumOfGenes_2",       # MAL31 + paralog, merged
    "MPH2/3": "X2689.YDL247W_NumOfGenes_2",      # MPH2 + MPH3, merged
}
pof_cols = {
    "FDC1": "X3229.YDR539W",                     # ferulic acid decarboxylase
    "PAD1 (merged)": "X3228.YDR538W_NumOfGenes_2",  # PAD1 + YDR537C, merged
}

# %%
# Step 1: load the copy-number matrix. I read only the columns I need because the full matrix
# has 7,796 ORF columns and I would not use the rest.
cn = pd.read_csv(cn_path, sep="\t",
                 usecols=["STD_name"] + list(gene_cols.values()) + list(pof_cols.values()))
cn = cn.rename(columns={v: k for k, v in {**gene_cols, **pof_cols}.items()})

# The total is the sum of the three transporter columns: all the maltose/maltotriose
# transporter copy number I can see in this matrix.
cn["Total"] = cn[list(gene_cols)].sum(axis=1)

# I expect 1,011 strains with one row each; duplicates would double-count strains in the tests.
assert cn.shape[0] == 1011, f"expected 1011 strains, got {cn.shape[0]}"
assert cn["STD_name"].is_unique, "duplicate strain IDs in the copy-number matrix"
assert cn[list(gene_cols) + list(pof_cols)].notna().all().all(), "missing copy-number values"
print("Copy-number matrix (strains, columns):", cn.shape)

# Step 2: load Table S1. The real header is on spreadsheet row 4 (header=3), and several
# column names and cells carry stray spaces, so I strip them before matching.
meta = pd.read_excel(meta_path, sheet_name="Table S1", header=3, engine="xlrd")
meta.columns = [str(c).strip() for c in meta.columns]
for col in ["Standardized name", "Ecological origins", "Clades"]:
    meta[col] = meta[col].astype("string").str.strip()

# The sheet ends with blank footer rows; I drop rows without a strain name.
meta = meta.dropna(subset=["Standardized name"])
assert meta["Standardized name"].is_unique, "duplicate strain names in Table S1"

# Step 3: join. An inner join keeps strains present in both files; every matrix strain should
# have metadata, otherwise group sizes would silently shrink.
df = cn.merge(meta[["Standardized name", "Ecological origins", "Clades"]],
              left_on="STD_name", right_on="Standardized name", how="inner")
assert len(df) == len(cn), "some copy-number strains have no Table S1 row"
print("Strains after join:", len(df))

# %%
# Raw ecological-origin labels, before I group them. I print these first so the mapping
# below can be checked against what the table really says.
eco_counts = df["Ecological origins"].fillna("(blank)").value_counts()
eco_counts

# %%
# Step 1: map each raw label to my group. Anything not in the map (Unknown, Lab strain, blank)
# becomes NaN and drops out of the tests.
eco_map = {
    "Beer": "Beer", "Wine": "Wine", "Sake": "Sake", "Bakery": "Bakery",
    "Distillery": "Distillery", "Bioethanol": "Bioethanol",
    "Palm wine": "Other fermentation", "Cider": "Other fermentation",
    "Fermentation": "Other fermentation", "Dairy": "Other fermentation",
    "Industrial": "Other fermentation", "Probiotic": "Other fermentation",
    "Tree": "Wild", "Soil": "Wild", "Insect": "Wild", "Fruit": "Wild",
    "Nature": "Wild", "Flower": "Wild", "Water": "Wild",
    "Human, clinical": "Clinical", "Human": "Clinical",
}
group_order = ["Beer", "Wine", "Sake", "Bakery", "Distillery", "Bioethanol",
               "Other fermentation", "Wild", "Clinical"]
domesticated = group_order[:7]
df["group"] = df["Ecological origins"].map(eco_map)

# Step 2: make sure I did not miss a label. Every raw label must be either mapped or on my
# deliberate exclusion list, otherwise strains would vanish without me noticing.
excluded_labels = {"Unknown", "Lab strain"}
unmapped = set(df["Ecological origins"].dropna()) - set(eco_map) - excluded_labels
assert not unmapped, f"labels with no group: {unmapped}"
df["domesticated"] = df["group"].isin(domesticated)
print("Strains in a group:", df["group"].notna().sum(), "of", len(df))
print("Strains excluded (unknown/lab/blank):", df["group"].isna().sum())

group_n = df["group"].value_counts().reindex(group_order)
group_n

# %%
# Clades: I keep the Table S1 clade labels, and merge the four "1. Wine/European (subclade k)"
# labels into one wine clade because they are one lineage in Peter 2018.
df["clade"] = df["Clades"].str.replace(r"\s*\(subclade \d\)", "", regex=True)

# Beer-related clades named in Table S1.
beer_clades = ["11. Ale beer", "7. Mosaic beer", "6. African beer"]

# I call a clade "wild" from the data, not from memory: a clade is wild when at least
# WILD_CLADE_FRAC of its strains have a wild ecological origin (my own cut-off, chosen so that
# a clade must be clearly wild, not merely mixed).
WILD_CLADE_FRAC = 0.75
clade_wild_frac = (df.dropna(subset=["clade"])
                     .groupby("clade")["group"].apply(lambda g: (g == "Wild").mean()))
wild_clades = sorted(clade_wild_frac[clade_wild_frac >= WILD_CLADE_FRAC].index)
assert all(c in set(df["clade"]) for c in beer_clades), "a beer clade name is not in Table S1"
assert len(wild_clades) > 0, "no clade passed the wild cut-off"
df["clade_group"] = np.where(df["clade"].isin(wild_clades), "Wild clades", df["clade"])
print("Wild clades:", ", ".join(wild_clades))
print("Strains in wild clades:", df["clade"].isin(wild_clades).sum())

clade_wild_frac.loc[wild_clades + beer_clades].round(2)

# %%
def cliffs_delta(x, y):
    """Cliff's delta of x vs y (= rank-biserial correlation, 2U/(n1*n2) - 1); positive means x
    tends to have more copies than y."""
    u = mannwhitneyu(x, y, alternative="two-sided").statistic
    return 2 * u / (len(x) * len(y)) - 1

def slug(s):
    """Lower-case key fragment for RESULTS, e.g. '11. Ale beer' -> 'ale_beer'."""
    return re.sub(r"[^a-z0-9]+", "_", re.sub(r"^[A-Z]?\d+\.\s*", "", s).lower()).strip("_")

genes = list(gene_cols) + ["Total"]

# Step 1: per-group summary (n, mean, median) for each gene, for ecological groups and for
# the clades I test. I save it as a table so the paper can quote it.
summary_rows = []
for grouping, col, levels in [("ecology", "group", group_order + ["Domesticated (all)"]),
                              ("clade", "clade_group", beer_clades + ["Wild clades"])]:
    for level in levels:
        mask = df["domesticated"] if level == "Domesticated (all)" else df[col] == level
        for gene in genes:
            x = df.loc[mask, gene]
            summary_rows.append({"grouping": grouping, "group": level, "gene": gene,
                                 "n": len(x), "mean": x.mean(), "median": x.median()})
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("data/results/mal_cnv_1011_summary.csv", index=False)
# I also keep the per-strain table (small: one row per strain) so the main figure can plot individual strains.
df.to_csv("data/results/mal_cnv_1011_per_strain.csv", index=False)
summary_df.pivot_table(index=["grouping", "group"], columns="gene",
                       values="mean", sort=False).round(2)

# %%
# Step 1: two-sided Mann-Whitney U tests. Copy number is a small, skewed, partly discrete
# count, so I use a rank test instead of a t-test. Each comparison is a family of 4 tests
# (3 transporter columns + total) and I correct within each family with Benjamini-Hochberg.
comparisons = {
    "beer_vs_wild": (df["group"] == "Beer", df["group"] == "Wild"),
    "domesticated_vs_wild": (df["domesticated"], df["group"] == "Wild"),
    "ale_beer_clade_vs_wild_clades": (df["clade"] == "11. Ale beer",
                                      df["clade"].isin(wild_clades)),
    "mosaic_beer_clade_vs_wild_clades": (df["clade"] == "7. Mosaic beer",
                                         df["clade"].isin(wild_clades)),
    "african_beer_clade_vs_wild_clades": (df["clade"] == "6. African beer",
                                          df["clade"].isin(wild_clades)),
}
test_rows = []
for name, (m1, m2) in comparisons.items():
    for gene in genes:
        x, y = df.loc[m1, gene], df.loc[m2, gene]
        test_rows.append({"comparison": name, "test": "Mann-Whitney U (two-sided)",
                          "gene": gene, "n1": len(x), "n2": len(y),
                          "median1": x.median(), "median2": y.median(),
                          "mean1": x.mean(), "mean2": y.mean(),
                          "statistic": mannwhitneyu(x, y, alternative="two-sided").statistic,
                          "p": mannwhitneyu(x, y, alternative="two-sided").pvalue,
                          "cliffs_delta": cliffs_delta(x, y)})

# Step 2: Kruskal-Wallis across the 9 ecological groups, the rank-based version of a one-way
# ANOVA, to see if ANY group differs (not just beer vs wild).
grouped = df.dropna(subset=["group"])
for gene in genes:
    samples = [grouped.loc[grouped["group"] == g, gene] for g in group_order]
    h = kruskal(*samples)
    test_rows.append({"comparison": "kruskal_ecology_groups", "test": "Kruskal-Wallis",
                      "gene": gene, "n1": len(grouped), "n2": np.nan,
                      "statistic": h.statistic, "p": h.pvalue, "cliffs_delta": np.nan})
tests_df = pd.DataFrame(test_rows)

# Step 3: BH correction within each comparison (4 genes per family).
tests_df["p_bh"] = tests_df.groupby("comparison")["p"].transform(
    lambda p: multipletests(p, method="fdr_bh")[1])
assert tests_df["p_bh"].ge(tests_df["p"] - 1e-12).all(), "BH should never lower a p-value"
tests_df.to_csv("data/results/mal_cnv_1011_tests.csv", index=False)
tests_df[["comparison", "gene", "n1", "n2", "median1", "median2", "p", "p_bh",
          "cliffs_delta"]].round(4)

# %%
# Step 1: bonus check on the flavour side. A strain with FDC1 (or the PAD1 column) at copy
# number 0 has lost the gene outright, which is another way to be POF- besides the point
# mutations studied in steps 3-5. ABSENT_CN = 0 because the matrix reports clean zeros
# (I checked the value counts: deletions are exactly 0.0, the next value is 0.5).
ABSENT_CN = 0
absent_rows = []
for level in group_order + ["Domesticated (all)"]:
    mask = df["domesticated"] if level == "Domesticated (all)" else df["group"] == level
    for gene in pof_cols:
        x = df.loc[mask, gene]
        absent_rows.append({"group": level, "gene": gene, "n": len(x),
                            "n_absent": int((x <= ABSENT_CN).sum()),
                            "frac_absent": (x <= ABSENT_CN).mean()})
absent_df = pd.DataFrame(absent_rows)
print("FDC1 absent, all strains:", int((df["FDC1"] <= ABSENT_CN).sum()), "of", len(df))
print("PAD1 (merged) absent, all strains:",
      int((df["PAD1 (merged)"] <= ABSENT_CN).sum()), "of", len(df))
print("Both absent:", int(((df["FDC1"] <= ABSENT_CN)
                           & (df["PAD1 (merged)"] <= ABSENT_CN)).sum()))
absent_df.pivot(index="group", columns="gene", values="frac_absent").loc[
    group_order + ["Domesticated (all)"]].round(3)

# %%
# Step 1: put every number the paper might quote into RESULTS as plain python types.
RESULTS["a4_n_strains_joined"] = int(len(df))
RESULTS["a4_n_strains_grouped"] = int(df["group"].notna().sum())
RESULTS["a4_n_domesticated"] = int(df["domesticated"].sum())
RESULTS["a4_wild_clade_frac_cutoff"] = WILD_CLADE_FRAC
RESULTS["a4_wild_clades"] = list(wild_clades)
RESULTS["a4_n_wild_clades_strains"] = int(df["clade"].isin(wild_clades).sum())
for g, n in group_n.items():
    RESULTS[f"a4_n_{slug(g)}"] = int(n)
for c in beer_clades:
    RESULTS[f"a4_n_clade_{slug(c)}"] = int((df["clade"] == c).sum())

# Step 2: summary means/medians and every test (p, BH p, Cliff's delta).
for r in summary_df.itertuples():
    key = f"a4_{r.grouping}_{slug(r.group)}_{slug(r.gene)}"
    RESULTS[key + "_mean"] = float(r.mean)
    RESULTS[key + "_median"] = float(r.median)
for r in tests_df.itertuples():
    key = f"a4_{r.comparison}_{slug(r.gene)}"
    RESULTS[key + "_p"] = float(r.p)
    RESULTS[key + "_p_bh"] = float(r.p_bh)
    if not np.isnan(r.cliffs_delta):
        RESULTS[key + "_cliffs_delta"] = float(r.cliffs_delta)
RESULTS["a4_min_p_bh_mannwhitney"] = float(tests_df.loc[
    tests_df["test"].str.startswith("Mann"), "p_bh"].min())

# Step 3: FDC1 / PAD1 absence counts and fractions.
for r in absent_df.itertuples():
    key = f"a4_fdc1_{'fdc1' if r.gene == 'FDC1' else 'pad1merged'}_absent_{slug(r.group)}"
    RESULTS[key + "_n"] = int(r.n_absent)
    RESULTS[key + "_frac"] = float(r.frac_absent)
RESULTS["a4_fdc1_fdc1_absent_all_n"] = int((df["FDC1"] <= ABSENT_CN).sum())
RESULTS["a4_fdc1_pad1merged_absent_all_n"] = int((df["PAD1 (merged)"] <= ABSENT_CN).sum())
print("a4_ keys in RESULTS:", sum(k.startswith("a4_") for k in RESULTS))

# %%
# Extended Data figure. Nature style: Arial 6-7 pt, 183 mm (double column) wide, bold lower-case
# panel letters, Okabe-Ito colours (colour-blind safe). I colour by role, not by group, so the
# colours carry one idea: beer (vermillion), other domesticated (orange), wild (bluish green),
# clinical (sky blue). Each point is one strain, jittered sideways so ties stay visible.
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica"],
                     "font.size": 6, "axes.labelsize": 7, "axes.titlesize": 7,
                     "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                     "ytick.major.width": 0.5, "pdf.fonttype": 42, "svg.fonttype": "none"})
okabe = {"beer": "#D55E00", "dom": "#E69F00", "wild": "#009E73", "clin": "#56B4E9"}
eco_colors = {g: okabe["dom"] for g in domesticated}
eco_colors.update({"Beer": okabe["beer"], "Wild": okabe["wild"], "Clinical": okabe["clin"]})
clade_colors = {c: okabe["beer"] for c in beer_clades}
clade_colors["Wild clades"] = okabe["wild"]
rng = np.random.default_rng(0)  # fixed seed so the jitter is the same on every run
MM = 1 / 25.4                   # inches per mm

def strip_box(ax, data, col, levels, colors, gene, labels=None):
    """Draw a box plot plus jittered strain points of `gene` per level of `col`, with n under
    each tick label."""
    vals = [data.loc[data[col] == lv, gene].values for lv in levels]
    ax.boxplot(vals, positions=range(len(levels)), widths=0.6, showfliers=False,
               medianprops={"color": "black", "linewidth": 0.8},
               boxprops={"linewidth": 0.5}, whiskerprops={"linewidth": 0.5},
               capprops={"linewidth": 0.5})
    for i, (lv, v) in enumerate(zip(levels, vals)):
        ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), v, s=1.5, color=colors[lv],
                   alpha=0.6, linewidths=0, zorder=3)
    names = labels or levels
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f"{nm}\nn={len(v)}" for nm, v in zip(names, vals)],
                       rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylabel(f"{gene} copy number")
    ax.spines[["top", "right"]].set_visible(False)

# Three rows of two panels, so the nine group labels per panel have room (183 x 210 mm).
fig = plt.figure(figsize=(183 * MM, 210 * MM))
gs = fig.add_gridspec(3, 2, hspace=0.75, wspace=0.25)
short = {g: ("Other ferm." if g == "Other fermentation" else g) for g in group_order}
gene_titles = {"MAL11": "MAL11 (AGT1)", "MAL31": "MAL31 (merged column)",
               "MPH2/3": "MPH2/MPH3 (merged column)", "Total": "Total transporter copies"}

# Panels a-d: each transporter column and the total, by ecological origin.
for i, gene in enumerate(genes):
    ax = fig.add_subplot(gs[i // 2, i % 2])
    strip_box(ax, grouped, "group", group_order, eco_colors, gene,
              labels=[short[g] for g in group_order])
    ax.set_title(gene_titles[gene], loc="left")
    ax.text(-0.12, 1.08, "abcd"[i], transform=ax.transAxes, fontsize=8, fontweight="bold")

# Panel e: total copies in the beer clades vs the pooled wild clades, the lineage-based view.
clade_levels = beer_clades + ["Wild clades"]
ax = fig.add_subplot(gs[2, 0])
strip_box(ax, df, "clade_group", clade_levels, clade_colors, "Total",
          labels=["Ale beer", "Mosaic beer", "African beer", "Wild clades"])
ax.set_title("Total transporter copies by clade", loc="left")
ax.text(-0.12, 1.08, "e", transform=ax.transAxes, fontsize=8, fontweight="bold")

# Panel f: fraction of strains with FDC1 deleted (copy number 0), by ecological origin.
ax = fig.add_subplot(gs[2, 1])
fdc1_frac = absent_df[absent_df["gene"] == "FDC1"].set_index("group").loc[group_order]
ax.bar(range(len(group_order)), fdc1_frac["frac_absent"], width=0.65,
       color=[eco_colors[g] for g in group_order])
ax.set_xticks(range(len(group_order)))
ax.set_xticklabels([f"{short[g]}\nn={n}" for g, n in zip(group_order, fdc1_frac["n"])],
                   rotation=45, ha="right", rotation_mode="anchor")
ax.set_ylabel("Fraction with FDC1 copy number 0")
ax.set_title("FDC1 deletion by ecological origin", loc="left")
ax.spines[["top", "right"]].set_visible(False)
ax.text(-0.12, 1.08, "f", transform=ax.transAxes, fontsize=8, fontweight="bold")

fig.savefig("figures/extended/ED_mal_cnv_1011.png", dpi=300, bbox_inches="tight")
fig.savefig("figures/extended/ED_mal_cnv_1011.pdf", bbox_inches="tight")
plt.show()
