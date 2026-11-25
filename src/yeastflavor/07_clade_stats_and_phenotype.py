# Technique source: Lecture 19 (enrichment / hypergeometric, scipy.stats) — notebook_notes/lecture_19.md;
# Lecture 22 (statsmodels) — notebook_notes/lecture_22.md. Fisher / CMH / Wilson CIs are beyond class.

# %%
# 7. Is POF loss a domestication signal once relatedness is accounted for? (A1) Does the genotype predict the measured phenotype? (A3)
import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.proportion import proportion_confint

RES, EXT = "data/results", "data/raw/external/gallone2016"
for _f in ("mmc1.xlsx", "mmc5.xlsx"):
    ensure_external(f"{EXT}/{_f}")      # cache-first download + checksum (step 00)
GENO_CSV = f"{RES}/pof_genotype_by_strain.csv"
assert os.path.exists(GENO_CSV), "run steps 01-03 first: pof_genotype_by_strain.csv is missing"

# Step 1: load the Gallone strain table. The first rows are section headers ("Strains sequenced in this
# study"), so I keep only rows whose Code looks like two letters + three digits (BE001, WI012, ...).
s1 = pd.read_excel(f"{EXT}/mmc1.xlsx")
s1 = s1[s1.Code.astype(str).str.fullmatch(r"[A-Z]{2}\d{3}")].copy()

# The Gallone tables use codes (BE001) but the NCBI BioSample records, and therefore our genotype table,
# use names like beer001. The prefix mapping below was checked one-to-one against the 157 BioSamples
#, so I can join on it safely.
PREFIX = {"BE": "beer", "WI": "wine", "SP": "spirits", "SA": "sake", "WL": "wild", "BI": "bioethanol",
          "BR": "bread", "LA": "laboratory"}
s1["strain"] = s1.Code.str[:2].map(PREFIX) + s1.Code.str[2:]
s1["clade"] = s1.Subpopulation.astype(str).str.strip().replace({"nan": np.nan, "NS": np.nan})
assert s1.strain.notna().all(), "an unexpected strain-code prefix appeared in Table S1"
print(f"Gallone strains in Table S1: {len(s1)}; with a clade label: {s1.clade.notna().sum()}")
s1.clade.value_counts()

# %%
# Step 2: join the clade labels onto our predicted genotypes. Only Gallone strains have clade labels, so
# everything in this section is restricted to that panel (one genome per strain, one sequencing study).
geno = pd.read_csv(GENO_CSV)
gal = s1.merge(geno, on="strain", how="left")
gal["called"] = gal.predicted_POF.isin(["POF-", "POF+"])
print(f"Gallone strains matched to a genotype row: {gal.predicted_POF.notna().sum()} / {len(gal)}")
print(f"... with a POF call (both genes callable, or one LOF): {gal.called.sum()}")

# I group strains the same way as section 3.6 so the numbers are comparable: domesticated beverage =
# beer, sake, spirits; reference = wine, wild, bioethanol. Bread and laboratory strains are left out of the
# 2x2 (same as before) because they fit neither side cleanly. I use the strain-name prefix, not the Niche
# column, because Niche has sub-labels ("Beer (Lager)", "Saké", "Wild (Cocoa)") that an exact match misses.
DOMESTIC_PREFIX, REFERENCE_PREFIX = ["beer", "sake", "spirits"], ["wine", "wild", "bioethanol"]
gal["prefix"] = gal.strain.str.extract(r"^([a-z]+)")[0]
gal["side"] = np.where(gal.prefix.isin(DOMESTIC_PREFIX), "domesticated",
                       np.where(gal.prefix.isin(REFERENCE_PREFIX), "reference", "excluded"))
gal["pof_neg"] = (gal.predicted_POF == "POF-").astype(int)
G = gal[gal.called & (gal.side != "excluded")]
pd.crosstab([G.clade.fillna("no label"), G.side], G.predicted_POF)

# %%
# 7.1 Strain-level Fisher test (Gallone panel only) and the clade-stratified CMH test
# Step 1: the uncorrected strain-level test on this panel, for comparison with section 3.6.
tab = [[((G.side == "domesticated") & (G.pof_neg == 1)).sum(), ((G.side == "domesticated") & (G.pof_neg == 0)).sum()],
       [((G.side == "reference") & (G.pof_neg == 1)).sum(), ((G.side == "reference") & (G.pof_neg == 0)).sum()]]
or_strain, p_strain = fisher_exact(tab)
print(f"strain-level 2x2 [[dom POF-, dom POF+], [ref POF-, ref POF+]]: {tab}")
print(f"Fisher OR: {or_strain:.2f}   p: {p_strain:.2e}")

# Step 2: build one 2x2 per clade and keep the informative ones (both sides present). I add 0.5 to empty
# cells only for the per-clade OR display; the CMH statistic itself uses the raw counts.
strata, rows = [], []
for clade, d in G.dropna(subset=["clade"]).groupby("clade"):
    t = np.array([[((d.side == "domesticated") & (d.pof_neg == 1)).sum(), ((d.side == "domesticated") & (d.pof_neg == 0)).sum()],
                  [((d.side == "reference") & (d.pof_neg == 1)).sum(), ((d.side == "reference") & (d.pof_neg == 0)).sum()]])
    informative = t[0].sum() > 0 and t[1].sum() > 0
    rows.append(dict(clade=clade, dom_POFneg=t[0, 0], dom_POFpos=t[0, 1], ref_POFneg=t[1, 0], ref_POFpos=t[1, 1],
                     informative=informative))
    if informative:
        strata.append(t)
cmh_tbl = pd.DataFrame(rows)
assert len(strata) >= 2, "CMH needs at least two clades that contain both domesticated and reference strains"

# Step 3: pooled Mantel-Haenszel OR with 95% CI, the CMH test of OR = 1, and Breslow-Day homogeneity.
st = StratifiedTable([t.astype(float) for t in strata])
cmh = st.test_null_odds(correction=True)
bd = st.test_equal_odds()
ci_lo, ci_hi = st.oddsratio_pooled_confint()
print(f"informative clades: {len(strata)} of {len(cmh_tbl)} ({', '.join(cmh_tbl[cmh_tbl.informative].clade)})")
print(f"strains in informative clades: {int(sum(t.sum() for t in strata))}")
print(f"MH pooled OR: {st.oddsratio_pooled:.2f} (95% CI {ci_lo:.2f}-{ci_hi:.2f}); CMH chi2: {cmh.statistic:.2f}, p: {cmh.pvalue:.2e}")
print(f"Breslow-Day homogeneity p: {bd.pvalue:.3f}")
# Step 3b: how fragile is the pooled OR? A clade only feeds the MH estimate if a*d + b*c > 0 (Britain, all
# POF-, adds nothing), and I count in how many of those the domesticated side is more often POF-.
contrib = [t for t in strata if t[0, 0] * t[1, 1] + t[0, 1] * t[1, 0] > 0]
dom_higher = sum(t[0, 0] / t[0].sum() > t[1, 0] / t[1].sum() for t in contrib)
print(f"clades that feed the MH estimate: {len(contrib)}; domesticated more often POF- in {dom_higher} of them")
cmh_tbl.to_csv(f"{RES}/cmh_strata_by_clade.csv", index=False)
cmh_tbl

# %%
# Step 4: clade-level view. I collapse to one row per clade so a clonal expansion counts once. A clade is
# called "POF− clade" when most of its called strains are POF−; that threshold is my own simple choice.
MAJORITY = 0.5
clade_view = (gal[gal.called].dropna(subset=["clade"]).groupby("clade")
              .agg(n=("strain", "size"), n_POFneg=("pof_neg", "sum"),
                   niches=("Niche", lambda x: ", ".join(f"{k} {v}" for k, v in x.value_counts().items()))))
clade_view["frac_POFneg"] = clade_view.n_POFneg / clade_view.n
clade_view["POFneg_clade"] = clade_view.frac_POFneg > MAJORITY
clade_view = clade_view.sort_values("frac_POFneg", ascending=False)
clade_view.to_csv(f"{RES}/pof_by_gallone_clade.csv")
print(f"clades that are mostly POF-: {clade_view.POFneg_clade.sum()} of {len(clade_view)}")
clade_view.round(2)

# %%
# Step 5: save the A1 numbers for the paper.
RESULTS.update({
    "a1_gallone_n_called": int(len(G)), "a1_gallone_table": [[int(x) for x in r] for r in tab],
    "a1_gallone_fisher_or": float(or_strain), "a1_gallone_fisher_p": float(p_strain),
    "a1_cmh_n_strata": len(strata), "a1_cmh_strata": list(cmh_tbl[cmh_tbl.informative].clade),
    "a1_cmh_n_strains": int(sum(t.sum() for t in strata)),
    "a1_cmh_n_contributing": len(contrib), "a1_cmh_n_dom_higher": int(dom_higher),
    "a1_cmh_or": float(st.oddsratio_pooled), "a1_cmh_or_ci": [float(ci_lo), float(ci_hi)],
    "a1_cmh_chi2": float(cmh.statistic), "a1_cmh_p": float(cmh.pvalue), "a1_breslow_day_p": float(bd.pvalue),
    "a1_n_clades": int(len(clade_view)), "a1_n_pofneg_clades": int(clade_view.POFneg_clade.sum()),
    "a1_pofneg_clades": list(clade_view[clade_view.POFneg_clade].index),
    "a1_clade_frac_pofneg": {k: float(v) for k, v in clade_view.frac_POFneg.items()},
    "a1_clade_n": {k: int(v) for k, v in clade_view.n.items()},
})

# %%
# 7.2 Does the predicted genotype match the measured phenotype? (A3)
# Step 1: read the phenotype and turn the z-score into a label.
s5 = pd.read_excel(f"{EXT}/mmc5.xlsx")
s5 = s5[s5["Sample Name"].astype(str).str.fullmatch(r"[A-Z]{2}\d{3}")].copy()
s5["strain"] = s5["Sample Name"].str[:2].map(PREFIX) + s5["Sample Name"].str[2:]
pof_num = pd.to_numeric(s5.POF, errors="coerce")
s5["measured_POF"] = np.select([np.isclose(pof_num, -0.921882, atol=1e-4), np.isclose(pof_num, 1.077693, atol=1e-4)],
                               ["POF-", "POF+"], default="ambiguous")
print(f"phenotyped strains: {len(s5)}; label counts: {s5.measured_POF.value_counts().to_dict()}")

# Step 2: join to the genotype calls and keep strains with both a clear phenotype and a POF call.
ph = s5.merge(gal[["strain", "Niche", "clade", "PAD1", "FDC1", "predicted_POF"]], on="strain", how="left")
both = ph[ph.measured_POF.isin(["POF-", "POF+"]) & ph.predicted_POF.isin(["POF-", "POF+"])].copy()
print(f"strains with both a clear phenotype and a genotype call: {len(both)}")
print(f"phenotyped but no genotype call ('?'): {(ph.measured_POF.isin(['POF-', 'POF+']) & ~ph.predicted_POF.isin(['POF-', 'POF+'])).sum()}")

# Step 3: confusion matrix and the usual diagnostic numbers, each with a Wilson 95% CI (better than the
# normal approximation when a proportion is close to 1, which it is here).
cm = pd.crosstab(both.predicted_POF, both.measured_POF).reindex(index=["POF-", "POF+"], columns=["POF-", "POF+"], fill_value=0)
tp, fp = cm.loc["POF-", "POF-"], cm.loc["POF-", "POF+"]
fn, tn = cm.loc["POF+", "POF-"], cm.loc["POF+", "POF+"]

def rate(k, n):
    """Proportion k/n with its Wilson 95% interval, as (p, lo, hi)."""
    lo, hi = proportion_confint(k, n, method="wilson")
    return float(k / n), float(lo), float(hi)

diag = {"sensitivity": rate(tp, tp + fn), "specificity": rate(tn, tn + fp),
        "ppv": rate(tp, tp + fp), "npv": rate(tn, tn + fn), "accuracy": rate(tp + tn, len(both))}
po = (tp + tn) / len(both)
pe = ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / len(both) ** 2
kappa = (po - pe) / (1 - pe)          # Cohen's kappa: agreement beyond chance
for k, (p, lo, hi) in diag.items():
    print(f"{k}: {p:.1%} (95% CI {lo:.1%}-{hi:.1%})")
print(f"Cohen's kappa: {kappa:.2f}")
cm.to_csv(f"{RES}/pof_genotype_vs_phenotype_confusion.csv")
cm

# %%
# Step 4: the discordant strains are the interesting ones, so I list them with their genotype detail.
# Predicted POF+ but measured POF− means a loss our caller missed (another mutation, a missing record, or
# regulation); predicted POF− but measured POF+ points to heterozygosity (a working copy hidden behind a
# consensus LOF call) or to a mislabelled record.
disc = both[both.predicted_POF != both.measured_POF].sort_values(["measured_POF", "strain"])
disc[["strain", "Niche", "clade", "PAD1", "FDC1", "predicted_POF", "measured_POF"]].to_csv(f"{RES}/pof_discordant_strains.csv", index=False)
disc[["strain", "Niche", "clade", "PAD1", "FDC1", "predicted_POF", "measured_POF"]]

# %%
# Step 5: save the A3 numbers.
RESULTS.update({
    "a3_n_phenotyped": int(len(s5)), "a3_n_ambiguous_pheno": int((s5.measured_POF == "ambiguous").sum()),
    "a3_n_compared": int(len(both)),
    "a3_confusion": {"TP": int(tp), "FP": int(fp), "FN": int(fn), "TN": int(tn)},
    **{f"a3_{k}": v for k, v in diag.items()}, "a3_kappa": float(kappa),
    "a3_discordant": disc[["strain", "Niche", "clade", "predicted_POF", "measured_POF"]].to_dict("records"),
})

# %%
# 7.3 Which alleles explain the disagreements?
# Step 1: one "allele signature" per strain = all nonsense/frameshift alleles it carries, sorted.
vc = pd.read_csv(f"{RES}/variant_calls.csv")
trunc = vc[vc.effect.isin(["nonsense", "frameshift"])]
sig = (trunc.assign(a=trunc.gene + " " + trunc.protein_change.str.replace("p.", "", regex=False))
       .groupby("strain").a.apply(lambda x: " + ".join(sorted(set(x)))))
ph2 = ph[ph.measured_POF.isin(["POF-", "POF+"])].copy()
ph2["alleles"] = ph2.strain.map(sig).fillna("no truncating allele")
allele_pheno = pd.crosstab(ph2.alleles, ph2.measured_POF).reindex(columns=["POF-", "POF+"], fill_value=0)
allele_pheno["n"] = allele_pheno.sum(axis=1)
allele_pheno = allele_pheno.sort_values("n", ascending=False)
allele_pheno.to_csv(f"{RES}/pof_phenotype_by_allele_signature.csv")
allele_pheno

# %%
# Step 2: sensitivity analysis. Rule A = the caller as is. Rule B counts the late FDC1 stops (W497*, Y500*)
# as LOF, which drops the 95%-length rule. Rule C additionally ignores c.501dupA when it is the strain's only
# truncating allele, because it sits in an A6 homopolymer where assembly errors are common.
LATE = {"FDC1 W497*", "FDC1 Y500*"}
def rescore(row, count_late, drop_lone_dup):
    """Return the predicted POF label for one strain under the chosen rule."""
    al = set() if row.alleles == "no truncating allele" else set(row.alleles.split(" + "))
    if drop_lone_dup and al == {"FDC1 K167fs"}:
        al = set()
    if not count_late:
        al -= LATE
    if al:
        return "POF-"
    # no truncating allele left: intact only if both genes had a callable record (PAD1/FDC1 are not NaN)
    return "POF+" if pd.notna(row.PAD1) and pd.notna(row.FDC1) else "?"

sens_rows = []
for name, late, dup in [("A: caller as published here", False, False), ("B: late stops count as LOF", True, False),
                        ("C: B + ignore lone c.501dupA", True, True)]:
    pred = ph2.apply(rescore, axis=1, count_late=late, drop_lone_dup=dup)
    ok = pred.isin(["POF-", "POF+"])
    m, q = ph2.measured_POF[ok], pred[ok]
    tp_, fn_ = ((q == "POF-") & (m == "POF-")).sum(), ((q == "POF+") & (m == "POF-")).sum()
    tn_, fp_ = ((q == "POF+") & (m == "POF+")).sum(), ((q == "POF-") & (m == "POF+")).sum()
    sens_rows.append(dict(rule=name, n=int(ok.sum()), TP=int(tp_), FN=int(fn_), TN=int(tn_), FP=int(fp_),
                          sensitivity=tp_ / (tp_ + fn_), specificity=tn_ / (tn_ + fp_),
                          accuracy=(tp_ + tn_) / ok.sum()))
sens = pd.DataFrame(sens_rows)
assert sens.loc[0, ["TP", "FN", "TN", "FP"]].tolist() == [tp, fn, tn, fp] or sens.loc[0, "n"] != len(both),     "rule A should reproduce the confusion matrix above"
sens.to_csv(f"{RES}/pof_caller_sensitivity.csv", index=False)
sens.round(3)

# %%
# Step 3: save the allele-level numbers the paper needs.
def counts_for(label):
    """Measured POF-/POF+ counts for one allele signature (0, 0 if absent)."""
    if label not in allele_pheno.index:
        return {"POF-": 0, "POF+": 0}
    return {k: int(allele_pheno.loc[label, k]) for k in ["POF-", "POF+"]}
RESULTS.update({
    "a3_allele_pheno": {k: {c: int(v) for c, v in r.items()} for k, r in allele_pheno.iterrows()},
    "a3_lone_dupA": counts_for("FDC1 K167fs"),
    "a3_w497": counts_for("FDC1 W497*"),
    "a3_sensitivity_rules": sens.to_dict("records"),
})
