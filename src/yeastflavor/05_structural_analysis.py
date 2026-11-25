# %%
# 5. Structural context of the PAD1 / FDC1 variants
import os, re, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, NeighborSearch
from Bio.SeqUtils import seq1

S = "data/structure"; R = "data/results"; FIG = "figures"
for d in (S, R, FIG):
    os.makedirs(d, exist_ok=True)

FDC1_PDB = f"{S}/4ZAC_A_CA_het.pdb"
PAD1_PDB = f"{S}/AF-P33751-F1_CA.pdb"
FEAT = f"{S}/uniprot_features.json"

def fetch_if_missing():
    """Only runs when the cached files are absent (needs internet + requests)."""
    import requests
    if not os.path.exists(FDC1_PDB):
        txt = requests.get("https://files.rcsb.org/download/4ZAC.pdb", timeout=60).text
        keep = [l for l in txt.splitlines()
                if l.startswith(("HEADER", "REMARK   2 RES"))
                or (l.startswith("HELIX ") and l[19] == "A") or (l.startswith("SHEET ") and l[21] == "A")
                or (l.startswith("HETATM") and l[21] == "A" and l[17:20] != "HOH")
                or (l.startswith("ATOM") and l[21] == "A" and l[12:16] == " CA " and l[16] in " A")]
        open(FDC1_PDB, "w").write("\n".join(keep) + "\nEND\n")
    if not os.path.exists(PAD1_PDB):
        url = requests.get("https://alphafold.ebi.ac.uk/api/prediction/P33751", timeout=60).json()[0]["pdbUrl"]
        txt = requests.get(url, timeout=60).text
        keep = [l for l in txt.splitlines() if l.startswith("ATOM") and l[12:16] == " CA "]
        open(PAD1_PDB, "w").write("\n".join(keep) + "\nEND\n")
    if not os.path.exists(FEAT):
        out = {}
        for acc in ("Q03034", "P33751"):
            j = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60).json()
            ss = {"Helix", "Beta strand", "Turn"}
            pf = requests.get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}", timeout=60).json()
            out[acc] = {"gene": j["genes"][0]["geneName"]["value"], "length": j["sequence"]["length"],
                        "sequence": j["sequence"]["value"],
                        "features": [{"type": f["type"], "start": f["location"]["start"]["value"],
                                      "end": f["location"]["end"]["value"], "description": f.get("description", ""),
                                      "ligand": f.get("ligand", {}).get("name", "")}
                                     for f in j["features"] if f["type"] not in ss | {"Chain"}],
                        "pfam": [[r["metadata"]["accession"], r["metadata"]["name"], fr["start"], fr["end"]]
                                 for r in pf["results"] for loc in r["proteins"][0]["entry_protein_locations"]
                                 for fr in loc["fragments"]],
                        "secondary_structure": [[f["type"], f["location"]["start"]["value"], f["location"]["end"]["value"]]
                                                for f in j["features"] if f["type"] in ss]}
        json.dump(out, open(FEAT, "w"), indent=1)

if not all(os.path.exists(p) for p in (FDC1_PDB, PAD1_PDB, FEAT)):
    fetch_if_missing()

# %%
# 5a. Load structures and map UniProt features onto chain A
feat = json.load(open(FEAT))
UNI = {"FDC1": feat["Q03034"], "PAD1": feat["P33751"]}
parser = PDBParser(QUIET=True)

def load_chain(path):
    chain = parser.get_structure(os.path.basename(path), path)[0]["A"]
    prot = {r.id[1]: r for r in chain if r.id[0] == " " and "CA" in r}
    het = [a for r in chain if r.id[0].startswith("H_") for a in r]
    return prot, het

STRUCT = {"FDC1": load_chain(FDC1_PDB), "PAD1": load_chain(PAD1_PDB)}
for g, (prot, het) in STRUCT.items():
    seq = UNI[g]["sequence"]
    mism = [n for n, r in prot.items() if seq1(r.get_resname()) != seq[n - 1]]
    print(f"{g}: {len(prot)} residues with CA (UniProt length {len(seq)}), "
          f"{len(seq) - len(prot)} missing, {len(mism)} sequence mismatches; HETATM groups: "
          f"{sorted({a.get_parent().get_resname() for a in het})}")

prot_f, het_f = STRUCT["FDC1"]
prfmn = [a for a in het_f if a.get_parent().get_resname() == "4LU"]
mn = [a for a in het_f if a.get_parent().get_resname() == "MN"]
ACTIVE = {"FDC1": {175: "R175", 280: "E280", 285: "E285 (proton donor)"}}

def binding_sites(g):
    """[(resnum, label)] for UniProt binding/active sites."""
    out = []
    for f in UNI[g]["features"]:
        if f["type"] in ("Binding site", "Active site"):
            lab = f["ligand"] or f["description"]
            out += [(i, lab) for i in range(f["start"], f["end"] + 1)]
    return out

BIND = {g: binding_sites(g) for g in UNI}
CA = {g: np.array([STRUCT[g][0][n]["CA"].coord for n in sorted(STRUCT[g][0])]) for g in UNI}
RESNUM = {g: np.array(sorted(STRUCT[g][0])) for g in UNI}
BFAC = {g: np.array([STRUCT[g][0][n]["CA"].bfactor for n in sorted(STRUCT[g][0])]) for g in UNI}

# prFMN contacts measured from the structure itself (CA within 8 Å of any prFMN atom)
ns_cof = NeighborSearch(prfmn + mn)
cof_contact = sorted(n for n, r in prot_f.items() if ns_cof.search(r["CA"].coord, 8.0))
print("Fdc1 residues with CA within 8 Å of prFMN/Mn:", cof_contact)
print("UniProt Fdc1 binding sites:", sorted({n for n, _ in BIND['FDC1']}))

# %%
# 5b. Variants: truncations (LOF alleles) and missense sites
AS, VC = f"{R}/allele_summary.csv", f"{R}/variant_calls.csv"
if os.path.exists(AS) and os.path.exists(VC):
    allele = pd.read_csv(AS); calls = pd.read_csv(VC); MOCK = False
else:  # shortcut: mock used only until section 3 has been run
    MOCK = True
    allele = pd.DataFrame({"gene": ["FDC1", "FDC1"], "strain": ["mock ale", "mock lager"],
                           "group": ["Ale", "Lager"], "protein_len": [153, 169], "partial": [False, False],
                           "LOF": [True, True], "LOF_allele": ["c.460C>T p.Q154*", "c.501dupA p.K167fs"]})
    calls = pd.DataFrame({"gene": ["FDC1", "FDC1", "PAD1"], "strain": ["mock"] * 3, "group": ["Wine"] * 3,
                          "effect": ["missense"] * 3, "protein_change": ["p.P249L", "p.S317L", "p.G220D"]})
print("MOCK variant data" if MOCK else "Using section-3 variant calls")

FULL = {g: UNI[g]["length"] for g in UNI}
first_pos = lambda s: min(int(m) for m in re.findall(r"p\.[A-Z](\d+)(?:\*|fs)", s))

ok = allele.callable.astype(bool) if "callable" in allele else ~allele.partial.astype(bool)  # B1 schema v2 / v1
trunc = allele[(allele.protein_len < allele.gene.map(FULL)) & ok].copy()
trunc["allele"] = trunc.apply(lambda r: r.LOF_allele if isinstance(r.LOF_allele, str) and r.LOF_allele
                              else "; ".join(calls[(calls.sacc == r.sacc) & (calls.gene == r.gene) &
                                                   calls.effect.isin(["nonsense", "frameshift"])].protein_change)
                              if "sacc" in calls else "", axis=1)
trunc["native_len"] = trunc.allele.map(first_pos) - 1
lof_rows = []
for (g, al), d in trunc.groupby(["gene", "allele"]):
    L, nat = FULL[g], int(d.native_len.iloc[0])
    lost_bind = sorted({f"{lab}:{n}" for n, lab in BIND[g] if n > nat})
    doms = [f"{name} {s}-{e} ({'lost' if s > nat else 'cut' if e > nat else 'kept'})"
            for _, name, s, e in UNI[g]["pfam"]]
    lof_rows.append({"gene": g, "allele": al, "native_residues": nat, "fraction_retained": round(nat / L, 3),
                     "LOF_called": bool(d.LOF.any()), "n_records": len(d),
                     "groups": "; ".join(sorted(set(d.group))), "strains": "; ".join(sorted(set(d.strain))),
                     "pfam_domains": "; ".join(doms),
                     "n_binding_residues_lost": len({n for n, _ in BIND[g] if n > nat}),
                     "n_binding_residues_total": len({n for n, _ in BIND[g]}),
                     "binding_lost": "; ".join(lost_bind)})
lof = pd.DataFrame(lof_rows).sort_values(["gene", "native_residues"])
lof["short"] = lof.allele.map(lambda a: re.findall(r"p\.[A-Z]\d+(?:\*|fs)", a)[0])
# one marker per truncation point for plotting (Q154* and Q154*+K167fs end at the same residue)
cuts = lof.groupby(["gene", "native_residues"]).agg(label=("short", lambda x: "/".join(dict.fromkeys(x))),
    LOF_called=("LOF_called", "any"), fraction_retained=("fraction_retained", "first")).reset_index()
lof.to_csv(f"{R}/structure_lof_context.csv", index=False)
print(lof[["gene", "allele", "native_residues", "fraction_retained", "n_binding_residues_lost",
           "n_binding_residues_total", "groups"]].to_string(index=False))

# %%
# Missense variants: distance to cofactor / active site, and local order (B-factor or pLDDT)
mis = calls[calls.effect == "missense"].copy()
mis[["ref_aa", "pos", "alt_aa"]] = mis.protein_change.str.extract(r"p\.([A-Z])(\d+)([A-Z])")
mis["pos"] = mis.pos.astype(int)
mis = (mis.groupby(["gene", "protein_change", "ref_aa", "pos", "alt_aa"])
          .agg(n_records=("strain", "size"), groups=("group", lambda x: "; ".join(sorted(set(x)))))
          .reset_index())

def ca_of(g, n):
    r = STRUCT[g][0].get(n)
    return None if r is None else r["CA"]

rows = []
for _, m in mis.iterrows():
    g, n = m.gene, m.pos
    assert UNI[g]["sequence"][n - 1] == m.ref_aa, f"{g} {m.protein_change}: ref aa mismatch"
    a = ca_of(g, n)
    row = m.to_dict()
    if a is None:
        rows.append(row); continue
    row["b_or_plddt"] = round(float(a.bfactor), 1)
    if g == "FDC1":
        row["d_prFMN_A"] = round(min(a - x for x in prfmn), 1)
        row["d_Mn_A"] = round(min(a - x for x in mn), 1)
        row["d_active_site_A"] = round(min(a - ca_of(g, k) for k in ACTIVE["FDC1"]), 1)
    else:
        sites = sorted({k for k, _ in BIND[g] if k != n})
        row["d_binding_site_CA_A"] = round(min(a - ca_of(g, k) for k in sites), 1)
    row["in_pfam_domain"] = "; ".join(name for _, name, s, e in UNI[g]["pfam"] if s <= n <= e) or "none"
    row["in_transit_peptide"] = g == "PAD1" and n <= 58
    rows.append(row)
mctx = pd.DataFrame(rows).sort_values(["gene", "pos"])
mctx.to_csv(f"{R}/structure_missense_context.csv", index=False)
print(mctx.drop(columns=["ref_aa", "alt_aa"]).to_string(index=False))

# %%
# 5c. Figures
# (a) feature tracks with truncation points and missense sites
fig, axes = plt.subplots(2, 1, figsize=(11, 5.2), gridspec_kw={"height_ratios": [1, 1]})
dom_col = ["#9ecae1", "#a1d99b", "#fdd0a2"]
for ax, g in zip(axes, ["PAD1", "FDC1"]):
    L = FULL[g]
    ax.add_patch(plt.Rectangle((1, -0.15), L, 0.3, color="#d9d9d9"))
    for i, (_, name, s, e) in enumerate(UNI[g]["pfam"]):
        ax.add_patch(plt.Rectangle((s, -0.15), e - s + 1, 0.3, color=dom_col[i % 3], label=name))
    for f in UNI[g]["features"]:
        if f["type"] == "Transit peptide":
            ax.add_patch(plt.Rectangle((f["start"], -0.15), f["end"] - f["start"] + 1, 0.3,
                                       fill=False, hatch="//", ec="grey", label="mito. transit peptide"))
    bpos = sorted({n for n, _ in BIND[g]})
    ax.plot(bpos, [0.22] * len(bpos), "v", color="#6a3d9a", ms=5, label="UniProt binding/active site")
    for k, lab in ACTIVE.get(g, {}).items():
        ax.annotate(lab.split()[0], (k, 0.22 if k != 280 else 0.36), xytext=(0, 6), textcoords="offset points",
                    ha="center", fontsize=7, color="#6a3d9a")
    mg = mctx[mctx.gene == g]
    ax.plot(mg.pos, [-0.24] * len(mg), "^", color="#1f78b4", ms=5, label="missense (any strain)")
    for i, (_, t) in enumerate(cuts[cuts.gene == g].iterrows()):
        x = t.native_residues + 0.5
        ax.axvline(x, color="#e31a1c" if t.LOF_called else "#fb9a99", lw=1.5)
        ax.text(x, 0.5 + 0.17 * (i % 2), f"{t.label} ({t.native_residues} aa)", color="#e31a1c",
                fontsize=6.5, ha="right" if x > 0.9 * L else "center", va="bottom")
    ax.set_xlim(0, L + 10); ax.set_ylim(-0.5, 0.95); ax.set_yticks([])
    ax.set_title(f"{g} ({L} aa)", loc="left", fontsize=10)
    ax.legend(fontsize=6.5, loc="lower right", ncol=3, frameon=False, bbox_to_anchor=(1, -0.62))
axes[-1].set_xlabel("residue")
fig.suptitle("Pad1 and Fdc1: Pfam domains, UniProt sites, truncation points (red) and missense sites", fontsize=10)
fig.tight_layout(); fig.savefig(f"{FIG}/05a_feature_tracks.png", dpi=200); plt.show()

# %%
# (b) PCA projection of the Fdc1 CA trace, retained vs lost for each Fdc1 truncation
X = CA["FDC1"]; mu = X.mean(0)
_, _, Vt = np.linalg.svd(X - mu, full_matrices=False)
P = (X - mu) @ Vt[:2].T
Pc = (np.array([a.coord for a in prfmn]) - mu) @ Vt[:2].T
Pm = (np.array([a.coord for a in mn]) - mu) @ Vt[:2].T
ft = cuts[(cuts.gene == "FDC1") & cuts.LOF_called].reset_index(drop=True)
nc = min(3, len(ft)); nr = -(-len(ft) // nc)
fig, axes = plt.subplots(nr, nc, figsize=(4.2 * nc, 4.3 * nr), squeeze=False)
for ax in axes.flat[len(ft):]:
    ax.axis("off")
for ax, (_, t) in zip(axes.flat, ft.iterrows()):
    keep = RESNUM["FDC1"] <= t.native_residues
    ax.plot(P[:, 0], P[:, 1], color="#d9d9d9", lw=0.6, zorder=1)
    ax.scatter(*P[~keep].T, s=8, c="#bdbdbd", label="lost", zorder=2)
    ax.scatter(*P[keep].T, s=8, c="#e31a1c", label=f"retained (1-{t.native_residues})", zorder=3)
    ax.scatter(*Pc.T, s=14, c="#ff7f00", marker="s", label="prFMN", zorder=4)
    ax.scatter(*Pm.T, s=40, c="#6a3d9a", marker="*", label="Mn2+", zorder=4)
    idx = np.isin(RESNUM["FDC1"], mctx[mctx.gene == "FDC1"].pos)
    ax.scatter(*P[idx].T, s=40, facecolors="none", edgecolors="#1f78b4", lw=1.2, label="missense", zorder=5)
    ax.set_title(f"{t.label}  ({t.fraction_retained:.0%} retained)", fontsize=9)
    ax.set_xlabel("PC1 (Å)"); ax.set_ylabel("PC2 (Å)"); ax.set_aspect("equal")
axes[0][0].legend(fontsize=7, loc="best")
fig.suptitle("Fdc1 (PDB 4ZAC chain A) CA trace projected on its first two principal axes", fontsize=10)
fig.tight_layout(); fig.savefig(f"{FIG}/05b_fdc1_ca_pca_truncations.png", dpi=200); plt.show()

# %%
# (c) CA-CA contact map (< 8 Å) with truncation lines
D = np.linalg.norm(CA["FDC1"][:, None] - CA["FDC1"][None], axis=-1)
C = D < 8.0
fig, ax = plt.subplots(figsize=(6.2, 6))
ax.imshow(C, cmap="Greys", origin="lower", extent=[0.5, FULL["FDC1"] + 0.5] * 2, interpolation="nearest")
for _, t in ft.iterrows():
    x = t.native_residues + 0.5
    ax.axvline(x, color="#e31a1c", lw=0.8); ax.axhline(x, color="#e31a1c", lw=0.8)
    ax.text(x, FULL["FDC1"] + 4, t.label, color="#e31a1c", fontsize=7, rotation=90, va="bottom")
ax.set_xlabel("residue"); ax.set_ylabel("residue")
ax.set_title("Fdc1 CA-CA contacts < 8 Å (4ZAC chain A); red = truncation points", fontsize=9, pad=30)
fig.tight_layout(); fig.savefig(f"{FIG}/05c_fdc1_contact_map.png", dpi=200); plt.show()

# fraction of each truncation's long-range contacts (|i-j| > 5) that cross the cut
lr = C & (np.abs(RESNUM["FDC1"][:, None] - RESNUM["FDC1"][None]) > 5)
for _, t in ft.iterrows():
    k = RESNUM["FDC1"] <= t.native_residues
    kept_lr = lr[np.ix_(k, k)].sum() / 2; cross = lr[np.ix_(k, ~k)].sum()
    print(f"{t.label}: retained segment has {int(kept_lr)} internal long-range contacts and "
          f"{int(cross)} contacts to the lost part ({cross / (kept_lr + cross):.0%} of its contacts)")

# %%
# (d) distance to cofactor for missense variants
fm = mctx[mctx.gene == "FDC1"].sort_values("d_prFMN_A", ascending=False)
pm = mctx[mctx.gene == "PAD1"].sort_values("d_binding_site_CA_A", ascending=False)
lab = lambda d: d.protein_change.str[2:] + " (n=" + d.n_records.astype(str) + ")"
fig, axes = plt.subplots(1, 2, figsize=(11, 0.2 * max(len(fm), len(pm)) + 1.5))
ax = axes[0]
ax.barh(lab(fm), fm.d_prFMN_A, color="#1f78b4", label="CA to nearest prFMN atom")
ax.scatter(fm.d_active_site_A, lab(fm), color="#e31a1c", s=12, zorder=3, label="CA to nearest R175/E280/E285 CA")
ax.axvline(8, ls="--", color="grey", lw=0.8)
ax.set_xlabel("distance (Å)"); ax.set_title("Fdc1 missense (PDB 4ZAC)", fontsize=9)
ax.legend(fontsize=7, loc="upper right"); ax.tick_params(axis="y", labelsize=7); ax.margins(y=0.01)
ax = axes[1]
cols = ["#fdbf6f" if t else "#33a02c" for t in pm.in_transit_peptide]
ax.barh(lab(pm), pm.d_binding_site_CA_A, color=cols)
for y, (d, q) in enumerate(zip(pm.d_binding_site_CA_A, pm.b_or_plddt)):
    ax.text(d + 0.5, y, f"{q:.0f}", va="center", fontsize=6)
ax.axvline(8, ls="--", color="grey", lw=0.8)
ax.set_title("Pad1 missense (AlphaFold model)\nnumber = pLDDT; orange = mito. transit peptide", fontsize=9)
ax.set_xlabel("CA to nearest FMN/DMAP-site CA (Å)"); ax.tick_params(axis="y", labelsize=7); ax.margins(y=0.01)
fig.tight_layout(); fig.savefig(f"{FIG}/05d_missense_distances.png", dpi=200); plt.show()
print("figures written:", sorted(f for f in os.listdir(FIG) if f.startswith("05")))
