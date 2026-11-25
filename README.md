# Clean Beer, Broken Genes

**Tracing how brewing yeasts lost the clove flavor**

<!-- DOI badge goes here once Zenodo mints it: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

Code, derived data and manuscript for a BCH 339N Bioinformatics (UT Austin, Fall 2026) group project on how
domesticated *Saccharomyces cerevisiae* lost the clove-like phenolic off-flavor (4-vinylguaiacol) through
loss-of-function mutations in *PAD1* and *FDC1*.

**Main findings** (every number is produced by the code and stored in `results/results.json`):
- Predicted POF− genotypes: 77 of 126 domesticated beverage strains vs 2 of 35 wild, wine and bioethanol strains; the enrichment holds within population clades (Cochran–Mantel–Haenszel test).
- 8 distinct loss-of-function mutations that add up to at least 4 separate losses (ale/Beer 1, saké/shochu, single wine and lab strains); Beer 2 may add one more through a late stop codon.
- Every common *FDC1* allele stops before all prFMN-, Mn²⁺- and active-site residues of Fdc1 (PDB 4ZAC).
- Genotype predicts the laboratory phenotype of Gallone *et al.* 2016 with 77% sensitivity and 86% specificity.
- Maltose-transporter copy number (1,011 genomes) shows gene-level trade-offs, not a net gain; the two lager assemblies carry no *S. eubayanus* copy of *PAD1*/*FDC1*.

The paper: [`paper/BCH339N_Sperling_brewing_yeast_flavor.pdf`](paper/BCH339N_Sperling_brewing_yeast_flavor.pdf).

## Reproduce

```bash
python -m venv .venv && .venv/Scripts/activate      # or: conda env create -f environment.yml
pip install -r requirements.txt
export NCBI_EMAIL=you@example.org                   # NCBI asks for a contact address; never hard-coded here
python run_all.py                                   # data -> analysis -> figures -> results/results.json (~4 min cached)
python -m pytest tests                              # unit checks + results consistency
python paper/build_paper.py                         # manuscript PDF (needs pandoc and Edge or Chrome)
```

Every acquisition step is **cache-first**: it reads the file in `data/` if present and otherwise fetches it from
NCBI / Europe PMC / the 1002 Genomes site, checking the SHA-256 listed in `data/MANIFEST.csv`. Large published
supplementary files live in `data/raw/` (git-ignored, re-downloaded on demand).

## Layout

| Path | Contents |
|---|---|
| `run_all.py` | runs every step in order (one shared namespace), writes `notebooks/*.ipynb` and `results/results.json` |
| `src/yeastflavor/NN_*.py` | analysis steps, notebook-style cells (`# %%`); each starts with the lecture its technique comes from |
| `data/` | cached inputs, `MANIFEST.csv` (source, URL/accession, date, checksum), `results/` derived tables |
| `figures/main/` | Figures 1–6 (PDF + 300-dpi PNG); `figures/` and `figures/extended/` hold Extended Data |
| `paper/` | `manuscript.md` (numbers as `{{tokens}}` from `results.json`), `references.bib` (all PubMed/DOI-verified), `build_paper.py` |
| `tests/` | LOF caller, phylogeny and structure self-checks; `results.json` cross-checks |

## Figure → script map

| Figure | Script |
|---|---|
| Fig. 1–6 | `src/yeastflavor/12_main_figures.py` |
| Extended Data 1–2 (variant tracks, LOF alignment windows) | `03_variant_calling.py` |
| Extended Data 3–5, 9 (identity heatmap, trees, PCA/k-means, assembly CN) | `04_phylogeny_clustering_cnv.py` |
| Extended Data 6–7 (contact map, missense distances) | `05_structural_analysis.py` |
| Extended Data 8 (1,011-genome CN) | `09_cnv_1011_genomes.py` |
| lager subgenome grid | `10_lager_eubayanus_subgenome.py` |

## Data sources

S288C reference (NCBI RefSeq); NCBI core_nt BLAST hits (*S. cerevisiae* txid4932, *S. pastorianus* txid27292);
Gallone *et al.* 2016 genomes (BioProject PRJNA323691) and supplementary Tables S1/S5; 1,011 Yeast Genomes copy-number
matrix and Table S1 (Peter *et al.* 2018); PDB 4ZAC; AlphaFold AF-P33751-F1; UniProt Q03034/P33751; *S. eubayanus*
FM1318 (GCF_001298625.1); lager assemblies GCA_011022315.1 and GCA_056824455.1. Full list with checksums in `data/MANIFEST.csv`.

## License and citation

Code: MIT (`LICENSE`). Generated data tables and figures: CC-BY-4.0. Please cite via `CITATION.cff`.
Generative-AI use is described in the paper's Methods (Claude, Anthropic).
