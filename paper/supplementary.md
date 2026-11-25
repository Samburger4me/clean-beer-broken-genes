
::: supp
# 7. Supplementary Material

Raw downloaded data are not reproduced here; every source is listed with URL, date and checksum in `data/MANIFEST.csv` and is fetched again by `run_all.py`. All derived tables are in `data/results/` in the code repository. Extended Data figures are produced by the pipeline steps named in each legend.

## Supplementary Tables

**Supplementary Table 1 | Strain genotypes.** Predicted *PAD1*/*FDC1* status and POF call for all {{R['core_n_strains']}} strains: `data/results/pof_genotype_by_strain.csv`; every variant call: `data/results/variant_calls.csv`; per-record allele summary: `data/results/allele_summary.csv`.

**Supplementary Table 2 | Distinct truncating alleles and their origins.** Strains, haplotypes, parsimony gains and Gallone clades per allele: `data/results/lof_independent_origins.csv`.

{{origins_table()}}

**Supplementary Table 3 | Sensitivity of the genotype–phenotype agreement to the two judgement calls in the LOF caller** (n = {{R['a3_n_compared']}} strains with a Gallone *et al.* phenotype and a genotype call; `data/results/pof_caller_sensitivity.csv`).

{{sens_table()}}

**Supplementary Table 4 | Clade strata used in the CMH test** (Gallone *et al.* clades; dom., beer/saké/spirits; ref., wine/wild/bioethanol; `data/results/cmh_strata_by_clade.csv`).

{{cmh_table()}}

**Supplementary Table 5 | Genotype–phenotype discordant strains** (`data/results/pof_discordant_strains.csv`) and the full allele-signature × phenotype table (`data/results/pof_phenotype_by_allele_signature.csv`). **Supplementary Table 6 | Copy-number summaries and tests** for the 1,011-genome panel (`data/results/mal_cnv_1011_summary.csv`, `mal_cnv_1011_tests.csv`). **Supplementary Table 7 | Lager subgenome calls** (`data/results/lager_subgenome_pad1_fdc1.csv`).

## Extended Data figures

::: figure
![](figures/03a_variant_tracks.png){width=92%}

**Extended Data Fig. 1 | Recurrent and LOF variants along *PAD1* and *FDC1* by ecological group.** Marker area is the fraction of the group carrying the variant (pipeline step 03).
:::

::: figure
![](figures/03b_lof_alignment_windows.png){width=92%}

**Extended Data Fig. 2 | Alignment windows around each LOF allele.** Nucleotide alignment of representative strains to S288C around each nonsense or frameshift site, showing the inserted A of c.501dupA in its A~6~ run and each stop codon (step 03).
:::

::: figure
![](figures/heatmap_identity_PAD1_FDC1.png){width=70%}

**Extended Data Fig. 3 | Pairwise identity of concatenated *PAD1*+*FDC1* sequences.** Hierarchically clustered heatmap (seaborn `clustermap`, as in Lecture 14) of percent identity between all strains (step 04).
:::

::: figure
![](figures/tree_PAD1_FDC1.png){width=92%}

**Extended Data Fig. 4 | Strain-level NJ and UPGMA trees of concatenated *PAD1*+*FDC1*.** Tip colors show the ecological group (step 04). The haplotype tree in Fig. 4 collapses identical sequences from these trees.
:::

::: figure
![](figures/pca_mds_kmeans_PAD1_FDC1.png){width=92%}

**Extended Data Fig. 5 | Strains in sequence space.** PCA (Lecture 16) and MDS of the one-hot-encoded variable alignment columns, colored by group, POF status and *k*-means cluster (Lecture 14; *k* chosen by silhouette score) (step 04). Clusters follow haplotype lineages rather than niches.
:::

::: figure
![](figures/05c_fdc1_contact_map.png){width=55%}

**Extended Data Fig. 6 | Fdc1 Cα–Cα contact map** (< 8 Å, PDB 4ZAC chain A) with truncation points in red: each common truncation removes most long-range contacts of the fold (step 05).
:::

::: figure
![](figures/05d_missense_distances.png){width=92%}

**Extended Data Fig. 7 | Missense variants are far from the catalytic site.** Distances from each missense site to prFMN, Mn^2+^ and the active site in Fdc1 (4ZAC), and positions in the Pad1 AlphaFold model with pLDDT (step 05).
:::

::: figure
![](figures/extended/ED_mal_cnv_1011.png){width=88%}

**Extended Data Fig. 8 | Maltose-transporter copy number across all ecological groups and clades of the 1,011-genome panel**, and frequency of *FDC1* deletion (copy number 0) by group (step 09).
:::

::: figure
![](figures/mal_copy_number.png){width=88%}

**Extended Data Fig. 9 | The superseded assembly-based *MAL/MPH* locus count** in {{R['core_mal_assembly'][3]['domesticated_n'] + R['core_mal_assembly'][3]['wild_lab_n']}} chromosome-level assemblies (step 04); assemblies collapse near-identical subtelomeric copies, so this count was replaced by read depth (Fig. 6).
:::

## Custom code

The complete custom code follows. Comments at the top of each step name the lecture its technique comes from.
:::
