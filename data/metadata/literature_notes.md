# Literature notes (A3) — verified from PubMed abstracts / PMC full text

Legend: [A] = verified from abstract, [FT] = verified from full text, [U] = not verified (memory / inference).

## 1. PAD1/FDC1 biology
- Mukai N, Masaki K, Fujii T, Kawamukai M, Iefuji H (2010) J Biosci Bioeng 109:564-569. PMID 20471595, doi:10.1016/j.jbiosc.2009.11.011
  - [A] Deletion/transformant work: adjacent PAD1 (YDR538W) and FDC1 (YDR539W) are BOTH essential for decarboxylation of ferulic, p-coumaric, cinnamic acid -> volatile phenols (4-VG = off-flavour in beer/wine). Pad1/Fdc1 are homologues of E. coli UbiX/UbiD; ubiquinone is still made in deletion mutants.
  - NOTE: the 2010 paper is the gene-function paper. The strain-panel SNP survey is the follow-up:
- Mukai N, Masaki K, Fujii T, Iefuji H (2014) J Biosci Bioeng 118:50-55. PMID 24507903, doi:10.1016/j.jbiosc.2013.12.017
  - [A] Panel: sake, shochu, wine, weizen, top-fermenting, bottom-fermenting and laboratory strains.
  - [A] Phenotype: most wine and weizen strains make 4-VG (POF+); sake, shochu, top-fermenting and bottom-fermenting strains do not (POF-).
  - [A] PAD1: nonsense mutation in standard top-fermenting yeast.
  - [A] FDC1: nonsense mutation in sake, shochu and standard top-fermenting yeasts; FRAMESHIFT in bottom-fermenting (lager) yeast.
  - [A] No nonsense/frameshift in laboratory, wine or weizen strains.
  - [A] Introducing FDC1 into sake/shochu strains restored ferulic-acid decarboxylation.
  - [U] The GenBank AB3687xx records (NCYC240, Weihenstephan 68, NCYC984, EC1118, K1-V1116, Kyokai wine 1/3/4, Kyokai 6, shochu...) are almost certainly this panel; cite Mukai 2014 for the alleles (could not open GenBank from here). Exact mutation positions not in the abstract -> derive them from our own alignments.
  - Expectation for our data: a 1513-bp FDC1 hit (vs 1512-bp S288C CDS) = +1 insertion frameshift, consistent with the lager (NCYC984) allele.

## 2. Canonical LOF alleles (what to look for in alignments)
- Mertens S, Gallone B, Steensels J, ... Verstrepen KJ (2019) PLoS One 14:e0209124. PMID 30625138, doi:10.1371/journal.pone.0209124 (correction PMID 31648239)
  - [FT] 73.80% of POF- strains of the Beer 1 lineage share FDC1 c.460C>T, CAA(Gln)->TAA stop = p.Q154*. Truncates before dimerisation domain/catalytic pocket.
  - [FT] Other natural LOF examples: beer strain BE011 homozygous FDC1 p.W497*; sake strain SP003 heterozygous p.Q154* and homozygous T>TA insertion at p.166 (frameshift).
  - [FT] POF is dominant: every PAD1 or FDC1 copy must be inactivated (ale ploidy averages 3.52).
  - [FT] W34/70 (S. pastorianus, lager, Germany) used as the POF- control; CBS1483 and W34/70 cited as industrial S. pastorianus.
  - [FT] Fdc1 decarboxylates ferulic acid; Pad1 makes the prenylated FMN cofactor for Fdc1.

## 3. Beer domestication
- Gallone B, Steensels J, ... Verstrepen KJ (2016) Cell 166:1397-1410. PMID 27610566, PMCID PMC5018251, doi:10.1016/j.cell.2016.08.020
  - [A] 157 industrial genomes; five industrial sublineages (Wine, Beer 1, Beer 2, Asia, Mixed) separated from wild strains.
  - [FT] Convergent selection: mutations and duplications of MAL genes; nonsense mutations in PAD1 and FDC1.
  - [FT] Many industrial strains have lost 4-VG production, while wild, bakery and bioethanol strains generally keep it. LOF mutations (SNPs and/or frameshift indels) in PAD1 and/or FDC1 are found especially in beer and sake strains, "never observed in strains from natural environments or bioethanol production". Different sublineages carry different disruptive mutations (convergence).
  - [FT] Ancestral-state reconstruction: an early Beer 1 subclade got LOF at the base of the clade; the Asian sake (not bioethanol) strains have their own LOF.
  - [FT] Hefeweizen strains (BE072, BE074, BE093) sit in Beer 1 but are mosaic; ~8-13% of the genome comes from Wine, including the subtelomeric ChrIV region with functional PAD1/FDC1 -> POF+ wheat beer yeast.
  - [FT] Maltotriose: beer yeasts use maltotriose much better. This tracks a specific allele of a high-affinity maltotriose transporter (AGT1/MAL11 [U: gene names were removed from the PMC text; by context it is AGT1]) present only in Beer 1 and some mosaics; the whole MAL locus is absent in the Wine subpopulation. Beer 2 strains ferment maltotriose but carry frameshifts in that transporter and lower MAL-locus CNV, so the trait evolved convergently.
  - [FT] Beer 1 last common ancestor dated to AD 1573-1604; Beer 2 to AD 1645-1671.
  - Fraction POF- in beer: no single number in the text I read. Use 73.8% of POF- Beer 1 strains carrying Q154* (Mertens 2019) and the qualitative "many/most beer strains POF-, except Hefeweizen".
- Gonçalves M, ... Sampaio JP (2016) Curr Biol 26:2750-2761. PMID 27720622, doi:10.1016/j.cub.2016.08.040
  - [A] Top-fermenting beer yeasts are polyphyletic. The main clade has at least three subgroups (German, British, wheat beer); other beer strains sit near sake, wine or bread yeasts. Beer domestication is distinct from wine and sake. The main beer clade's nucleotide diversity is more than double wine's. POF is a keyword (signature gene).
- Peter J, ... Liti G, Schacherer J (2018) Nature 556:339-344. PMID 29643504, doi:10.1038/s41586-018-0030-5
  - [A] 1,011 genomes; single "out-of-China" origin followed by several independent domestications. Domesticated isolates vary widely in ploidy, aneuploidy and gene content. Copy-number changes are the largest class of GWAS variants and have bigger effects than SNPs (supports a CNV analysis of MAL/MPH).

## 4. Maltose/maltotriose transporters in ale vs lager
- Vidgren V, Ruohonen L, Londesborough J (2005) Appl Environ Microbiol 71:7846-7857. PMID 16332759, doi:10.1128/AEM.71.12.7846-7857.2005
  - [A] Two ale and three lager strains: all had complete MAL1-4 loci (one ale lacked MAL2), plus AGT1 and MAL11 alleles.
  - [A] **MPH genes were present in the LAGER but not the ALE strains** -> directly relevant to our MPH2/MPH3 counts.
  - [A] Ale AGT1 = full length (616 aa); lager AGT1 = truncated (394 aa, premature stop) and likely non-functional. Lager strains rely on MALx1 transporters, ale strains on Agt1.
- Vidgren V, ... Londesborough J (2009) Appl Environ Microbiol 75:2333-2345. PMID 19181838, doi:10.1128/AEM.01558-08
  - [A] Repairing the lager AGT1 premature stop with ale sequence improved VHG wort fermentation (more maltose and maltotriose use).
- Magalhães F, Vidgren V, Ruohonen L, Gibson B (2016) FEMS Yeast Res 16:fow053. PMID 27364826, doi:10.1093/femsyr/fow053
  - [A] Group I (Saaz) lager strains vary. Maltotriose users carry SeAGT1 and MTT1. Strain A58 has higher MTT1 copy number. CBS1513 = group I type strain.
- Brouwers N, ... Daran JMG (2019) Appl Environ Microbiol 85:e01516-19. PMID 31519660, doi:10.1128/AEM.01516-19
  - [A] Himalayan (Holarctic) S. eubayanus, the proposed S. pastorianus parent, carries several copies of an alpha-oligoglucoside transporter but cannot grow on maltose or maltotriose without a functional MAL regulator. Lager maltotriose use is heterotic (cross-talk between subgenomes).

## 5. Other relevant facts
- Badura J, ... Wendland J (2025) Appl Microbiol Biotechnol 109:17. PMID 39841271, doi:10.1007/s00253-024-13397-8: [A] Kveik strain "Odin" shares inactivating PAD1 and FDC1 mutations with lager yeasts (POF-); POF- is "a key feature of lager yeasts".
- Nasuti C, ... Krogerus K, Solieri L (2026) Int J Food Microbiol 456:111786. PMID 42025082, doi:10.1016/j.ijfoodmicro.2026.111786: [A] a POF-reduced lager-like hybrid lost the S. eubayanus PAD1-FDC1 cluster and had an S. cerevisiae FDC1 frameshift.
- Genome sources for strains: Borneman 2011 PLoS Genet PMID 21304888 (4 wine + 2 brewing strains; [U] list believed AWRI796, QA23, VIN7, VIN13, Foster's B, Foster's O; AWRI1631 is from an earlier 2008 paper); Akao 2011 DNA Res PMID 21900213 (sake K7); Zhang 2018 DNA Res PMID 29415277 (N85, Chinese rice wine); Ulaganathan 2015 Genome Announc PMID 25676764 (NCIM3107, bioethanol); Shao 2018 Nature PMID 30069045 (SY14 single-chromosome); Liti 2009 Nature PMID 19212322 (SGRP: YPS128, DBVPG6044, DBVPG6765, UWOPS03-461.4, Y12, SK1...).

## Expected results (for the report)
1. Lab, wine, wild, bioethanol and Hefeweizen (W68) strains: intact PAD1/FDC1 ORFs (POF+).
2. Sake (Kyokai 6/7), shochu and awamori strains: FDC1 premature stop (POF-).
3. Standard ale (NCYC240): PAD1 nonsense + FDC1 nonsense. Beer 1 ales are expected to carry FDC1 c.460C>T / p.Q154*.
4. Lager (NCYC984, W34/70): FDC1 frameshift on the S. cerevisiae subgenome (length 1513 instead of 1512 fits a +1 indel).
5. MAL/MPH: lager strains carry MPH genes that ale strains lack (Vidgren 2005); wine strains can lack whole MAL loci (Gallone 2016); lager AGT1 is truncated (394 aa).

## Surprises / caveats
- The "wine yeast Kyokai No.1/3/4" records are Brewing Society of Japan WINE yeasts, not the historical sake Kyokai No.1-5. Classify them as Wine (Mukai found wine strains intact). Do not merge them with sake.
- Mukai 2010 does not contain the strain survey; the SNP/strain data are in Mukai 2014 (PMID 24507903).
- Hefeweizen is the expected exception: a beer yeast that is POF+.
