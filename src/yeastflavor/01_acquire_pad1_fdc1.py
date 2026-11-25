# %%
# 1. Data acquisition: S288C reference genes and PAD1/FDC1 alleles
import os, io, time
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Entrez.email = os.environ.get("NCBI_EMAIL", "anonymous@example.org")
Entrez.tool = "bch339n_project"
os.makedirs("data/reference", exist_ok=True); os.makedirs("data/blast", exist_ok=True); os.makedirs("data/seqs", exist_ok=True)

def efetch_region(acc, start, stop, strand=1):
    """FASTA text for acc:start-stop; strand=2 returns the reverse complement."""
    h = Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text",
                      seq_start=start, seq_stop=stop, strand=strand)
    txt = h.read(); h.close(); time.sleep(0.34)  # NCBI limit: 3 requests/s
    return txt

# %%
# 1a. S288C reference CDS
# making a reference sequence dictionary to make sure that later NCBI querys are matches to the protein strains we want
# the genes we want to analyze
genes = pd.DataFrame([
    ("PAD1", "YDR538W", "NC_001136.10", 1510902, 1511630, "+"),
    ("FDC1", "YDR539W", "NC_001136.10", 1512094, 1513605, "+"),
    ("MAL11", "YGR289C", "NC_001139.9", 1073963, 1075813, "-"),
    ("MAL31", "YBR298C", "NC_001134.8", 802631, 804475, "-"),
    ("MPH2", "YDL247W", "NC_001136.10", 5985, 7814, "+"),
    ("MPH3", "YJR160C", "NC_001142.9", 738008, 739816, "-"),
], columns=["gene", "systematic", "chr_acc", "start", "stop", "strand"])
# then storing genes to a data file to save them, we don't have to query every time we run this notebook
genes.to_csv("data/reference/genes.csv", index=False)

# dictionary of sequences
ref = {}
for g in genes.itertuples():
    path = f"data/reference/{g.gene}_cds.fasta"
    # check if this file is not in data to save time
    if not os.path.exists(path):
        rec = SeqIO.read(io.StringIO(efetch_region(g.chr_acc, g.start, g.stop, 1 if g.strand == "+" else 2)), "fasta")
        open(path, "w").write(f">{g.gene}_S288C {rec.description}\n{rec.seq}\n")
    # add gene sequence to dictionary as {g.gene: seq}
    ref[g.gene] = SeqIO.read(path, "fasta").seq
    # for printing below, translate sequence to protein
    prot = ref[g.gene].translate()
    # check that the first three aa of the sequence should be ATG, the sequence should be a multiple of 3, the protein should end as '*' and not have any '*' within it
    assert ref[g.gene][:3] == "ATG" and len(ref[g.gene]) % 3 == 0 and prot[-1] == "*" and "*" not in prot[:-1], g.gene
    # print results nicely
    print(f"{g.gene:6s} {len(ref[g.gene]):5d} bp  {len(prot) - 1:4d} aa  OK")

# %%
# 1b. BLAST PAD1 and FDC1 against core_nt (yeast only)
import requests
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
COLS = ["qacc", "sacc", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

def blast_tabular(name, seq, db="core_nt", entrez_query="txid4932[ORGN] OR txid27292[ORGN]"):
    put = requests.post(BLAST, data=dict(CMD="Put", PROGRAM="blastn", MEGABLAST="on", DATABASE=db,
                        QUERY=f">{name}\n{seq}", ENTREZ_QUERY=entrez_query, HITLIST_SIZE=500, EXPECT=1e-50,
                        EMAIL=Entrez.email, TOOL=Entrez.tool)).text
    rid = put.split("RID = ")[1].split()[0]
    while "Status=READY" not in requests.get(BLAST, params=dict(CMD="Get", FORMAT_OBJECT="SearchInfo", RID=rid)).text:
        time.sleep(10)
    txt = requests.get(BLAST, params=dict(CMD="Get", FORMAT_TYPE="Text", ALIGNMENT_VIEW="Tabular",
                       DESCRIPTIONS=500, ALIGNMENTS=500, RID=rid)).text
    rows = [l.split("\t")[:12] for l in txt.splitlines() if l.count("\t") >= 11 and not l.startswith("#")]
    return pd.DataFrame(rows, columns=COLS)

hits = {}
for gene in ["PAD1", "FDC1"]:
    path = f"data/blast/{gene}_core_nt.tsv"
    if not os.path.exists(path):
        blast_tabular(gene, ref[gene]).to_csv(path, sep="\t", index=False)
    hits[gene] = pd.read_csv(path, sep="\t")
    print(gene, len(hits[gene]), "hits;", "pident range", hits[gene].pident.min(), "-", hits[gene].pident.max())

# %%
# 1c. Beer, wine, sake, wild... genomes with known ecology (Gallone et al. 2016, Cell 166:1397)
from Bio import Align
aligner = Align.PairwiseAligner(mode="global", match_score=1, mismatch_score=-1, open_gap_score=-5, extend_gap_score=-1)

def best_match(S, k, start=0, max_mm=3):
    """Leftmost position in S (>= start) where k matches with the fewest mismatches (<= max_mm)."""
    best = (max_mm + 1, -1)
    for i in range(start, len(S) - len(k) + 1):
        mm = sum(a != b for a, b in zip(S[i:i + len(k)], k))
        if mm < best[0]:
            best = (mm, i)
            if mm == 0: break
    return best[1]

def locate_cds(S, cds):
    """Return (start, end) 0-based half-open of the CDS in S using start/stop 24-mer anchors, or None."""
    a = best_match(S, str(cds[:24]))
    if a < 0: return None
    b = best_match(S[a:a + len(cds) + 400], str(cds[-24:]), start=int(0.85 * len(cds)))
    return (a, a + b + 24) if b >= 0 else None

# self-check on a synthetic chromosome: flank + PAD1 with a 1-bp insertion + flank
_toy = "ACGT" * 50 + str(ref["PAD1"][:300]) + "A" + str(ref["PAD1"][300:]) + "TTGA" * 50
assert locate_cds(_toy, ref["PAD1"]) == (200, 200 + len(ref["PAD1"]) + 1)

if not all(os.path.exists(f"data/blast/{g}_gallone_chrIV.tsv") for g in ["PAD1", "FDC1"]):
    ids = Entrez.read(Entrez.esearch(db="nuccore", term='PRJNA323691[BioProject] AND "chromosome IV"[title]', retmax=500))["IdList"]
    chr4 = [(d["AccessionVersion"], int(d["Length"])) for d in Entrez.read(Entrez.esummary(db="nuccore", id=",".join(ids)))]
    rows = {"PAD1": [], "FDC1": []}
    for acc, L in chr4:
        st = max(1, L - 250000)
        S = str(SeqIO.read(io.StringIO(efetch_region(acc, st, L)), "fasta").seq).upper()
        for gene in rows:
            loc = locate_cds(S, ref[gene])
            if loc is None:
                print("not found:", acc, gene); continue
            c = aligner.align(ref[gene], S[loc[0]:loc[1]])[0].counts()
            n = c.identities + c.mismatches + c.gaps
            rows[gene].append([gene, acc, round(100 * c.identities / n, 3), n, c.mismatches, c.gaps,
                               1, len(ref[gene]), st + loc[0], st + loc[1] - 1, "NA", "NA"])
    for gene, rr in rows.items():
        pd.DataFrame(rr, columns=COLS).to_csv(f"data/blast/{gene}_gallone_chrIV.tsv", sep="\t", index=False)
for gene in hits:
    g2 = pd.read_csv(f"data/blast/{gene}_gallone_chrIV.tsv", sep="\t")
    print(gene, "Gallone chrIV alleles:", len(g2))
    hits[gene] = pd.concat([hits[gene], g2], ignore_index=True)

# %%
# 1d. Hit titles (Entrez esummary) and allele sequences (Entrez efetch)
# time to get the actual nt sequences
# check if we have it first
path = "data/blast/hit_titles.csv"
if not os.path.exists(path):
    # esummary first to see what hits we find from our hits dictionary
    accs = sorted(set(pd.concat(hits.values()).sacc))
    summ = Entrez.read(Entrez.esummary(db="nuccore", id=",".join(accs)))
    pd.DataFrame([(d["AccessionVersion"], d["Title"]) for d in summ], columns=["sacc", "title"]).to_csv(path, index=False)
titles = pd.read_csv(path).set_index("sacc").title

# now we can do the full efetch for each hit
for gene, df in hits.items():
    path = f"data/seqs/{gene}_alleles.fasta"
    if not os.path.exists(path):
        recs, seen = [SeqRecord(ref[gene], id="S288C_REF", description="S288C reference CDS")], {}
        for h in df.itertuples():
            lo, hi, minus = min(h.sstart, h.send), max(h.sstart, h.send), h.send < h.sstart
            seq = SeqIO.read(io.StringIO(efetch_region(h.sacc, lo, hi, 2 if minus else 1)), "fasta").seq
            seen[h.sacc] = seen.get(h.sacc, 0) + 1
            rid = h.sacc if seen[h.sacc] == 1 else f"{h.sacc}_{seen[h.sacc]}"
            # if less than 95% of sequence is present in hit, flag as partial, but keep it
            flag = (" PARTIAL" if (h.qend - h.qstart + 1) < 0.95 * len(ref[gene]) else "") + (" CONTAINS_N" if "N" in seq else "")
            recs.append(SeqRecord(seq, id=rid, description=f"{h.sstart}-{h.send} strand={'-' if minus else '+'}{flag} | {titles.get(h.sacc, '')}"))
        SeqIO.write(recs, path, "fasta")

# %%
# 1e. Sanity checks and PAD1/FDC1 pairing by strain
import re
# dictionary w FASTAs of pad1 and fdc1
alleles = {g: list(SeqIO.parse(f"data/seqs/{g}_alleles.fasta", "fasta")) for g in ["PAD1", "FDC1"]}
for g, recs in alleles.items():
    # match record id to length
    L = pd.Series({r.id: len(r.seq) for r in recs})
    # count how many records are x amount of nucleotides long
    print(f"{g}: {len(recs)} records; lengths {L.value_counts().to_dict()}")
    # not divisible by 3 (frameshift mutation candidate) and "partials" (premature stop candidate)
    print("   not %3:", list(L[L % 3 != 0].index), "| partial:", [r.id for r in recs if "PARTIAL" in r.description])

def strain_of(title):
    # find strain name by searching for "strain:" and finding the words up until chromosome or complete
    m = re.search(r"strain:? ([^,]+?)(?: chromosome.*|, complete.*)?$", title)
    if m: return m.group(1).strip()
    # if that fails, search for "cerevisiae", words that follow that until chromosome or genes
    m = re.search(r"cerevisiae (\S+) (?:chromosome|genes)", title)
    return m.group(1) if m else None

def strain_table(gene):
    # one row per strain; for RIB6002 (several alleles) keep all accessions joined by ';'
    # make table with one row per strain
    t = pd.DataFrame([(r.id, strain_of(r.description.split(" | ", 1)[-1])) for r in alleles[gene][1:]], columns=["sacc", "strain"])
    return t.dropna().groupby("strain").sacc.apply(";".join)

# make file with strain table including pad1 and fdc1, if it doesn't already exist
path = "data/seqs/paired_strains.csv"
if not os.path.exists(path):
    pd.concat([strain_table("PAD1").rename("pad1_sacc"), strain_table("FDC1").rename("fdc1_sacc")], axis=1, join="inner") \
      .rename_axis("strain").reset_index().to_csv(path, index=False)
paired = pd.read_csv(path)
paired
