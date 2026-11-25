# %%
# 2. Maltose / maltotriose transporter copy number in complete genomes
import os, io, re, time
import pandas as pd

EMAIL, TOOL = os.environ.get('NCBI_EMAIL', 'anonymous@example.org'), 'bch339n_project'
# file names into the new cnv folder
CNV = 'data/cnv'
STRAINS_CSV = f'{CNV}/complete_genome_strains.csv'
HITS_TSV = f'{CNV}/mal_blast_hits.tsv'
COPY_CSV = f'{CNV}/mal_copy_number.csv'
QFASTA = f'{CNV}/mal_query_cds.fasta'
GENES = {  # S288C CDS: chromosome record, start, stop, strand (1 = +, 2 = -)
    'MAL11': ('NC_001139.9', 1073963, 1075813, 2),
    'MAL31': ('NC_001134.8', 802631, 804475, 2),
    'MPH2':  ('NC_001136.10', 5985, 7814, 1),
    'MPH3':  ('NC_001142.9', 738008, 739816, 2)}
os.makedirs(CNV, exist_ok=True)

# %%
# Query CDSs (Entrez efetch of the S288C coordinates) - cached
if not os.path.exists(QFASTA):
    from Bio import Entrez
    Entrez.email, Entrez.tool = EMAIL, TOOL
    with open(QFASTA, 'w') as out:
        for g, (acc, s, e, strand) in GENES.items():
            fa = Entrez.efetch(db='nuccore', id=acc, rettype='fasta', retmode='text',
                               seq_start=s, seq_stop=e, strand=strand).read()
            out.write(f'>{g}_S288C {acc}:{s}-{e}\n' + ''.join(fa.split('\n')[1:]).strip() + '\n')
            time.sleep(0.4)
from Bio import SeqIO
QUERIES = {r.id.split('_')[0]: str(r.seq) for r in SeqIO.parse(QFASTA, 'fasta')}
QLEN = {g: len(s) for g, s in QUERIES.items()}
print(QLEN)

# %%
# Strain set: chromosome-level assemblies -> nuccore records (Entrez esearch / esummary / elink)
if not os.path.exists(STRAINS_CSV):
    from Bio import Entrez
    Entrez.email, Entrez.tool = EMAIL, TOOL
    term = ('(txid4932[orgn] OR txid27292[orgn]) AND ("complete genome"[filter] OR '
            '"chromosome level"[filter]) AND latest[filter]')
    ids = Entrez.read(Entrez.esearch(db='assembly', term=term, retmax=1000))['IdList']
    # shortcut: we hand-picked strains spanning niches from these summaries (contig N50 >= ~400 kb,
    # i.e. long-read assemblies); re-curate the list below if NCBI adds better ale assemblies.
    raise SystemExit(f'{len(ids)} assemblies found - curate {STRAINS_CSV} from their esummary/biosample records')
strains = pd.read_csv(STRAINS_CSV)
acc2strain = {a: row.strain for row in strains.itertuples() for a in row.nuccore_accessions.split(';')}
print(len(strains), 'assemblies;', strains.group.value_counts().to_dict())

# %%
# BLAST each assembly (records packed into <= 9.5 Mb batches: the bl2seq subject limit is 10 Mb)
def blast2seq(query_fasta, subject_accs):
    import requests
    url = 'https://blast.ncbi.nlm.nih.gov/'
    put = requests.post(url + 'BlastAlign.cgi', data=dict(
        CMD='Put', PROGRAM='blastn', MEGABLAST='on', BLAST_SPEC='blast2seq', QUERY=query_fasta,
        SUBJECTS='\n'.join(subject_accs), EXPECT='1e-50', HITLIST_SIZE=500, EMAIL=EMAIL, TOOL=TOOL)).text
    rid = re.search(r'name="RID" value="(\w+)"', put).group(1)
    while True:
        time.sleep(15)
        info = requests.get(url + 'Blast.cgi', params=dict(CMD='Get', FORMAT_OBJECT='SearchInfo', RID=rid)).text
        if 'Status=READY' in info: break
        if 'Status=WAITING' not in info: raise RuntimeError(rid + ' failed')
    txt = requests.get(url + 'Blast.cgi', params=dict(CMD='Get', FORMAT_TYPE='Text', ALIGNMENT_VIEW='Tabular',
                                                      DESCRIPTIONS=5000, ALIGNMENTS=5000, RID=rid)).text
    return [l.split('\t') for l in txt.splitlines() if re.match(r'^\S+\t', l) and not l.startswith('#')]

if not os.path.exists(HITS_TSV):
    from Bio import Entrez
    Entrez.email, Entrez.tool = EMAIL, TOOL
    qfa = ''.join(f'>{g}\n{s}\n' for g, s in QUERIES.items())
    rows = []
    for row in strains.itertuples():
        accs = row.nuccore_accessions.split(';')
        summ = Entrez.read(Entrez.esummary(db='nuccore', id=','.join(accs)))
        size = {d['AccessionVersion']: int(d['Length']) for d in summ}
        title = {d['AccessionVersion']: d['Title'] for d in summ}
        batches = []
        for a in sorted(accs, key=size.get, reverse=True):          # first-fit packing
            b = next((b for b in batches if b[0] + size[a] <= 9.5e6), None)
            if b is None: b = [0, []]; batches.append(b)
            b[0] += size[a]; b[1].append(a)
        for _, batch in batches:
            for h in blast2seq(qfa, batch):
                rows.append(h[:12] + [title.get(h[1], '')])
            time.sleep(10)
    cols = 'gene sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore title'.split()
    hits = pd.DataFrame(rows, columns=cols)
    for c in cols[2:12]: hits[c] = pd.to_numeric(hits[c])
    hits['qcov'] = (hits.qend - hits.qstart + 1) / hits.gene.map(QLEN)
    hits['strain'] = hits.sacc.map(acc2strain)
    hits[['gene', 'sacc', 'pident', 'length', 'qcov', 'sstart', 'send', 'qstart', 'qend', 'evalue',
          'bitscore', 'strain', 'title']].to_csv(HITS_TSV, sep='\t', index=False)
hits = pd.read_csv(HITS_TSV, sep='\t')
print(len(hits), 'HSPs (all, unfiltered)')

# %%
# Filter, then count copies: raw HSP counts and best-query assignment of merged loci
# filter and count best hits for genes of interest, 80% sequence identity and 70% of length
keep = hits[(hits.pident >= 80) & (hits.length >= 0.7 * hits.gene.map(QLEN))].copy()
# smallest chromosome coordinate is lo, larger coordinate is hi
keep['lo'] = keep[['sstart', 'send']].min(axis=1)
keep['hi'] = keep[['sstart', 'send']].max(axis=1)

# Score = bitscore / query length: raw bitscore favours the longer MPH2 query even at the MPH3 locus
# (MPH2 and MPH3 are 99.9% identical), so we normalise by the query's maximum possible length.
# normalize bitscores by dividing by gene total length
keep['score'] = keep.bitscore / keep.gene.map(QLEN)

def best_loci(df):
    """Greedy: best-scoring HSP first; an HSP overlapping an existing locus by >50% of the shorter joins it."""
    loci = []
    # sort hits by score
    for h in df.sort_values('score', ascending=False).itertuples():
        for L in loci:
            # if h overlaps any genes already saved in that chromosome, then don't save it
            if L['sacc'] == h.sacc and min(L['hi'], h.hi) - max(L['lo'], h.lo) > 0.5 * min(L['hi'] - L['lo'], h.hi - h.lo):
                break
        else:
            # otherwise, save it
            loci.append(dict(sacc=h.sacc, lo=h.lo, hi=h.hi, gene=h.gene, strain=h.strain))
    return pd.DataFrame(loci)

loci = best_loci(keep)
# MPH2/MPH3 differ at a handful of sites, so sequence cannot tell them apart outside S288C;
# name them by synteny instead: chr IV-left copy = MPH2, chr X-right copy = MPH3 (S288C positions).
title = keep.drop_duplicates('sacc').set_index('sacc').title
chrom = loci.sacc.map(title).str.extract(r'chromosome (?:S[ce])?(IV|X|4|10)\b', expand=False)
mph = loci.gene.isin(['MPH2', 'MPH3'])
loci.loc[mph & chrom.isin(['IV', '4']), 'gene'] = 'MPH2'
loci.loc[mph & chrom.isin(['X', '10']), 'gene'] = 'MPH3'
# count BLAST hits per gene of interest in total, then by loci to exclude duplicates
raw = keep.groupby(['strain', 'gene']).size().unstack(fill_value=0)
best = loci.groupby(['strain', 'gene']).size().unstack(fill_value=0)
genes = list(GENES)
# strains df from earlier being edited
cn = (strains[['strain', 'assembly_accession', 'group', 'species', 'source_niche']]
      .rename(columns={'assembly_accession': 'assembly_or_record_set'}).set_index('strain'))
# filling new cn df with raw hits and best hits by strain
for g in genes:
    cn[g] = best.get(g, pd.Series(dtype=int)).reindex(cn.index).fillna(0).astype(int)   # best = main count
    cn[g + '_raw'] = raw.get(g, pd.Series(dtype=int)).reindex(cn.index).fillna(0).astype(int)
# MPH2 vs MPH3 difference is near-arbitrary, especially when genes duplicate and jump
# better to have a total count of copies of the gene
cn['MPH2_3'] = cn.MPH2 + cn.MPH3
# sum all genes to measure overall MAL transporter capacity between strains
cn['total_loci'] = cn[genes].sum(axis=1)
cn.reset_index().to_csv(COPY_CSV, index=False)

# self-checks: S288C reference has exactly one locus of each gene; best <= raw everywhere
assert cn.loc['S288C', genes].tolist() == [1, 1, 1, 1], cn.loc['S288C', genes]
assert all((cn[g] <= cn[g + '_raw']).all() for g in genes)
# print for each gene of interest, raw and best, sort by group
cn[genes + [g + '_raw' for g in genes] + ['group']].sort_values('group')
