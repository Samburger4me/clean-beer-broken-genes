# Technique source: Notebook 08, BLAST and Biopython (PairwiseAligner, Seq.translate); Notebook 09, Entrez — notebook_notes/lecture_08.md; notebook_notes/lecture_09.md
# %%
# 10. The S. eubayanus half of lager yeast: are its PAD1 / FDC1 copies intact?
# Libraries used in this section
import os, re, ast, time, hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data.CodonTable import standard_dna_table

# NCBI asks every script to say who it is; the e-mail comes from the environment.
EMAIL, TOOL = os.environ.get('NCBI_EMAIL', 'anonymous@example.org'), 'bch339n_project'
Entrez.email, Entrez.tool = EMAIL, TOOL
LAGER = 'data/lager'
SE_CDS_FA = f'{LAGER}/Seub_FM1318_PAD1_FDC1_cds.fasta'   # Se reference CDSs
HITS_TSV = f'{LAGER}/lager_pad1_fdc1_blast_hits.tsv'     # raw blast2seq HSPs
BATCH_TSV = f'{LAGER}/lager_blast_batches.tsv'          # which records went in which BLAST job
LOCI_FA = f'{LAGER}/lager_pad1_fdc1_loci.fasta'          # efetched locus sequences
OUT_CSV = 'data/results/lager_subgenome_pad1_fdc1.csv'
FIG_DIR = 'figures/extended'
for d in (LAGER, FIG_DIR, 'data/results'):
    os.makedirs(d, exist_ok=True)
SE_ASSEMBLY = 'GCF_001298625.1'   # S. eubayanus FM1318 SEUB3.0 (Baker et al. 2015), RefSeq reference
SE_TAXID = 1080349                # Saccharomyces eubayanus
GENES = ['PAD1', 'FDC1']

# Thresholds (our choice). A copy belongs to a subgenome if it matches that subgenome's reference at
# >= 98 % and the other reference at < 90 %. Sc and Se orthologs are only ~78-82 % identical
# (step 2 below), so the gap between 90 and 98 leaves room for strain-level SNPs - for example a
# lager Sc copy need not be 99 % identical to the lab strain S288C. MIN_QCOV is section 2's filter.
SC_MIN_ID = 98.0
SE_MIN_ID = 98.0
OTHER_MAX_ID = 90.0
MIN_QCOV = 0.7
BLAST_EVALUE = '1e-5'   # sensitive on purpose: a missing copy must not be a search artefact
PAD_NT = 30        # extra flank fetched on each side so indels near the CDS ends are not cut off
DATE = '2026-11-18'

# %%
def add_manifest(source, acc_or_url, path, notes):
    """Append one row to data/MANIFEST.csv for a cached file, unless that path is already logged."""
    man = pd.read_csv('data/MANIFEST.csv')
    if path in set(man.local_path):
        return
    sha = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    row = dict(source=source, accession_or_url=acc_or_url, local_path=path,
               date_downloaded=DATE, sha256=sha, notes=notes)
    pd.concat([man, pd.DataFrame([row])]).to_csv('data/MANIFEST.csv', index=False)

# Step 0: the S288C reference CDSs (same files section 3 uses).
sc_cds = {g: str(next(SeqIO.parse(f'data/reference/{g}_cds.fasta', 'fasta')).seq).upper()
          for g in GENES}
print('S288C CDS lengths (nt):', {g: len(s) for g, s in sc_cds.items()})

# %%
# 10.1 S. eubayanus reference CDSs
# Step 1: find the Se PAD1/FDC1 genes in NCBI Gene and fetch their CDS (cached).
if not os.path.exists(SE_CDS_FA):
    term = f'txid{SE_TAXID}[orgn] AND (PAD1[sym] OR FDC1[sym])'
    ids = Entrez.read(Entrez.esearch(db='gene', term=term, retmax=10))['IdList']
    summ = Entrez.read(Entrez.esummary(db='gene', id=','.join(ids)), validate=False)
    recs = []
    for d in summ['DocumentSummarySet']['DocumentSummary']:
        gi = d['GenomicInfo'][0]
        a, b = int(gi['ChrStart']), int(gi['ChrStop'])
        strand = 1 if b > a else 2             # NCBI Gene lists minus-strand genes as start > stop
        lo, hi = min(a, b) + 1, max(a, b) + 1  # 0-based -> 1-based
        fa = Entrez.efetch(db='nuccore', id=gi['ChrAccVer'], rettype='fasta', retmode='text',
                           seq_start=lo, seq_stop=hi, strand=strand).read()
        seq = ''.join(fa.split('\n')[1:]).strip()
        desc = (f"{gi['ChrAccVer']}:{lo}-{hi} strand={strand} GeneID:{d.attributes['uid']} "
                f"locus_tag:{d['OtherAliases']} chromosome {d['Chromosome']} {d['Description']}")
        recs.append(f">Se_{d['Name']} {desc}\n{seq}\n")
        time.sleep(0.4)
    open(SE_CDS_FA, 'w').write(''.join(recs))
add_manifest('NCBI Gene + nuccore (S. eubayanus FM1318, RefSeq ' + SE_ASSEMBLY + ', Baker et al. 2015)',
             'Entrez gene esearch txid1080349 PAD1/FDC1 -> efetch nuccore CDS ranges',
             SE_CDS_FA, 'S. eubayanus PAD1 (DI49_4415) and FDC1 (DI49_4416) CDS')
se_recs = {r.id.split('_')[1]: r for r in SeqIO.parse(SE_CDS_FA, 'fasta')}
se_cds = {g: str(se_recs[g].seq).upper() for g in GENES}
for g in GENES:
    print(f'Se {g}:', se_recs[g].description)

# Check that what I fetched really is a CDS: a multiple of 3, starts ATG, ends in a stop and has no
# internal stop. If the 0/1-based conversion were wrong, these would fail.
for g, s in se_cds.items():
    p = str(Seq(s).translate())
    assert len(s) % 3 == 0 and s.startswith('ATG') and p.endswith('*') and p.count('*') == 1, g

# %%
# Step 2: how far apart are the Sc and Se orthologs? I use the class-style global aligner
# (attribute style, BLOSUM62 / open -10 / extend -1 as in the student's HW04, lecture_08.md).
prot_class = PairwiseAligner()
prot_class.mode = 'global'
prot_class.substitution_matrix = substitution_matrices.load('BLOSUM62')
prot_class.open_gap_score = -10
prot_class.extend_gap_score = -1

def pct_id(aln):
    """Percent identity over aligned (non-gap) columns of a 2-row alignment."""
    cols = [(x, y) for x, y in zip(aln[0], aln[1]) if x != '-' and y != '-']
    return 100 * sum(x == y for x, y in cols) / len(cols)

# Nucleotide aligner: same scoring as section 3 (match 2, mismatch -3, open -10, extend -1).
nt_class = PairwiseAligner(mode='global', match_score=2, mismatch_score=-3,
                           open_gap_score=-10, extend_gap_score=-1)
sc_se = []
for g in GENES:
    nt = pct_id(nt_class.align(sc_cds[g], se_cds[g])[0])
    aa = pct_id(prot_class.align(str(Seq(sc_cds[g]).translate(to_stop=True)),
                                 str(Seq(se_cds[g]).translate(to_stop=True)))[0])
    sc_se.append(dict(gene=g, sc_len_nt=len(sc_cds[g]), se_len_nt=len(se_cds[g]),
                      nt_identity=round(nt, 2), aa_identity=round(aa, 2)))
    RESULTS[f'a5_sc_vs_se_{g}_nt_identity_pct'] = round(nt, 2)
    RESULTS[f'a5_sc_vs_se_{g}_aa_identity_pct'] = round(aa, 2)
# The two orthologs should be clearly different (~80-90 %), otherwise BLAST could not tell
# the subgenomes apart by identity.
assert all(r['nt_identity'] < SE_MIN_ID for r in sc_se), 'Sc and Se references are too similar'

# Step 3: where do the Se genes sit on their chromosome? A gene close to a chromosome end lives in a
# subtelomeric region, which is more often lost or rearranged, so I note the distance to the end.
se_chr = se_recs['FDC1'].description.split(' ')[1].split(':')[0]
se_chr_len = int(Entrez.read(Entrez.esummary(db='nuccore', id=se_chr))[0]['Length'])
se_end = max(int(se_recs[g].description.split(':')[1].split('-')[1].split(' ')[0]) for g in GENES)
print(f'Se reference chromosome {se_chr} length (nt): {se_chr_len}')
print(f'Distance from Se PAD1-FDC1 cluster to the chromosome end (nt): {se_chr_len - se_end}')
RESULTS['a5_se_ref_cluster_to_chr_end_nt'] = se_chr_len - se_end
pd.DataFrame(sc_se)

# %%
# 10.2 BLAST both references against each lager assembly
strains_cnv = pd.read_csv('data/cnv/complete_genome_strains.csv')
lagers = strains_cnv[strains_cnv.group == 'Lager'].reset_index(drop=True)
# The task names exactly these two assemblies; stop if the strain table ever changes under us.
assert sorted(lagers.assembly_accession) == ['GCA_011022315.1', 'GCA_056824455.1'], lagers
acc2strain = {a: r.strain for r in lagers.itertuples() for a in r.nuccore_accessions.split(';')}

def blast2seq(query_fasta, subject_accs):
    """Run one blastn align-two-sequences job on NCBI and return the tabular HSP rows."""
    import requests
    url = 'https://blast.ncbi.nlm.nih.gov/'
    put = requests.post(url + 'BlastAlign.cgi', data=dict(
        CMD='Put', PROGRAM='blastn', BLAST_SPEC='blast2seq', QUERY=query_fasta,
        SUBJECTS='\n'.join(subject_accs), EXPECT=BLAST_EVALUE, HITLIST_SIZE=500,
        EMAIL=EMAIL, TOOL=TOOL)).text
    rid = re.search(r'name="RID" value="(\w+)"', put).group(1)
    # NCBI asks clients to poll no faster than every ~10 s
    while True:
        time.sleep(15)
        info = requests.get(url + 'Blast.cgi', params=dict(
            CMD='Get', FORMAT_OBJECT='SearchInfo', RID=rid)).text
        if 'Status=READY' in info:
            break
        if 'Status=WAITING' not in info:
            raise RuntimeError(rid + ' failed')
    txt = requests.get(url + 'Blast.cgi', params=dict(
        CMD='Get', FORMAT_TYPE='Text', ALIGNMENT_VIEW='Tabular', DESCRIPTIONS=5000,
        ALIGNMENTS=5000, RID=rid)).text
    return [l.split('\t') for l in txt.splitlines()
            if re.match(r'^\S+\t', l) and not l.startswith('#')]

queries = {f'Sc_{g}': sc_cds[g] for g in GENES} | {f'Se_{g}': se_cds[g] for g in GENES}
# Positive control: section 2's MAL/MPH queries hit chromosomes of both subgenomes in these lagers.
# If a batch returns no control hit, that batch was not really searched and "absent" means nothing.
ctrl = {f'CTRL_{r.id.split("_")[0]}': str(r.seq)
        for r in SeqIO.parse('data/cnv/mal_query_cds.fasta', 'fasta')}
QLEN = {q: len(s) for q, s in (queries | ctrl).items()}

# Step 3: BLAST (cached). If the TSV exists I skip the network entirely.
if not os.path.exists(HITS_TSV):
    qfa = ''.join(f'>{q}\n{s}\n' for q, s in (queries | ctrl).items())
    rows, batch_rows = [], []
    for r in lagers.itertuples():
        accs = r.nuccore_accessions.split(';')
        summ = Entrez.read(Entrez.esummary(db='nuccore', id=','.join(accs)))
        size = {d['AccessionVersion']: int(d['Length']) for d in summ}
        title = {d['AccessionVersion']: d['Title'] for d in summ}
        # first-fit packing, largest records first (bl2seq subject limit is 10 Mb)
        batches = []
        for a in sorted(accs, key=size.get, reverse=True):
            b = next((b for b in batches if b[0] + size[a] <= 9.5e6), None)
            if b is None:
                b = [0, []]
                batches.append(b)
            b[0] += size[a]
            b[1].append(a)
        for bi, (_, batch) in enumerate(batches):
            bid = f"{r.strain.replace(' ', '')}_b{bi}"
            batch_rows.append(dict(strain=r.strain, batch=bid, accessions=';'.join(batch)))
            for h in blast2seq(qfa, batch):
                rows.append(h[:12] + [title.get(h[1], ''), bid])
            time.sleep(10)
    cols = 'query sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore title batch'
    pd.DataFrame(rows, columns=cols.split()).to_csv(HITS_TSV, sep='\t', index=False)
    pd.DataFrame(batch_rows).to_csv(BATCH_TSV, sep='\t', index=False)
add_manifest('NCBI BLAST URL API (blastn, BLAST_SPEC=blast2seq)',
             'https://blast.ncbi.nlm.nih.gov/BlastAlign.cgi; subjects = nuccore records of '
             'GCA_011022315.1 (CBS 1483) and GCA_056824455.1 (CBS 1513)', HITS_TSV,
             'S288C + S. eubayanus FM1318 PAD1/FDC1 CDS + MAL/MPH control CDS vs lager '
             'chromosomes; blastn E<=1e-5')
add_manifest('NCBI BLAST URL API (job log)', 'https://blast.ncbi.nlm.nih.gov/BlastAlign.cgi',
             BATCH_TSV, 'which nuccore records were packed into each blast2seq job (<=9.5 Mb)')
hits_all = pd.read_csv(HITS_TSV, sep='\t')
batches_df = pd.read_csv(BATCH_TSV, sep='\t')
# every record of both assemblies must have been in exactly one batch
assert sorted(';'.join(batches_df.accessions).split(';')) == sorted(acc2strain), 'records missing'
ctrl_per_batch = (hits_all[hits_all['query'].str.startswith('CTRL_')].groupby('batch').size()
                  .reindex(batches_df.batch, fill_value=0))
print('Control (MAL/MPH) HSPs per BLAST batch:', ctrl_per_batch.to_dict())
assert (ctrl_per_batch > 0).all(), 'a BLAST batch returned no control hit - search not trustworthy'
RESULTS['a5_blast_batches_with_control_hit'] = f'{(ctrl_per_batch > 0).sum()}/{len(ctrl_per_batch)}'
hits = hits_all[~hits_all['query'].str.startswith('CTRL_')].copy()
hits['strain'] = hits.sacc.map(acc2strain)
hits['gene'] = hits['query'].str.split('_').str[1]
hits['qcov'] = (hits.qend - hits.qstart + 1) / hits['query'].map(QLEN)
print('HSPs (all):', len(hits))
hits[['query', 'strain', 'sacc', 'pident', 'length', 'qcov', 'sstart', 'send', 'title']].round(3)

# %%
# 10.3 Assign each locus to a subgenome
# Step 1: keep near-full-length HSPs (same coverage cut as section 2).
keep = hits[hits.qcov >= MIN_QCOV].copy()
keep['lo'] = keep[['sstart', 'send']].min(axis=1)
keep['hi'] = keep[['sstart', 'send']].max(axis=1)

# Step 2: merge overlapping HSPs of the same gene on the same record into loci.
loci = []
for h in keep.sort_values('bitscore', ascending=False).itertuples():
    L = next((L for L in loci if L['sacc'] == h.sacc and L['gene'] == h.gene
              and min(L['hi'], h.hi) > max(L['lo'], h.lo)), None)
    if L is None:
        L = dict(strain=h.strain, gene=h.gene, sacc=h.sacc, title=h.title, lo=h.lo, hi=h.hi,
                 best=h, id_vs_Sc=np.nan, id_vs_Se=np.nan)
        loci.append(L)
    sp = h.query.split('_')[0]
    L[f'id_vs_{sp}'] = np.nanmax([L[f'id_vs_{sp}'], h.pident])

# Step 3: the subgenome call from identity alone, plus the title's chromosome label as a check.
for L in loci:
    # a missing identity (no HSP of that query at this locus) counts as 0.
    sc_id, se_id = np.nan_to_num(L['id_vs_Sc']), np.nan_to_num(L['id_vs_Se'])
    sc_ok = sc_id >= SC_MIN_ID and se_id < OTHER_MAX_ID
    se_ok = se_id >= SE_MIN_ID and sc_id < OTHER_MAX_ID
    L['subgenome'] = 'Sc' if sc_ok else 'Se' if se_ok else 'ambiguous'
    L['chromosome_label'] = re.search(r'chromosome (\S+)', L['title']).group(1)
loci_df = pd.DataFrame(loci)
print('Loci found:', len(loci_df))
# An ambiguous locus would mean the identity rule cannot separate the subgenomes - stop and look.
assert (loci_df.subgenome != 'ambiguous').all(), loci_df[['strain', 'gene', 'id_vs_Sc', 'id_vs_Se']]
loci_df[['strain', 'gene', 'subgenome', 'sacc', 'chromosome_label', 'lo', 'hi',
         'id_vs_Sc', 'id_vs_Se']]

# %%
# 10.4 Extract each copy and call intact vs loss-of-function
# Step 1: pull the allele-caller definitions out of section 3's source (functions + aligner setup).
WANT = {'aa_of', 'right_shift_ins', 'right_shift_del', 'call_allele'}
tree3 = ast.parse(open('src/yeastflavor/03_variant_calling.py', encoding='utf-8').read())
nodes = [n for n in tree3.body
         if (isinstance(n, ast.FunctionDef) and n.name in WANT)
         or (isinstance(n, ast.Assign) and re.match(r'(nt_aligner|prot_aligner|CODON)\b',
                                                    ast.unparse(n.targets[0])))]
exec(compile(ast.Module(body=nodes, type_ignores=[]), '03_variant_calling.py', 'exec'), globals())
assert all(k in globals() for k in WANT | {'nt_aligner', 'prot_aligner', 'CODON'})

# call_allele looks its reference up in ref_cds / ref_prot by name, so I key them per subgenome.
ref_cds = {f'{sp}_{g}': s for sp, d in (('Sc', sc_cds), ('Se', se_cds)) for g, s in d.items()}
ref_prot = {k: str(Seq(s).translate(to_stop=True)) for k, s in ref_cds.items()}

# Step 2: fetch each locus (cached as one FASTA, ids = strain|gene|subgenome).
def locus_range(L):
    """Full-CDS range (1-based, padded) and efetch strand for a locus from its best HSP."""
    h = L['best']
    ql = QLEN[h.query]
    if h.send > h.sstart:   # plus strand
        return h.sstart - (h.qstart - 1) - PAD_NT, h.send + (ql - h.qend) + PAD_NT, 1
    return h.send - (ql - h.qend) - PAD_NT, h.sstart + (h.qstart - 1) + PAD_NT, 2

loci_df['key'] = loci_df.strain + '|' + loci_df.gene + '|' + loci_df.subgenome
if not os.path.exists(LOCI_FA):
    out = []
    for L in loci_df.itertuples():
        lo, hi, strand = locus_range(L._asdict())
        fa = Entrez.efetch(db='nuccore', id=L.sacc, rettype='fasta', retmode='text',
                           seq_start=lo, seq_stop=hi, strand=strand).read()
        out.append(f">{L.key.replace(' ', '_')} {L.sacc}:{lo}-{hi} strand={strand}\n"
                   + ''.join(fa.split('\n')[1:]).strip() + '\n')
        time.sleep(0.4)
    open(LOCI_FA, 'w').write(''.join(out))
add_manifest('NCBI nuccore efetch (lager PAD1/FDC1 loci)',
             'Entrez efetch nuccore seq_start/seq_stop/strand on CP0489xx / CP1500xx records',
             LOCI_FA, 'Sc and Se PAD1/FDC1 loci +-30 nt from CBS 1483 and CBS 1513')
loc_seqs = {r.id.replace('_', ' ', 1) if r.id.startswith('CBS_') else r.id: r
            for r in SeqIO.parse(LOCI_FA, 'fasta')}
assert set(loc_seqs) == set(loci_df.key), 'cached loci FASTA does not match the BLAST loci'

# Step 3: run the section-3 caller on each copy. I strip the description so the PARTIAL/N flags of
# section 3 (which come from GenBank titles) do not apply - these are complete-genome records.
calls = []
for L in loci_df.itertuples():
    rec = loc_seqs[L.key]
    variants, sm, protein, aln = call_allele(f'{L.subgenome}_{L.gene}', SeqRecord(rec.seq, id=rec.id,
                                                                                  description=''))
    fs_or_stop = [v['hgvs_c'] + ' ' + v['protein_change'] for v in variants
                  if v.get('effect') in ('frameshift', 'nonsense')]
    calls.append(dict(strain=L.strain, gene=L.gene, subgenome=L.subgenome, sacc=L.sacc,
                      chromosome_label=L.chromosome_label, region=rec.description.split(' ')[1],
                      cds_start=L.lo, cds_end=L.hi,
                      blast_id_vs_Sc_ref=L.id_vs_Sc, blast_id_vs_Se_ref=L.id_vs_Se,
                      reference=f'{L.subgenome}_{L.gene}', ref_len_aa=len(ref_prot[f'{L.subgenome}_{L.gene}']),
                      nt_identity_to_ref=sm['pct_identity'], protein_len_aa=sm['protein_len'],
                      frameshift=sm['frameshift'], n_missense=sm['n_missense'],
                      protein_identity_to_ref=sm['protein_pct_identity'],
                      lof_variants='; '.join(fs_or_stop), call=sm['LOF_call']))
calls_df = pd.DataFrame(calls)

# Step 4: make the table complete - one row per strain x gene x subgenome, 'absent' if no locus.
grid = pd.MultiIndex.from_product([lagers.strain, GENES, ['Sc', 'Se']],
                                  names=['strain', 'gene', 'subgenome']).to_frame(index=False)
status = grid.merge(calls_df, how='left', on=['strain', 'gene', 'subgenome'])
status['call'] = status['call'].fillna('absent')
status.to_csv(OUT_CSV, index=False)
print('Saved:', OUT_CSV)
status[['strain', 'gene', 'subgenome', 'chromosome_label', 'cds_start', 'cds_end',
        'nt_identity_to_ref', 'protein_len_aa', 'ref_len_aa', 'lof_variants', 'call']]

# %%
# Step 5: for a strain with no Sc copy, I list its chromosome-IV records (Sc PAD1/FDC1 sit near the
# right end of chromosome IV in S288C), because a translocation in the title would explain a loss.
for st in status[(status.subgenome == 'Sc') & (status.call == 'absent')].strain.unique():
    accs = lagers.set_index('strain').loc[st, 'nuccore_accessions'].split(';')
    for d in Entrez.read(Entrez.esummary(db='nuccore', id=','.join(accs))):
        if re.search(r'chromosome \S*ScIV\b', d['Title']):
            print(f"{st} chromosome IV record: {d['AccessionVersion']} {d['Title']} "
                  f"({int(d['Length'])} nt)")

# %%
# Step 6: checks and the citable numbers.
# 1) Section 3 found FDC1 c.501dupA in lager alleles; the Sc FDC1 copy here should carry it too.
sc_fdc1 = status[(status.gene == 'FDC1') & (status.subgenome == 'Sc')]
for r in sc_fdc1[sc_fdc1.call != 'absent'].itertuples():
    has_dup = isinstance(r.lof_variants, str) and 'c.501dupA' in r.lof_variants
    print(f'{r.strain} Sc FDC1 carries c.501dupA:', has_dup)
    RESULTS[f"a5_{r.strain.replace(' ', '')}_Sc_FDC1_c501dupA"] = has_dup
# 2) one number per strain x gene x subgenome
for r in status.itertuples():
    k = f"a5_{r.strain.replace(' ', '')}_{r.subgenome}_{r.gene}"
    RESULTS[k + '_call'] = r.call
    if r.call != 'absent':
        RESULTS[k + '_locus'] = f'{r.sacc}:{int(r.cds_start)}-{int(r.cds_end)}'
        RESULTS[k + '_protein_len_aa'] = int(r.protein_len_aa)
        RESULTS[k + '_nt_identity_to_ref_pct'] = float(r.nt_identity_to_ref)
RESULTS['a5_se_reference_assembly'] = SE_ASSEMBLY
RESULTS['a5_se_reference_loci'] = {g: se_recs[g].description.split(' ')[1] for g in GENES}
RESULTS['a5_n_se_copies_intact'] = int(((status.subgenome == 'Se') & (status.call == 'intact')).sum())
RESULTS['a5_n_se_copies_lof'] = int(((status.subgenome == 'Se') & (status.call == 'LOF')).sum())
RESULTS['a5_n_se_copies_absent'] = int(((status.subgenome == 'Se') & (status.call == 'absent')).sum())
print('Se copies intact / LOF / absent:', RESULTS['a5_n_se_copies_intact'],
      RESULTS['a5_n_se_copies_lof'], RESULTS['a5_n_se_copies_absent'])

# %%
# Extended Data figure: subgenome status grid
# Nature Extended Data style: Arial 6-7 pt, 89 mm single column, bold lowercase panel label.
plt.rcParams.update({'font.family': 'Arial', 'font.size': 6, 'axes.linewidth': 0.5,
                     'pdf.fonttype': 42})
CALL_COL = {'intact': '#0072B2', 'LOF': '#D55E00', 'absent': '#BBBBBB'}   # Okabe-Ito blue / vermillion
cols = [(g, sp) for g in GENES for sp in ('Sc', 'Se')]
strains_order = list(lagers.strain)
mm = 1 / 25.4
fig, ax = plt.subplots(figsize=(89 * mm, 38 * mm))
for i, st in enumerate(strains_order):
    for j, (g, sp) in enumerate(cols):
        r = status[(status.strain == st) & (status.gene == g) & (status.subgenome == sp)].iloc[0]
        ax.add_patch(plt.Rectangle((j, i), 0.94, 0.9, color=CALL_COL[r.call], lw=0))
        txt = r.call if r.call == 'absent' else f'{r.call}\n{int(r.protein_len_aa)}/{int(r.ref_len_aa)} aa'
        if r.call == 'LOF' and isinstance(r.lof_variants, str) and r.lof_variants:
            txt += '\n' + r.lof_variants.split(';')[0].split(' ')[0]
        # white text is unreadable on the light grey 'absent' cells, so those get dark text
        ax.text(j + 0.47, i + 0.45, txt, ha='center', va='center', fontsize=5.5,
                color='#333333' if r.call == 'absent' else 'white')
ax.set_xlim(0, len(cols))
ax.set_ylim(len(strains_order), 0)
ax.set_xticks([j + 0.47 for j in range(len(cols))])
ax.set_xticklabels([f'{g}\n{sp} copy' for g, sp in cols], fontstyle='italic')
ax.set_yticks([i + 0.45 for i in range(len(strains_order))])
ax.set_yticklabels(strains_order)
ax.set_ylabel(r'$\it{S.\ pastorianus}$ strain')   # mathtext so only the species name is italic
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title('PAD1 / FDC1 status by lager subgenome', fontsize=7)
ax.text(-0.32, 1.08, 'a', transform=ax.transAxes, fontsize=8, fontweight='bold')
plt.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'{FIG_DIR}/ED_lager_subgenome.{ext}', dpi=300)
plt.show()
