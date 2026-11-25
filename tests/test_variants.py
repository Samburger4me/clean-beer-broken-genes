"""Self-check for src/yeastflavor/03_variant_calling.py on a mock allele set with known variants.
Run from the repo root: python -m pytest tests/  (or python tests/test_variants.py)"""
import os, shutil, subprocess, sys, tempfile
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tmp = tempfile.mkdtemp(prefix='vc_mock_')
D = f'{tmp}/data'
for sub in ('reference', 'metadata', 'seqs'):
    os.makedirs(f'{D}/{sub}')
for g in ('PAD1', 'FDC1'):
    shutil.copy(f'{ROOT}/data/reference/{g}_cds.fasta', f'{D}/reference/')
shutil.copy(f'{ROOT}/data/metadata/strain_lookup.csv', f'{D}/metadata/')
ref = {g: str(next(SeqIO.parse(f'{D}/reference/{g}_cds.fasta', 'fasta')).seq) for g in ('PAD1', 'FDC1')}

def sub(s, pos, b):  # 1-based substitution
    assert s[pos - 1] != b
    return s[:pos - 1] + b + s[pos:]

F = ref['FDC1']
assert F[459:462] == 'CAA'
# find a homopolymer run (>=3) in FDC1 around codon ~166 to test dup + 3' shifting
hp = next(i for i in range(480, 1400) if F[i] == F[i + 1] == F[i + 2] and F[i - 1] != F[i] and F[i + 3] != F[i])
hp_base, hp_last = F[hp], hp + 3          # run occupies 1-based hp+1..hp+3
codon10 = F[27:30]
mis = next(b for b in 'ACGT' if b != F[27] and str(Seq(b + F[28:30]).translate()) not in ('*', str(Seq(codon10).translate())))
c20 = F[57:60]  # third-position change: find one that is synonymous
syn = next((b for b in 'ACGT' if b != F[59] and Seq(F[57:59] + b).translate() == Seq(c20).translate()), None)
assert syn, 'pick another codon'
alleles = {
    'MOCK_Q154.1': ('wine yeast Kyokai No.6 FDC1 gene', sub(F, 460, 'T')),
    'MOCK_FS.1': ('Saccharomyces pastorianus FDC1 gene, strain: NCYC984', F[:hp + 1] + hp_base + F[hp + 1:]),  # +1 in run
    'MOCK_MIS.1': ('Saccharomyces cerevisiae EC1118 chromosome IV', sub(sub(F, 28, mis), 60, syn)),
    'MOCK_PART.1': ('Saccharomyces cerevisiae strain YPS128 FDC1 gene, partial cds PARTIAL', sub(F, 460, 'T')[300:1200]),
    'MOCK_DEL3.1': ('Saccharomyces cerevisiae strain Foo99 FDC1 gene', F[:600] + F[603:]),
    'MOCK_GAL.1': ('Saccharomyces cerevisiae strain beer042 chromosome IV, whole genome shotgun sequence', sub(F, 460, 'T')),
}
with open(f'{D}/seqs/FDC1_alleles.fasta', 'w') as fh:
    fh.write(f'>S288C_REF 1-1512 strand=+ | S288C\n{F}\n')
    for k, (t, s) in alleles.items():
        fh.write(f'>{k} 1-{len(s)} strand=+ | {t}\n{s}\n')
    fh.write(f'>MOCK_N.1 1-1512 strand=+ CONTAINS_N | Saccharomyces cerevisiae strain wild005 chromosome IV\n{sub(F, 460, "T")[:700]}NNNN{F[704:]}\n')
P = ref['PAD1']
pad_stop = next(i for i in range(30, 400, 3) if P[i:i + 3] in ('CAA', 'CAG', 'CGA'))  # C->T gives stop
with open(f'{D}/seqs/PAD1_alleles.fasta', 'w') as fh:
    fh.write(f'>S288C_REF | S288C\n{P}\n>MOCK_PADX.1 | Saccharomyces cerevisiae NCYC240 PAD1 gene\n{sub(P, pad_stop + 1, "T")}\n'
             f'>MOCK_PADOK.1 | Saccharomyces cerevisiae EC1118 chromosome IV\n{P}\n')

env = dict(os.environ, VC_DATA=D, VC_FIG=f'{tmp}/fig', MPLBACKEND='Agg')
r = subprocess.run([sys.executable, f'{ROOT}/src/yeastflavor/03_variant_calling.py'], env=env, capture_output=True, text=True)
assert r.returncode == 0, r.stderr[-3000:]
V = pd.read_csv(f'{D}/results/variant_calls.csv')
S = pd.read_csv(f'{D}/results/allele_summary.csv').set_index('sacc')
M = pd.read_csv(f'{D}/metadata/strains.csv').set_index('sacc')

def calls(sacc):
    return sorted(zip(V[V.sacc == sacc].hgvs_c, V[V.sacc == sacc].effect, V[V.sacc == sacc].protein_change.fillna('')))

assert calls('MOCK_Q154.1') == [('c.460C>T', 'nonsense', 'p.Q154*')], calls('MOCK_Q154.1')
assert S.loc['MOCK_Q154.1', 'protein_len'] == 153 and S.loc['MOCK_Q154.1', 'LOF']
fsc = calls('MOCK_FS.1')
assert fsc == [(f'c.{hp_last}dup{hp_base}', 'frameshift', fsc[0][2])] and fsc[0][2].endswith('fs'), fsc
assert S.loc['MOCK_FS.1', 'frameshift'] and S.loc['MOCK_FS.1', 'LOF']
mc = calls('MOCK_MIS.1')
assert [c[:2] for c in mc] == [(f'c.28{F[27]}>{mis}', 'missense'), (f'c.60{F[59]}>{syn}', 'synonymous')], mc
assert mc[0][2].startswith('p.') and mc[0][2][2:].startswith(str(Seq(codon10).translate()) + '10')
assert S.loc['MOCK_MIS.1', 'LOF'] == 0 and S.loc['MOCK_MIS.1', 'protein_len'] == 503
assert calls('MOCK_PART.1') == [('c.460C>T', 'nonsense', 'p.Q154*')]       # no end-gap 'deletions'
assert S.loc['MOCK_PART.1', 'ref_covered'] == '301-1200' and S.loc['MOCK_PART.1', 'LOF']
assert S.loc['MOCK_PART.1', 'protein_len'] == 153
d3 = calls('MOCK_DEL3.1')
assert len(d3) == 1 and d3[0][1] == 'inframe_deletion' and 'del' in d3[0][0], d3
assert not S.loc['MOCK_DEL3.1', 'frameshift'] and S.loc['MOCK_DEL3.1', 'LOF'] == 0 and S.loc['MOCK_DEL3.1', 'protein_len'] == 502
assert calls('MOCK_PADX.1')[0][1] == 'nonsense' and S.loc['MOCK_PADX.1', 'protein_len'] == pad_stop // 3
assert calls('MOCK_PADOK.1') == [] and S.loc['MOCK_PADOK.1', 'pct_identity'] == 100
# strain parsing
assert M.loc['MOCK_FS.1', 'strain'] == 'NCYC984' and M.loc['MOCK_FS.1', 'group'] == 'Lager'
assert M.loc['MOCK_MIS.1', 'strain'] == 'EC1118' and M.loc['MOCK_Q154.1', 'strain'] == 'Kyokai wine yeast No.6' or \
       M.loc['MOCK_Q154.1', 'group'] == 'Other/Unknown'  # 'wine yeast No.6' is not in lookup -> unknown
assert M.loc['MOCK_PADX.1', 'group'] == 'Ale' and M.loc['MOCK_PART.1', 'strain'] == 'YPS128'
assert M.loc['MOCK_DEL3.1', 'group'] == 'Other/Unknown'
assert 'Q154' in str(M.loc['MOCK_Q154.1', 'notes']) or M.loc['MOCK_Q154.1', 'group'] != 'Other/Unknown'
G = pd.read_csv(f'{D}/results/pof_genotype_by_strain.csv').set_index('strain')
assert M.loc['MOCK_GAL.1', 'group'] == 'Beer (Gallone)' and M.loc['MOCK_GAL.1', 'dataset'] == 'Gallone2016'
assert M.loc['MOCK_FS.1', 'dataset'] == 'NCBI_core_nt'
assert S.loc['MOCK_N.1', 'LOF_call'] == 'not callable' and pd.isna(S.loc['MOCK_N.1', 'LOF'])
assert G.loc['wild005', 'predicted_POF'] == '?' and G.loc['beer042', 'predicted_POF'] == 'POF-'
R = pd.read_csv(f'{D}/results/recurrent_lof_alleles.csv').set_index('hgvs_c')
assert R.loc['c.460C>T', 'n_strains'] == 3  # Q154 mock, PART (flag only in title, so callable), beer042; MOCK_N excluded
assert os.path.exists(f'{D}/results/lof_enrichment_tests.csv')
assert G.loc['EC1118', 'predicted_POF'] == 'POF+' and G.loc['NCYC240', 'predicted_POF'] == 'POF-'
for f in ('03a_variant_tracks', '03b_lof_alignment_windows', '03c_lof_fraction_by_group', '03d_lof_matrix'):
    assert os.path.exists(f'{tmp}/fig/{f}.png'), f
shutil.rmtree(tmp)
print('test_variants: all checks passed')


def test_module_checks_passed():
    """The asserts above run at import; reaching here means they all passed (lets pytest report it)."""
    assert True
