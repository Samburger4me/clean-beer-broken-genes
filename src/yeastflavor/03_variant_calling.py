# %%
# 3. Variant calling in PAD1 and FDC1 alleles
import os, re
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data.CodonTable import standard_dna_table

DATA = os.environ.get('VC_DATA', 'data')      # VC_* env vars only used by tests/test_variants.py
FIG = os.environ.get('VC_FIG', 'figures')
RES = f'{DATA}/results'
for d in (RES, FIG):
    os.makedirs(d, exist_ok=True)
VC_GENES = ['PAD1', 'FDC1']
GROUP_ORDER = ['Laboratory', 'Wild', 'Wine', 'Clinical', 'Bioethanol', 'Baking', 'Wheat/Hefeweizen',
               'Beer (Gallone)', 'Ale', 'Lager', 'Sake', 'Shochu/Awamori', 'Spirits', 'Other/Unknown']
ref_cds = {g: str(next(SeqIO.parse(f'{DATA}/reference/{g}_cds.fasta', 'fasta')).seq).upper() for g in VC_GENES}
ref_prot = {g: str(Seq(s).translate(to_stop=True)) for g, s in ref_cds.items()}
alleles = {g: [r for r in SeqIO.parse(f'{DATA}/seqs/{g}_alleles.fasta', 'fasta')
               if not r.id.startswith('S288C_REF')] for g in VC_GENES}
for g in VC_GENES:
    print(f'{g}: reference {len(ref_cds[g])} nt / {len(ref_prot[g])} aa; {len(alleles[g])} alleles')

# %%
# 3.1 Strain metadata
lookup = pd.read_csv(f'{DATA}/metadata/strain_lookup.csv')

def norm_strain(s):
    s = re.sub(r'\bno\b\.?', '', str(s).lower())
    s = re.sub(r'[\s\-_.:/\'’,()]', '', s)
    # DDBJ (Mukai) titles prefix the niche: 'Wine Yeast Kyokai No. 1', 'Top Fermentation Yeast NCYC240'
    s = s.replace('wineyeastkyokai', 'kyokaiwineyeast').replace('shochuyeastkyokai', 'shochukyokai')
    s = re.sub(r'^(?:top|bottom)fermentationyeast|^(?:sake|shochu|wine)yeast', '', s)
    return s.replace('kagishima', 'kagoshima')  # GenBank typo, see strain_lookup notes

lookup_by_norm = {norm_strain(s): s for s in lookup.strain}
STOPWORDS = r'(?:chromosome|genomic|DNA|PAD1|FDC1|YDR53[89]W|gene|genes|complete|mitochondri\w*|partial|cds|sequence|strain|isolate|clone|contig|scaffold)'

def parse_strain(title, sacc):
    t = re.sub(r'^\S+\s+', '', title) if re.match(r'^[A-Z]{1,2}_?\d+\.\d\s', title) else title
    for pat in [r'strain[:\s]+', r'isolate[:\s]+', r'Saccharomyces (?:cerevisiae|pastorianus|carlsbergensis)\s+(?:x\s+\S+\s+)?']:
        m = re.search(pat + r'(.+?)(?=\s+' + STOPWORDS + r'\b|[,;]|$)', t, flags=re.I)
        name = m.group(1).strip() if m else ''
        if name and (pat.startswith(('strain', 'isolate')) or (len(name.split()) <= 2 and re.search(r'\d', name))):
            return name  # bare 'S. cerevisiae X ...' only counts if X looks like a strain code (N85, Y7)
    return None

def match_strain(name):
    if name is None:
        return None
    n = norm_strain(name)
    if n in lookup_by_norm:
        return lookup_by_norm[n]
    hits = [k for k in lookup_by_norm if len(k) >= 4 and (k in n or n in k)]  # e.g. 'NCYC 240 haploid'
    return lookup_by_norm[max(hits, key=len)] if hits else None

titles = {}
for g in VC_GENES:
    for r in alleles[g]:
        titles[r.id] = r.description.split('|', 1)[1].strip() if '|' in r.description else r.description
GALLONE = {'beer': 'Beer (Gallone)', 'wine': 'Wine', 'sake': 'Sake', 'spirits': 'Spirits', 'wild': 'Wild',
           'bioethanol': 'Bioethanol', 'bread': 'Baking', 'laboratory': 'Laboratory'}
GALLONE_CIT = 'Gallone et al. 2016 Cell PMID:27610566 (BioProject PRJNA323691)'
S288C_RECORDS = {'U43834.1'}  # yeast genome project chromosome IV clone (S288C); NM_/NC_ = RefSeq S288C
rows = []
for sacc, title in titles.items():
    base = re.sub(r'(\.\d+)_\d+$', r'\1', sacc)
    raw = parse_strain(title, base)
    gm = re.fullmatch(r'([a-z]+)\d{3}', raw or '')
    hit = None if gm else match_strain(raw)
    if gm and gm.group(1) in GALLONE:
        rows.append(dict(sacc=sacc, strain=raw, species='S. cerevisiae', group=GALLONE[gm.group(1)],
                         source_niche=f'{gm.group(1)} (Gallone 2016 niche label)', citation=GALLONE_CIT, notes='',
                         dataset='Gallone2016', parsed_name=raw, title=title))
    elif base.startswith(('NM_', 'NC_')) or base in S288C_RECORDS:
        L = lookup[lookup.strain == 'S288C'].iloc[0]
        rows.append(dict(sacc=sacc, strain=f'S288C ({sacc})', species=L.species, group='Laboratory',
                         source_niche=L.source_niche, citation=L.citation, notes='S288C-derived reference record',
                         parsed_name=raw, title=title))
    elif hit:
        L = lookup[lookup.strain == hit].iloc[0]
        rows.append(dict(sacc=sacc, strain=hit, species=L.species, group=L.group, source_niche=L.source_niche,
                         citation=L.citation, notes=L.notes, parsed_name=raw, title=title))
    else:
        sp = 'S. pastorianus' if re.search(r'pastorianus|carlsbergensis', title, re.I) else 'S. cerevisiae'
        rows.append(dict(sacc=sacc, strain=raw or f'unnamed ({sacc})', species=sp, group='Other/Unknown', source_niche='',
                         citation='', notes='strain not in lookup table', parsed_name=raw, title=title))
strains = pd.DataFrame(rows)
strains['dataset'] = strains.reindex(columns=['dataset'])['dataset'].fillna('NCBI_core_nt')
strains.drop(columns=['parsed_name', 'title']).to_csv(f'{DATA}/metadata/strains.csv', index=False)
pd.set_option('display.width', 200, 'display.max_colwidth', 60)
print(pd.crosstab(strains.group, strains.dataset, margins=True).to_string())
print('\nNCBI core_nt records:')
print(strains[strains.dataset == 'NCBI_core_nt'][['sacc', 'strain', 'species', 'group']].to_string(index=False))

# %%
# 3.2 Pairwise alignment and variant calling
nt_aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-3,
                             open_gap_score=-10, extend_gap_score=-1)
nt_aligner.target_end_gap_score = nt_aligner.query_end_gap_score = 0.0
prot_aligner = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                               open_gap_score=-11, extend_gap_score=-1)
CODON = dict(standard_dna_table.forward_table, **{c: '*' for c in standard_dna_table.stop_codons})

def aa_of(codon):
    return CODON.get(codon, 'X')

def right_shift_ins(ref, p, ins):
    """Insertion of `ins` after 1-based ref position p -> 3'-most equivalent position (HGVS)."""
    while p < len(ref) and ref[p] == ins[0]:
        ins, p = ins[1:] + ins[0], p + 1
    return p, ins

def right_shift_del(ref, a, b):
    """Deletion of ref[a:b] (0-based) -> 3'-most equivalent."""
    while b < len(ref) and ref[b] == ref[a]:
        a, b = a + 1, b + 1
    return a, b

def call_allele(gene, rec):
    ref, alt = ref_cds[gene], str(rec.seq).upper()
    aln = nt_aligner.align(ref, alt)[0]
    R, A = aln[0], aln[1]
    both = [i for i in range(len(R)) if R[i] != '-' and A[i] != '-']
    i0, i1 = both[0], both[-1]
    rpos, apos = i0 - R[:i0].count('-'), i0 - A[:i0].count('-')
    r_start, a_start = rpos, apos            # 0-based first covered ref / allele positions
    variants, snv_alt, ident, i = [], {}, 0, i0
    while i <= i1:
        if R[i] != '-' and A[i] != '-':
            if R[i] == A[i]:
                ident += 1
            elif A[i] in 'ACGT':
                snv_alt[rpos] = A[i]
            else:
                variants.append(dict(cds_pos=rpos + 1, ref=R[i], alt=A[i], type='ambiguous', effect='unknown',
                                     protein_change='', hgvs_c=f'c.{rpos + 1}{R[i]}>{A[i]}'))
            rpos, apos, i = rpos + 1, apos + 1, i + 1
        else:
            j = i
            while j <= i1 and (R[j] == '-') == (R[i] == '-') and (R[j] == '-' or A[j] == '-'):
                j += 1
            if R[i] == '-':                                  # insertion after ref position rpos (1-based)
                p, ins = right_shift_ins(ref, rpos, A[i:j])
                n = len(ins)
                if ref[p - n:p] == ins:
                    kind, h = 'duplication', (f'c.{p}dup{ins}' if n == 1 else f'c.{p - n + 1}_{p}dup{ins}')
                    pos = p - n + 1
                else:
                    kind, h, pos = 'insertion', f'c.{p}_{p + 1}ins{ins}', p
                codon = p // 3 + 1 if kind == 'insertion' else (p - n) // 3 + 1
                variants.append(dict(cds_pos=pos, ref='-', alt=ins, type=kind, indel_len=n, hgvs_c=h, codon=codon))
                apos += j - i
            else:                                            # deletion of ref[rpos:rpos+len]
                a, b = right_shift_del(ref, rpos, rpos + (j - i))
                h = f'c.{a + 1}del{ref[a:b]}' if b - a == 1 else f'c.{a + 1}_{b}del{ref[a:b]}'
                variants.append(dict(cds_pos=a + 1, ref=ref[a:b], alt='-', type='deletion', indel_len=-(b - a),
                                     hgvs_c=h, codon=a // 3 + 1))
                rpos += j - i
            i = j
    r_end = rpos                                                  # 0-based exclusive
    for p, b in snv_alt.items():                                  # SNV effect with all SNVs of that codon applied
        c0 = p - p % 3
        rc = ref[c0:c0 + 3]
        mc = ''.join(snv_alt.get(c0 + k, rc[k]) for k in range(3))
        ra, ma, cn = aa_of(rc), aa_of(mc), c0 // 3 + 1
        if c0 == 0 and rc == 'ATG' and mc != 'ATG':
            eff = 'start_lost'
        elif ra == ma:
            eff = 'synonymous'
        elif ma == '*':
            eff = 'nonsense'
        elif ra == '*':
            eff = 'stop_lost'
        else:
            eff = 'missense'
        variants.append(dict(cds_pos=p + 1, ref=ref[p], alt=b, type='SNV', effect=eff,
                             protein_change='' if eff == 'synonymous' else f'p.{ra}{cn}{ma}',
                             hgvs_c=f'c.{p + 1}{ref[p]}>{b}'))
    net = sum(v.get('indel_len', 0) for v in variants)
    frameshift = net % 3 != 0
    for v in variants:
        if 'indel_len' in v:
            fs = v['indel_len'] % 3 != 0
            v['effect'] = 'frameshift' if fs else ('inframe_' + ('insertion' if v['indel_len'] > 0 else 'deletion'))
            raa = ref_prot[gene][v['codon'] - 1] if v['codon'] <= len(ref_prot[gene]) else '*'
            v['protein_change'] = f"p.{raa}{v['codon']}{'fs' if fs else ('ins' if v['indel_len'] > 0 else 'del')}"
    # translate allele in the reference frame, starting at the first full codon it covers
    k = (-r_start) % 3
    s = alt[a_start + k:]
    s = s[:len(s) // 3 * 3]
    aa_offset = (r_start + k) // 3
    prot = str(Seq(s).translate()) if s else ''
    stop = prot.find('*')
    has_stop = stop >= 0
    protein = prot[:stop] if has_stop else prot
    protein_len = aa_offset + len(protein)
    reflen = len(ref_prot[gene])
    premature = has_stop and protein_len < 0.95 * reflen
    lof = bool(premature or frameshift)
    lof_allele = ''
    if lof:
        causal = sorted([v for v in variants if v['effect'] in ('frameshift', 'nonsense')], key=lambda v: v['cds_pos'])
        what = '; '.join(f"{v['hgvs_c']} {v['protein_change']}" for v in causal) or 'premature stop'
        lof_allele = f"{what} -> {protein_len} aa" + ('' if has_stop else ' (no stop within record)')
    pa = prot_aligner.align(ref_prot[gene], protein)[0] if protein else None
    prot_id = (sum(x == y for x, y in zip(pa[0], pa[1]) if x != '-' and y != '-') / len(protein) * 100) if pa else np.nan
    for v in variants:
        v.pop('codon', None)
    flags = ' '.join(f for f in ('PARTIAL', 'CONTAINS_N') if f in rec.description.split('|')[0])
    callable_ = not flags  # partial or N-containing records cannot support an 'intact' or 'LOF' call
    summary = dict(length_nt=len(alt), ref_covered=f'{r_start + 1}-{r_end}', pct_identity=round(100 * ident / (i1 - i0 + 1), 3),
                   n_snv=sum(v['type'] == 'SNV' for v in variants),
                   n_missense=sum(v['effect'] == 'missense' for v in variants),
                   n_synonymous=sum(v['effect'] == 'synonymous' for v in variants),
                   n_nonsense=sum(v['effect'] == 'nonsense' for v in variants),
                   n_indel=sum('indel_len' in v for v in variants), frameshift=frameshift,
                   protein_len=protein_len, stop_found=has_stop, protein_pct_identity=round(prot_id, 2),
                   flags=flags, callable=callable_)
    summary.update(LOF=lof if callable_ else np.nan, LOF_allele=lof_allele if callable_ else '',
                   LOF_call=('LOF' if lof else 'intact') if callable_ else 'not callable')
    return variants, summary, protein, aln

meta = strains.set_index('sacc')
var_rows, sum_rows, ALN = [], [], {}
for g in VC_GENES:
    prot_recs = []
    for rec in alleles[g]:
        vs, sm, protein, aln = call_allele(g, rec)
        ALN[(g, rec.id)] = aln
        st, gr = meta.loc[rec.id, 'strain'], meta.loc[rec.id, 'group']
        for v in vs:
            v.pop('indel_len', None)
            var_rows.append(dict(gene=g, sacc=rec.id, strain=st, group=gr, **v))
        sum_rows.append(dict(gene=g, sacc=rec.id, strain=st, group=gr, **sm))
        prot_recs.append(f'>{rec.id} {st} | {gr} | {sm["protein_len"]} aa | {sm["LOF_call"]}\n{protein}\n')
    with open(f'{RES}/protein_alleles_{g}.fasta', 'w') as fh:
        fh.write(f'>S288C_REF {g} {len(ref_prot[g])} aa\n{ref_prot[g]}\n' + ''.join(prot_recs))
variant_calls = pd.DataFrame(var_rows, columns=['gene', 'sacc', 'strain', 'group', 'cds_pos', 'ref', 'alt', 'type',
                                                'effect', 'protein_change', 'hgvs_c']).sort_values(['gene', 'sacc', 'cds_pos'])
allele_summary = pd.DataFrame(sum_rows)
variant_calls.to_csv(f'{RES}/variant_calls.csv', index=False)
KNOWN_LOF = {('FDC1', 'c.160A>T'): 'carries sake-type K54*', ('FDC1', 'c.460C>T'): 'carries Beer 1 ale-type Q154*',
             ('FDC1', 'c.501dupA'): 'carries ale/lager FDC1 c.501dupA frameshift', ('PAD1', 'c.305G>A'): 'carries ale-type PAD1 W102*'}
for (g, h), note in KNOWN_LOF.items():
    carriers = set(variant_calls[(variant_calls.gene == g) & (variant_calls.hgvs_c == h)].sacc)
    m = strains.sacc.isin(carriers) & (strains.group == 'Other/Unknown')
    strains.loc[m, 'notes'] = strains.loc[m, 'notes'].fillna('').str.cat([note] * m.sum(), sep='; ').str.strip('; ')
for r in allele_summary[allele_summary.pct_identity < 95].itertuples():  # e.g. introgressed alleles
    strains.loc[strains.sacc == r.sacc, 'notes'] = (strains.loc[strains.sacc == r.sacc, 'notes'].fillna('') +
                                                    f'; divergent {r.gene} allele ({r.pct_identity:.1f}% nt identity to S288C) - possible introgression').str.strip('; ')
strains.drop(columns=['parsed_name', 'title']).to_csv(f'{DATA}/metadata/strains.csv', index=False)
allele_summary.to_csv(f'{RES}/allele_summary.csv', index=False)
print(allele_summary.groupby(['gene', 'LOF_call']).size().unstack(fill_value=0).to_string())
print('\nNot callable:', allele_summary[~allele_summary.callable][['gene', 'sacc', 'strain', 'flags']].to_string(index=False))
print('\nDivergent alleles (<95 % nt identity):')
print(allele_summary[allele_summary.pct_identity < 95][['gene', 'sacc', 'strain', 'group', 'pct_identity', 'LOF_call']].to_string(index=False))

# %%
# Variant counts by effect (full list in data/results/variant_calls.csv)
print(variant_calls.groupby(['gene', 'effect']).size().unstack(fill_value=0).to_string())

# %%
# 3.3 Variant tracks along each CDS, aggregated by group
EFFECT_COL = {'synonymous': '#b0b0b0', 'missense': '#4c78a8', 'inframe_insertion': '#72b7b2',
              'inframe_deletion': '#72b7b2', 'nonsense': '#e45756', 'frameshift': '#b279a2',
              'start_lost': '#f58518', 'stop_lost': '#f58518', 'unknown': '#eeca3b'}
LOFEFF = ['nonsense', 'frameshift']
AS_c = allele_summary[allele_summary.callable]
VC = variant_calls[variant_calls.sacc.isin(AS_c.sacc)].merge(AS_c[['gene', 'sacc']], on=['gene', 'sacc'])
grp_order = [g for g in GROUP_ORDER if g in set(allele_summary.group)]
n_grp = AS_c.groupby(['gene', 'group']).size()
freq = VC.groupby(['gene', 'group', 'hgvs_c', 'cds_pos', 'effect', 'protein_change']).sacc.nunique().rename('n').reset_index()
freq['frac'] = freq.n / [n_grp[(g, gr)] for g, gr in zip(freq.gene, freq.group)]
n_tot = VC.groupby(['gene', 'hgvs_c']).sacc.nunique()
freq['n_total'] = [n_tot[(g, h)] for g, h in zip(freq.gene, freq.hgvs_c)]
shown = freq[freq.effect.isin(LOFEFF) | ((freq.effect != 'synonymous') & (freq.n_total > 1))]

fig, axes = plt.subplots(1, 2, figsize=(17, 0.5 * len(grp_order) + 3))
for ax, g in zip(axes, VC_GENES):
    S = shown[shown.gene == g]
    ypos = {gr: i for i, gr in enumerate(grp_order)}
    for x in S[S.effect.isin(LOFEFF)].drop_duplicates('hgvs_c').itertuples():
        ax.axvline(x.cds_pos, color=EFFECT_COL[x.effect], lw=0.6, ls=':', zorder=0)
        ax.text(x.cds_pos, -0.8, f'{x.protein_change}', rotation=60, fontsize=7, ha='left', va='bottom')
    for eff, sub in S.groupby('effect'):
        ax.scatter(sub.cds_pos, sub.group.map(ypos), s=15 + 400 * sub.frac, color=EFFECT_COL.get(eff, 'k'),
                   marker='X' if eff in LOFEFF else 'o', alpha=0.85, edgecolor='k', linewidth=0.3, label=eff, zorder=3)
    ax.set_yticks(range(len(grp_order)))
    ax.set_yticklabels([f'{gr} (n={n_grp.get((g, gr), 0)})' for gr in grp_order], fontsize=8)
    ax.set_ylim(len(grp_order) - 0.5, -1.0)
    ax.set_xlim(0, len(ref_cds[g]) + 1)
    ax.set_xlabel(f'{g} CDS position (nt, S288C coordinates; {len(ref_cds[g])} nt)')
    ax.set_title(f'{g}: recurrent / LOF variants by group (mark area = fraction of group)', pad=45)
    ax.grid(axis='y', color='#eeeeee')
    ax.legend(fontsize=7, loc='lower right', framealpha=0.9, markerscale=0.6)
plt.tight_layout()
plt.savefig(f'{FIG}/03a_variant_tracks.png', dpi=200)
plt.show()

# %%
# 3.4 Recurrent LOF alleles and alignment windows
BASE_COL = {'A': '#59a14f', 'C': '#4e79a7', 'G': '#f28e2b', 'T': '#e15759', '-': '#ffffff'}

def ref_projection(aln):
    """dict ref_pos(0-based) -> (base, inserted_after) for one allele alignment."""
    R, A = aln[0], aln[1]
    out, rpos = {}, -1
    for x, y in zip(R, A):
        if x == '-':
            if rpos >= 0:
                b, ins = out[rpos]
                out[rpos] = (b, ins + y)
        else:
            rpos += 1
            out[rpos] = (y, '')
    return out

def window_plot(ax, g, center, sacc_list, half=12):
    lo, hi = max(0, center - 1 - half), min(len(ref_cds[g]), center + half)
    proj = {s: ref_projection(ALN[(g, s)]) for s in sacc_list}
    max_ins = {p: max([len(proj[s].get(p, ('', ''))[1]) for s in sacc_list] + [0]) for p in range(lo, hi)}
    cols = []  # (ref_pos, ins_index or None)
    for p in range(lo, hi):
        cols.append((p, None))
        cols += [(p, k) for k in range(max_ins[p])]
    rows = [('S288C', None)] + [(s, proj[s]) for s in sacc_list]
    for y, (name, pr) in enumerate(rows):
        for x, (p, k) in enumerate(cols):
            if pr is None:
                b = ref_cds[g][p] if k is None else '-'
                refb = b
            else:
                b0, ins = pr.get(p, (' ', ''))
                b = b0 if k is None else (ins[k] if k < len(ins) else '-')
                refb = ref_cds[g][p] if k is None else '-'
            same = b == refb
            ax.add_patch(plt.Rectangle((x, y), 1, 1, color=BASE_COL.get(b, '#dddddd'), alpha=0.25 if same and pr is not None else 1))
            ax.text(x + 0.5, y + 0.55, b, ha='center', va='center', fontsize=6.5, fontweight='normal' if same else 'bold')
    ticks = [x for x, (p, k) in enumerate(cols) if k is None and p % 3 == 0]
    ax.set_xticks([t + 0.05 for t in ticks])
    ax.set_xticklabels([f'c.{cols[t][0] + 1}\n(codon {cols[t][0] // 3 + 1})' for t in ticks], fontsize=6)
    ax.set_yticks([y + 0.5 for y in range(len(rows))])
    ax.set_yticklabels(['S288C'] + [f"{meta.loc[s, 'strain']} [{meta.loc[s, 'group']}]" for s in sacc_list], fontsize=6.5)
    ax.set_xlim(0, len(cols)); ax.set_ylim(len(rows), 0)
    for x, (p, k) in enumerate(cols):
        if (k is None and p == center - 1) or k is not None:  # mark the site and any insertion column
            ax.add_patch(plt.Rectangle((x, 0), 1, len(rows), fill=False, ec='k', lw=1.2))

lof_v = VC[VC.effect.isin(LOFEFF)].merge(AS_c[['gene', 'sacc', 'LOF']], on=['gene', 'sacc'])
recurrent = (lof_v.groupby(['gene', 'hgvs_c', 'cds_pos', 'effect', 'protein_change'])
             .agg(n_strains=('strain', 'nunique'), n_records=('sacc', 'nunique'), n_called_LOF=('LOF', lambda x: int(np.nansum(x))),
                  groups=('group', lambda x: '; '.join(f'{k} {v}' for k, v in x.value_counts().items())),
                  carriers=('sacc', list))
             .reset_index().sort_values(['n_strains', 'gene'], ascending=[False, True]))
recurrent.drop(columns='carriers').to_csv(f'{RES}/recurrent_lof_alleles.csv', index=False)
print(recurrent.drop(columns=['carriers', 'cds_pos']).to_string(index=False))

def pick(saccs, k):
    """Up to k records, spread across groups (round-robin)."""
    by = meta.loc[saccs].groupby('group').groups
    lists, out = [list(v) for v in by.values()], []
    while len(out) < k and any(lists):
        for l in lists:
            if l and len(out) < k:
                out.append(l.pop(0))
    return out

if len(recurrent):
    show = recurrent.head(6)
    fig, axes = plt.subplots(len(show), 1, figsize=(12, 2.9 * len(show)), squeeze=False)
    for ax, r in zip(axes[:, 0], show.itertuples()):
        others = pick([s for s in AS_c[AS_c.gene == r.gene].sacc if s not in r.carriers], 4)
        window_plot(ax, r.gene, r.cds_pos, pick(r.carriers, 9) + others)
        ax.set_title(f'{r.gene} {r.hgvs_c} {r.protein_change} - {r.n_strains} strain(s): {r.groups}', fontsize=8)
    plt.tight_layout()
    plt.savefig(f'{FIG}/03b_lof_alignment_windows.png', dpi=200)
    plt.show()
else:
    print('No nonsense/frameshift variants found - no alignment windows drawn.')

# %%
# 3.5 LOF frequency by ecological group and the predicted POF phenotype
lof_mat = allele_summary.pivot_table(index=['group', 'strain'], columns='gene', values='LOF', aggfunc='mean').reindex(columns=VC_GENES)
for st_ in set(allele_summary.set_index(['group', 'strain']).index) - set(lof_mat.index):  # only non-callable records
    lof_mat.loc[st_, :] = np.nan

def pof_call(r):
    if (r == 1).any():
        return 'POF-'
    if r.isna().any() or (r > 0).any():
        return '?'
    return 'POF+'

pof = lof_mat.apply(pof_call, axis=1).rename('predicted_POF')
genotype = lof_mat.assign(predicted_POF=pof).reset_index()
genotype['g_rank'] = genotype.group.map({g: i for i, g in enumerate(GROUP_ORDER)}).fillna(99)
genotype = genotype.sort_values(['g_rank', 'strain']).drop(columns='g_rank').reset_index(drop=True)
genotype.to_csv(f'{RES}/pof_genotype_by_strain.csv', index=False)

frac = (AS_c.groupby(['group', 'gene']).LOF.agg(['mean', 'size']).reset_index())
gp = genotype[genotype.predicted_POF != '?'].groupby('group').predicted_POF
frac = pd.concat([frac, pd.DataFrame({'group': gp.size().index, 'gene': 'either (strain POF-)',
                                      'mean': gp.apply(lambda x: (x == 'POF-').mean()).values, 'size': gp.size().values})])
hue_order = VC_GENES + ['either (strain POF-)']
frac = (frac.set_index(['group', 'gene']).reindex(pd.MultiIndex.from_product([grp_order, hue_order], names=['group', 'gene']))
        .fillna({'mean': 0, 'size': 0}).reset_index())
fig, ax = plt.subplots(figsize=(12, 4.8))
sns.barplot(data=frac, x='group', y='mean', hue='gene', order=grp_order, hue_order=hue_order,
            palette=['#4c78a8', '#e45756', '#555555'], ax=ax)
for cont, g in zip(ax.containers, hue_order):
    sub = frac[frac.gene == g].set_index('group')['size']
    ax.bar_label(cont, labels=[f'{int(sub[t.get_text()])}' for t in ax.get_xticklabels()], fontsize=6.5)
ax.set_ylim(0, 1.12); ax.set_ylabel('fraction LOF (callable records / strains)'); ax.set_xlabel('')
ax.set_title('Predicted loss-of-function PAD1 / FDC1 by ecological group (numbers above bars = n)')
plt.xticks(rotation=30, ha='right'); plt.tight_layout()
plt.savefig(f'{FIG}/03c_lof_fraction_by_group.png', dpi=200)
plt.show()
print(frac.pivot(index='group', columns='gene', values=['mean', 'size']).reindex(grp_order).round(3).to_string())

# %%
def geno_class(r):
    p, f = r['PAD1'], r['FDC1']
    if p == 1 and f == 1: return 'both LOF'
    if f == 1: return 'FDC1 LOF only'
    if p == 1: return 'PAD1 LOF only'
    if p == 0 and f == 0: return 'both intact'
    return 'incomplete / not callable'

CLASS_COL = {'both intact': '#93c47d', 'PAD1 LOF only': '#6fa8dc', 'FDC1 LOF only': '#e06666', 'both LOF': '#990000',
             'incomplete / not callable': '#d9d9d9'}
genotype['genotype_class'] = genotype.apply(geno_class, axis=1)
genotype.to_csv(f'{RES}/pof_genotype_by_strain.csv', index=False)
counts = pd.crosstab(genotype.group, genotype.genotype_class).reindex(index=grp_order, columns=list(CLASS_COL), fill_value=0)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, max(5, 0.42 * len(grp_order) + 2)), gridspec_kw={'width_ratios': [3, 1.2]})
(counts.div(counts.sum(axis=1), axis=0)).plot.barh(stacked=True, color=list(CLASS_COL.values()), ax=a1, width=0.8, edgecolor='white')
a1.set_yticklabels([f'{g} (n={counts.loc[g].sum()})' for g in grp_order], fontsize=8); a1.invert_yaxis()
a1.set_xlabel('fraction of strains'); a1.set_ylabel(''); a1.legend(fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=5, frameon=False)
a1.set_title('PAD1 / FDC1 genotype class per group')
hm = genotype[VC_GENES].astype(float)
sns.heatmap(hm, cmap=sns.color_palette(['#93c47d', '#f6b26b', '#cc0000'], as_cmap=True), vmin=0, vmax=1, ax=a2,
            cbar_kws={'label': '0 intact / 1 LOF (grey = no call)'}, yticklabels=False)
a2.set_facecolor('#d9d9d9')
bounds = genotype.index[genotype.group.ne(genotype.group.shift())].tolist() + [len(genotype)]
for b0, b1 in zip(bounds[:-1], bounds[1:]):
    a2.axhline(b0, color='k', lw=0.6)
    a2.text(-0.08, (b0 + b1) / 2, genotype.group[b0], ha='right', va='center', fontsize=6.5, transform=a2.get_yaxis_transform())
a2.set_title(f'all {len(genotype)} strains, sorted by group', fontsize=9); a2.set_xlabel('')
plt.tight_layout()
plt.savefig(f'{FIG}/03d_lof_matrix.png', dpi=200)
plt.show()
print(counts.to_string())

# %%
# 3.6 Is LOF enriched in domesticated beverage yeasts? (Fisher's exact test)
from scipy.stats import fisher_exact, chi2_contingency
Gt = genotype[genotype.predicted_POF != '?']
DOMESTIC = ['Beer (Gallone)', 'Ale', 'Lager', 'Sake', 'Shochu/Awamori', 'Spirits']
stat_rows = []
for label, A_, B_ in [('Domesticated beverage vs Wild+Wine+Bioethanol', DOMESTIC, ['Wild', 'Wine', 'Bioethanol']),
                      ('Beer (Gallone)+Ale vs Wild+Wine', ['Beer (Gallone)', 'Ale'], ['Wild', 'Wine'])]:
    a, b = Gt[Gt.group.isin(A_)], Gt[Gt.group.isin(B_)]
    tab = [[(a.predicted_POF == 'POF-').sum(), (a.predicted_POF == 'POF+').sum()],
           [(b.predicted_POF == 'POF-').sum(), (b.predicted_POF == 'POF+').sum()]]
    odds, pval = fisher_exact(tab)
    stat_rows.append(dict(test=label, a_POFneg=tab[0][0], a_POFpos=tab[0][1], b_POFneg=tab[1][0], b_POFpos=tab[1][1],
                          odds_ratio=odds, fisher_p=pval))
ct = pd.crosstab(Gt.group, Gt.predicted_POF)
chi2, chi_p, dof, _ = chi2_contingency(ct)
stat_rows.append(dict(test=f'chi-square, all groups x POF ({dof} df)', fisher_p=np.nan, chi2=chi2, chi2_p=chi_p))
stats_tbl = pd.DataFrame(stat_rows)
stats_tbl.to_csv(f'{RES}/lof_enrichment_tests.csv', index=False)
print(ct.to_string()); print()
print(stats_tbl.to_string(index=False))

# %%
# Computed summary used by the interpretation below, cross-checked against the literature (Mukai 2014; Mertens 2019)
def strain_list(df, cap=12):
    names = sorted(set(df.strain if isinstance(df, pd.DataFrame) else df))
    return (', '.join(names[:cap]) + (f' ... (+{len(names) - cap})' if len(names) > cap else '')) or 'none'

print(f"Records: PAD1 {sum(allele_summary.gene == 'PAD1')}, FDC1 {sum(allele_summary.gene == 'FDC1')} "
      f"({(~allele_summary.callable).sum()} not callable); strains: {len(genotype)}")
for g in VC_GENES:
    A = AS_c[AS_c.gene == g]
    print(f"\n{g}: {int(A.LOF.sum())}/{len(A)} callable records LOF")
    for grp in grp_order:
        sub = A[A.group == grp]
        if len(sub):
            print(f"  {grp:18s} {int(sub.LOF.sum()):3d}/{len(sub):<3d} LOF: {strain_list(sub[sub.LOF == 1])}")
for h, lab in [('c.460C>T', 'FDC1 Q154*'), ('c.501dupA', 'FDC1 c.501dupA (= c.495_496insA left-aligned)'),
               ('c.160A>T', 'FDC1 K54*'), ('c.305G>A', 'PAD1 W102*')]:
    r = recurrent[recurrent.hgvs_c == h]
    print(f"{lab}: " + (f"{r.n_strains.iloc[0]} strains - {r.groups.iloc[0]}" if len(r) else 'none'))
bneg = genotype[(genotype.group == 'Beer (Gallone)') & (genotype.predicted_POF == 'POF-')].strain
q154 = set(VC[(VC.hgvs_c == 'c.460C>T') & (VC.gene == 'FDC1')].strain)
print(f"Q154* among POF- Beer (Gallone) strains: {len(set(bneg) & q154)}/{len(bneg)} "
      f"({100 * len(set(bneg) & q154) / max(1, len(bneg)):.1f} %; Mertens 2019: 73.8 % of POF- Beer 1)")
late = VC[VC.effect == 'nonsense'].merge(AS_c[['gene', 'sacc', 'LOF', 'protein_len']], on=['gene', 'sacc'])
late = late[late.LOF == 0]
print('\nNonsense SNVs NOT called LOF (stop at >= 95 % of reference length):')
print(late.groupby(['gene', 'hgvs_c', 'protein_change', 'protein_len']).strain.apply(strain_list).to_string() if len(late) else '  none')
EXPECT = {'Laboratory': 'POF+', 'Wine': 'POF+', 'Wild': 'POF+', 'Bioethanol': 'POF+', 'Wheat/Hefeweizen': 'POF+',
          'Ale': 'POF-', 'Lager': 'POF-', 'Sake': 'POF-', 'Shochu/Awamori': 'POF-'}
genotype['expected'] = genotype.group.map(EXPECT).fillna('?')
mism = genotype[(genotype.expected != '?') & (genotype.predicted_POF != '?') & (genotype.expected != genotype.predicted_POF)]
print('\nStrains whose predicted POF disagrees with the expectation for their group:')
print(mism[['group', 'strain', 'PAD1', 'FDC1', 'predicted_POF', 'expected']].to_string(index=False) if len(mism) else '  none')
inc = genotype[genotype.predicted_POF == '?']
print(f"\nStrains with no POF call ('?'): {strain_list(inc, 30)}")
amb = variant_calls[variant_calls.type == 'ambiguous']
print('Ambiguous (IUPAC) base calls:', ', '.join(f'{r.gene} {r.strain} {r.hgvs_c}' for r in amb.itertuples()) or 'none')
