"""Preprocess ImmunoVerse per-cancer CSVs into compact JSON for the portal."""
import pandas as pd
import json, ast, re, os
from pathlib import Path

SRC = Path('/sessions/clever-zen-mayer/mnt/working files/ImmunoVerse/table')
OUT = Path('/sessions/clever-zen-mayer/work/portal/data')
OUT.mkdir(parents=True, exist_ok=True)

CANCER_META = {
    'AML':  ('Acute Myeloid Leukemia',           'Heme'),
    'BLCA': ('Bladder Urothelial Carcinoma',      'GU'),
    'BRCA': ('Breast Invasive Carcinoma',         'Breast'),
    'CESC': ('Cervical Squamous Cell Carcinoma',  'Gyn'),
    'CHOL': ('Cholangiocarcinoma',                'GI'),
    'COAD': ('Colon Adenocarcinoma',              'GI'),
    'DLBC': ('Diffuse Large B-cell Lymphoma',     'Heme'),
    'ESCA': ('Esophageal Carcinoma',              'GI'),
    'GBM':  ('Glioblastoma Multiforme',           'CNS'),
    'HNSC': ('Head & Neck Squamous Cell Carc.',   'H&N'),
    'KIRC': ('Kidney Renal Clear Cell Carc.',     'GU'),
    'LIHC': ('Liver Hepatocellular Carcinoma',    'GI'),
    'LUAD': ('Lung Adenocarcinoma',               'Lung'),
    'LUSC': ('Lung Squamous Cell Carcinoma',      'Lung'),
    'MESO': ('Mesothelioma',                      'Thoracic'),
    'NBL':  ('Neuroblastoma',                     'Pediatric'),
    'OV':   ('Ovarian Serous Cystadenocarc.',     'Gyn'),
    'PAAD': ('Pancreatic Adenocarcinoma',         'GI'),
    'RT':   ('Rhabdoid Tumor',                    'Pediatric'),
    'SKCM': ('Skin Cutaneous Melanoma',           'Skin'),
    'STAD': ('Stomach Adenocarcinoma',            'GI'),
}

CLASS_LABELS = {
    'self_gene':             'Canonical self-antigen',
    'splicing':              'Alternative splicing',
    'variant':               'Variant / neoantigen',
    'nuORF':                 'Cryptic / non-canonical ORF',
    'ERV':                   'Endogenous retrovirus',
    'TE_chimeric_transcript':'TE chimeric transcript',
    'intron_retention':      'Intron retention',
    'fusion':                'Gene fusion',
    'pathogen':              'Pathogen-derived',
}


def safe_lit(v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    if not isinstance(v, str):
        return v
    s = v.strip()
    if not s or s.lower() == 'nan':
        return None
    try:
        return ast.literal_eval(s)
    except Exception:
        return s


def clean_gene(row):
    """Return a clean, human-readable gene symbol string.

    Priority: gene_symbol column (dedup, join with '/') → parsed from source →
    first ENSG from ensgs → empty string.
    """
    g = safe_lit(row.get('gene_symbol'))
    names = []
    if isinstance(g, (list, tuple)):
        for x in g:
            if x is None: continue
            x = str(x).strip()
            if x and x.lower() != 'nan' and x not in names:
                names.append(x)
    elif isinstance(g, str):
        for x in re.split(r'[;,/]+', g):
            x = x.strip()
            if x and x.lower() != 'nan' and x not in names:
                names.append(x)
    if names:
        return '/'.join(names[:3])

    # Parse from source
    src = row.get('source')
    if isinstance(src, str) and src:
        # Look for symbol-like tokens at end of pipe-delimited fields
        # e.g. chr12:...|ENST...|PMEL,PMEL   or   ENSG...|ENST...|HLA-C;ENSG...|...|HLA-B
        found = []
        for piece in re.split(r'[;,]', src):
            m = re.findall(r'\|([A-Za-z0-9\-_.]+)\s*$', piece)
            for hit in m:
                if hit and not hit.startswith('ENS') and hit not in found and hit not in ('nuORF','ERV','TE'):
                    found.append(hit)
        if found:
            return '/'.join(found[:3])

    return ''


def clean_hlas(row):
    """Return de-duplicated list of HLA allele strings (HLA-A*02:01 → A*02:01)."""
    d = safe_lit(row.get('presented_by_each_sample_hla')) or {}
    alleles = set()
    if isinstance(d, dict):
        for _, hits in d.items():
            if not isinstance(hits, (list, tuple)):
                continue
            for h in hits:
                if isinstance(h, (list, tuple)) and h:
                    allele = str(h[0])
                else:
                    allele = str(h)
                allele = allele.replace('HLA-','').strip()
                if allele and allele.lower() != 'none':
                    alleles.add(allele)
    # Sort canonical: A, B, C then numeric
    def k(a):
        return (a[:1], a)
    return sorted(alleles, key=k)


def num_or_none(v):
    if v is None:
        return None
    if isinstance(v, float) and pd.isna(v):
        return None
    try:
        f = float(v)
        if pd.isna(f): return None
        return f
    except Exception:
        return None


def first_num_list(v):
    """Parse a list-like and return first numeric value."""
    x = safe_lit(v)
    if isinstance(x, (list, tuple)) and x:
        for y in x:
            n = num_or_none(y)
            if n is not None:
                return n
    return num_or_none(x)


def min_num_list(v):
    x = safe_lit(v)
    if isinstance(x, (list, tuple)) and x:
        nums = [num_or_none(y) for y in x]
        nums = [n for n in nums if n is not None]
        if nums: return min(nums)
    return num_or_none(x)


def homogeneity(row):
    """A simple score from sc_pert list if present (mean)."""
    x = safe_lit(row.get('sc_pert'))
    if isinstance(x, (list, tuple)) and x:
        nums = [num_or_none(y) for y in x if num_or_none(y) is not None]
        if nums: return round(sum(nums)/len(nums), 3)
    return None


def process_cancer(code):
    fp = SRC / f'{code}_final_enhanced_all.csv'
    if not fp.exists():
        print(f'[skip] {code}: no file')
        return None
    df = pd.read_csv(fp, low_memory=False)
    df = df.fillna(value={'relative_abundance': 0, 'highest_score': 0, 'n_psm': 0})

    rows = []
    class_count = {}
    gene_set = set()

    for _, r in df.iterrows():
        pep   = str(r.get('pep', '')).strip()
        if not pep:
            continue
        cls   = str(r.get('typ', '')).strip() or 'self_gene'
        psm   = int(num_or_none(r.get('n_psm')) or 0)
        score = num_or_none(r.get('highest_score'))
        abund = num_or_none(r.get('relative_abundance')) or 0
        qval  = num_or_none(r.get('best_pep'))
        tumor = first_num_list(r.get('median_tumor'))
        norm  = first_num_list(r.get('max_median_gtex'))
        dep   = first_num_list(r.get('depmap_median'))
        homo  = homogeneity(r)
        uniq  = 1 if str(r.get('unique','')).strip().lower() in ('true','1') else 0
        gene  = clean_gene(r)
        hlas  = clean_hlas(r)

        rows.append([
            pep,          # 0
            gene,         # 1
            cls,          # 2
            psm,          # 3
            round(score, 1) if score is not None else None,  # 4
            round(abund, 4) if abund is not None else None,  # 5
            qval,         # 6
            round(tumor, 2) if tumor is not None else None,  # 7
            round(norm, 2) if norm is not None else None,    # 8
            round(dep, 3) if dep is not None else None,      # 9
            homo,         # 10
            uniq,         # 11
            hlas,         # 12
        ])
        class_count[cls] = class_count.get(cls, 0) + 1
        if gene:
            for g in gene.split('/'):
                gene_set.add(g)

    # Sort by relative abundance desc, then spectral score desc
    rows.sort(key=lambda r: (-(r[5] or 0), -(r[4] or 0)))

    name, group = CANCER_META.get(code, (code, ''))
    out = {
        'code': code,
        'name': name,
        'group': group,
        'count': len(rows),
        'classes': class_count,
        'genes': len(gene_set),
        'cols': ['pep','gene','class','psm','score','abund','qval','tumorExp','normalExp','depmap','homo','unique','hlas'],
        'rows': rows,
    }
    with open(OUT / f'{code}.json', 'w') as f:
        json.dump(out, f, separators=(',', ':'), default=str)
    print(f'{code}: {len(rows):>6,} rows, {len(gene_set):>5} genes')
    return {
        'code': code, 'name': name, 'group': group,
        'count': len(rows), 'classes': class_count, 'genes': len(gene_set),
    }


def main():
    cancers = []
    for code in CANCER_META.keys():
        meta = process_cancer(code)
        if meta:
            cancers.append(meta)

    # Summary
    total = sum(c['count'] for c in cancers)
    class_tot = {}
    for c in cancers:
        for k, v in c['classes'].items():
            class_tot[k] = class_tot.get(k, 0) + v

    summary = {
        'cancers': cancers,
        'stats': {
            'cancers': len(cancers),
            'immunopeptidomes': 1823,
            'rnaSeq': 7188,
            'normalControls': 17384,
            'scRNA': 594,
            'classes': len(class_tot),
        },
        'classes': class_tot,
        'totals': {'targets': total, 'genes': sum(c['genes'] for c in cancers)},
        'labels': CLASS_LABELS,
    }
    with open(OUT / '_summary.json', 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f'\nTOTAL: {total:,} rows across {len(cancers)} cancers')
    print('Classes:', class_tot)


if __name__ == '__main__':
    main()
