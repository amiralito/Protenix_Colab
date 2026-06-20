#!/usr/bin/env python3
"""Generate Protenix batch input JSON(s).

A Protenix input file is a top-level LIST of job dicts; `protenix pred`/`protenix msa`
iterate over every entry. This builds that list from FASTA files for common screen designs.

Examples
--------
# one job per sequence (monomers)
python make_protenix_batch.py monomer --fasta seqs.fasta -o batch_inputs

# homo-oligomer: each sequence as an N-mer (e.g. a hexameric resistosome cap)
python make_protenix_batch.py homomer --fasta nlrs.fasta --copies 6 -o batch_inputs

# every effector x NLR pair (cartesian product of two FASTAs)
python make_protenix_batch.py all_pairs --fasta_a effectors.fasta --fasta_b nlrs.fasta -o batch_inputs

# explicit pairs from a TSV/CSV (columns: idA,idB) with sequences pulled from the FASTAs
python make_protenix_batch.py pairs --pairs pairs.tsv --fasta_a effectors.fasta --fasta_b nlrs.fasta -o batch_inputs

Optional add-ons applied to every job:  --ligand CCD_ATP   --ion MG
Split a big screen into runnable chunks:  --chunk 50   (writes batch_001.json, batch_002.json, ...)
"""
import argparse, csv, json, os, re, sys
from itertools import product


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name = line[1:].split()[0]          # first whitespace-delimited token
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    if not seqs:
        sys.exit(f'No sequences parsed from {path}')
    return seqs


_SAN = re.compile(r'[^A-Za-z0-9_.-]+')
def sanitize(s):
    return _SAN.sub('_', s).strip('_') or 'job'

def clean_seq(s):
    return re.sub(r'\s+', '', s).upper()

def protein(seq, count=1):
    return {'proteinChain': {'sequence': clean_seq(seq), 'count': int(count)}}

def ligand(spec, count=1):
    return {'ligand': {'ligand': spec, 'count': int(count)}}

def ion(code, count=1):
    return {'ion': {'ion': code, 'count': int(count)}}

def add_extras(ents, a):
    if a.ligand:
        ents.append(ligand(a.ligand, a.ligand_copies))
    if a.ion:
        ents.append(ion(a.ion, a.ion_copies))
    return ents


def build_jobs(a):
    jobs = []
    if a.mode == 'monomer':
        for nm, sq in read_fasta(a.fasta).items():
            jobs.append({'name': sanitize(nm), 'sequences': add_extras([protein(sq, 1)], a)})

    elif a.mode == 'homomer':
        for nm, sq in read_fasta(a.fasta).items():
            jobs.append({'name': sanitize(f'{nm}_x{a.copies}'),
                         'sequences': add_extras([protein(sq, a.copies)], a)})

    elif a.mode == 'all_pairs':
        A = read_fasta(a.fasta_a)
        if a.fasta_b:
            pairs = list(product(A.items(), read_fasta(a.fasta_b).items()))
        else:
            items = list(A.items()); pairs = []
            for i in range(len(items)):
                lo = i if a.include_self else i + 1
                for j in range(lo, len(items)):
                    pairs.append((items[i], items[j]))
        for (na, sa), (nb, sb) in pairs:
            ents = add_extras([protein(sa, a.copies_a), protein(sb, a.copies_b)], a)
            jobs.append({'name': sanitize(f'{na}__{nb}'), 'sequences': ents})

    elif a.mode == 'pairs':
        A = read_fasta(a.fasta_a)
        B = read_fasta(a.fasta_b) if a.fasta_b else A
        with open(a.pairs) as f:
            sniff = f.read(4096); f.seek(0)
            delim = '\t' if sniff.count('\t') >= sniff.count(',') else ','
            rows = [r for r in csv.reader(f, delimiter=delim) if r and not r[0].lstrip().startswith('#')]
        if rows and rows[0][0].strip().lower() in ('ida', 'id_a', 'a', 'effector'):
            rows = rows[1:]
        for r in rows:
            ida, idb = r[0].strip(), r[1].strip()
            if ida not in A:
                sys.exit(f'pair id "{ida}" not found in fasta_a')
            if idb not in B:
                sys.exit(f'pair id "{idb}" not found in fasta_b')
            ents = add_extras([protein(A[ida], a.copies_a), protein(B[idb], a.copies_b)], a)
            jobs.append({'name': sanitize(f'{ida}__{idb}'), 'sequences': ents})

    # make job names unique
    seen = {}
    for jb in jobs:
        base = jb['name']; n = seen.get(base, 0)
        if n:
            jb['name'] = f'{base}_{n + 1}'
        seen[base] = n + 1
    if not jobs:
        sys.exit('No jobs generated - check your inputs.')
    return jobs


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('mode', choices=['monomer', 'homomer', 'all_pairs', 'pairs'])
    p.add_argument('--fasta'); p.add_argument('--fasta_a'); p.add_argument('--fasta_b'); p.add_argument('--pairs')
    p.add_argument('--copies', type=int, default=2, help='homomer copy number')
    p.add_argument('--copies_a', type=int, default=1)
    p.add_argument('--copies_b', type=int, default=1)
    p.add_argument('--include_self', action='store_true', help='all_pairs within one FASTA: also pair i with itself')
    p.add_argument('--ligand', help='CCD_XXX, a SMILES string, or FILE_/path.sdf')
    p.add_argument('--ligand_copies', type=int, default=1)
    p.add_argument('--ion', help='bare CCD code, e.g. MG, ZN, NA')
    p.add_argument('--ion_copies', type=int, default=1)
    p.add_argument('-o', '--out_dir', default='batch_inputs')
    p.add_argument('--name', default='batch', help='base name for the output JSON file(s)')
    p.add_argument('--chunk', type=int, default=0, help='split into files of N jobs (0 = single file)')
    a = p.parse_args()

    # minimal arg sanity per mode
    if a.mode in ('monomer', 'homomer') and not a.fasta:
        p.error(f'{a.mode} needs --fasta')
    if a.mode == 'all_pairs' and not a.fasta_a:
        p.error('all_pairs needs --fasta_a (and optionally --fasta_b)')
    if a.mode == 'pairs' and not (a.pairs and a.fasta_a):
        p.error('pairs needs --pairs and --fasta_a (and --fasta_b if cross-set)')

    jobs = build_jobs(a)
    os.makedirs(a.out_dir, exist_ok=True)
    chunks = [jobs] if a.chunk <= 0 else [jobs[i:i + a.chunk] for i in range(0, len(jobs), a.chunk)]

    manifest = []
    for ci, ch in enumerate(chunks):
        suffix = '' if len(chunks) == 1 else f'_{ci + 1:03d}'
        fp = os.path.join(a.out_dir, f'{a.name}{suffix}.json')
        with open(fp, 'w') as f:
            json.dump(ch, f, indent=2)
        manifest += [(os.path.basename(fp), jb['name'], len(jb['sequences'])) for jb in ch]
        print(f'wrote {fp}  ({len(ch)} jobs)')

    with open(os.path.join(a.out_dir, f'{a.name}_manifest.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['json_file', 'job_name', 'n_entities']); w.writerows(manifest)

    print(f'\nTotal {len(jobs)} jobs across {len(chunks)} file(s) in {a.out_dir}/  (+ manifest).')
    print('Next, per file:')
    print('  protenix msa  --input <file>.json --out_dir ./msa --msa_server_mode colabfold')
    print('  protenix pred -i <file>-update-msa.json -o ./outputs -n protenix-v2 --use_msa true')


if __name__ == '__main__':
    main()
