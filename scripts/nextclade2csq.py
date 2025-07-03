#!/usr/bin/env python3
"""
Convert a Nextclade-generated GFF3 into a minimal GFF3 for bcftools csq.

Requirements:
  - gene lines:      ID=gene:<gene_id>;biotype=<biotype>;Name=<gene_name> (optional)
  - transcript line: ID=transcript:<tx_id>;Parent=gene:<gene_id>;biotype=<biotype>
  - CDS lines:       Parent=transcript:<tx_id>

Usage:
  chmod +x nextclade2csq.py
  ./nextclade2csq.py input.nextclade.gff3 > output.csq.gff3
"""
import sys, urllib.parse
from collections import defaultdict

if len(sys.argv) != 2:
    sys.exit("Usage: nextclade2csq.py <nextclade.gff3>")

infile = sys.argv[1]

# Collect gene intervals and CDS records
genes = {}  # gene_id -> dict {seqid, start, end, biotype, gene_name}
cds_records = []  # list of dicts for each CDS

for line in open(infile):
    if line.startswith('##'):
        # preserve only version and sequence-region
        if line.startswith('##gff-version') or line.startswith('##sequence-region'):
            print(line.rstrip())
        continue
    cols = line.rstrip().split('\t')
    if len(cols) < 9:
        continue
    seqid, _, feat_type, start, end, _, strand, phase, attrs = cols
    # We're only interested in CDS features
    if feat_type != 'CDS':
        continue
    A = {k: urllib.parse.unquote(v) for k,v in (x.split('=',1) for x in attrs.split(';') if '=' in x)}
    gene_id   = A.get('gene_name') or A.get('gene') or A.get('Name')
    biotype   = 'protein_coding'
    gene_name = A.get('gene_name') or A.get('gene')
    tx_id     = gene_id + '_t1'
    prot_id   = A.get('protein_id') or tx_id + '_cds'
    # record gene interval
    s,e = int(start), int(end)
    if gene_id not in genes:
        genes[gene_id] = dict(seqid=seqid, start=s, end=e, biotype=biotype, gene_name=gene_name)
    else:
        g = genes[gene_id]
        g['start'] = min(g['start'], s)
        g['end']   = max(g['end'],   e)
    # store CDS for later
    cds_records.append(dict(seqid=seqid, start=s, end=e,
                            strand=strand if strand in ('+','-') else '+',
                            phase=phase, tx_id=tx_id, prot_id=prot_id))

# Output gene and transcript lines
for gene_id, info in genes.items():
    seqid = info['seqid']
    s,e   = info['start'], info['end']
    biotype = info['biotype']
    gene_name = info.get('gene_name')
    # gene line
    attrs = [f"ID=gene:{gene_id}", f"biotype={biotype}"]
    if gene_name:
        attrs.append(f"Name={gene_name}")
    print('\t'.join([seqid,'.', 'gene', str(s), str(e), '.', '+', '.', ';'.join(attrs)]))
    # transcript line
    tx_id = gene_id + '_t1'
    attrs = [f"ID=transcript:{tx_id}", f"Parent=gene:{gene_id}", f"biotype={biotype}"]
    print('\t'.join([seqid,'.','transcript', str(s), str(e), '.', '+', '.', ';'.join(attrs)]))

# Output CDS lines
for rec in cds_records:
    attrs = [f"Parent=transcript:{rec['tx_id']}"
             ]
    print('\t'.join([rec['seqid'], '.', 'CDS',
                       str(rec['start']), str(rec['end']),
                       '.', rec['strand'], rec['phase'],
                       ';'.join(attrs)]))
