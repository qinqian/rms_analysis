#!/usr/bin/env python

import sys

fa = sys.argv[1]
gtf = open(fa.replace('.fa', '.gtf'), 'w')

with open(fa) as fa:
    for line in fa:
        if line.startswith('>'):
            name = line.strip().replace('>', '')
            gene = '_'.join(name.split('.')[0:2])
        else: 
            seq_len = len(line.strip())
            gtf.write(f'{name}\texternal\tgene\t1\t{seq_len}\t.\t-\t.\tgene_id "{gene}"; gene_version "1"; gene_name "{gene}"; gene_source "external"; gene_biotype "external";\n')
            gtf.write(f'{name}\texternal\ttranscript\t1\t{seq_len}\t.\t-\t.\tgene_id "{gene}"; gene_version "1"; transcript_id "{name}"; transcript_version "1"; gene_name "{gene}"; gene_source "external"; gene_biotype "external"; transcript_name "{name}"; transcript_source "external"; transcript_biotype "external";\n')
            gtf.write(f'{name}\texternal\texon\t1\t{seq_len}\t.\t-\t.\tgene_id "{gene}"; gene_version "1"; transcript_id "{name}"; transcript_version "1"; exon_number "1"; gene_name "{gene}"; gene_source "external"; gene_biotype "external"; transcript_name "{name}"; transcript_source "external"; transcript_biotype "external"; exon_id "{gene}_exon1"; exon_version "1";\n')
