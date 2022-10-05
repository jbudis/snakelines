import sys
import pandas as pd
import csv

blast = pd.read_csv(sys.argv[1], sep='\t', index_col=None, encoding='utf-8', quoting=csv.QUOTE_NONE)

min_query_coverage = float(sys.argv[2])
max_target_seqs = int(sys.argv[3])
reference = sys.argv[4]
seqids, info_list = [], []

for seqid, aligns in blast.groupby('qseqid'):
    aligns['qcov'] = ((aligns['qend'] - aligns['qstart']) / aligns['qlen']).abs()
    aligns = aligns[aligns['qcov'] >= min_query_coverage]
    aligns = aligns.iloc[:max_target_seqs]

    if len(aligns) == 0:
        continue

    def create_link(align):
        return '%3.2f%% - <a href="https://www.ncbi.nlm.nih.gov/nuccore/%s", title="%s">%s</a>' % \
                            (align['qcov'] * 100, align['sacc'], align['taxonomy'], align['stitle'])

    seqids.append(seqid)
    info_list.append({
        'Homologue title (%s)' % reference: ';'.join(aligns['stitle']),
        'Homologue accession (%s)' % reference: ';'.join(aligns['stitle']),
        'Homologue link (%s)' % reference: '<br />'.join(aligns.apply(create_link, axis=1))
    })

infos = pd.DataFrame(info_list, index=seqids)
infos.to_csv(sys.argv[5], sep='\t')
