import pandas as pd
import numpy as np
import sys
import os

def format_list(array):
    return ['{:.1%}'.format(x) for x in array]

table_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'covid_table-21K_21L.tsv')
df = pd.read_csv(table_path, sep='\t')
summary = []
pokrytie = 100
min_variants = 0

input_file = sys.argv[1]
output_file = sys.argv[2]

genomic_positions = []
mutations = []
weights_delta = []
weights_21K = []
weights_21L = []
delta_AF = []
omicron21K_AF = []
omicron21L_AF = []
with open(input_file, 'r') as infile:
    for line in infile:
        if line.startswith('#'):
            continue
        line = line.strip().split()
        info = line[7]
        info = info.split(';')
        position = int(line[1])
        df_annotation = df.loc[df.genomic_position == position]
        variant_type = ', '.join(map(str, df_annotation['variant_type'].unique()))
        
        if variant_type != '':
            save = False
            alt = line[4]
            index_DP = [idx for idx, s in enumerate(info) if 'DP=' in s][0]
            DP = int(info[index_DP].split('=')[1])

            if alt == '.' and DP > pokrytie:
                save = True
                if 'delta' in variant_type:
                    delta_AF.append(0)
                    omicron21K_AF.append('')
                    omicron21L_AF.append('')
                    weights_delta.append(DP)
                elif 'omicron (21K)' in variant_type:
                    omicron21K_AF.append(0)
                    delta_AF.append('')
                    omicron21L_AF.append('')
                    weights_21K.append(DP)
                elif 'omicron (21L)' in variant_type:
                    omicron21L_AF.append(0)
                    delta_AF.append('')
                    omicron21K_AF.append('')
                    weights_21L.append(DP)
                else:
                    save = False

            elif df_annotation['ref:alt'].values[0].split(':')[1] in alt and DP > pokrytie:
                save = True
                index_DP4 = [idx for idx, s in enumerate(info) if 'DP4=' in s][0]
                DP4 = info[index_DP4].split('=')[1].split(',')
                ref_sum = int(DP4[0]) + int(DP4[1])
                alt_sum = int(DP4[2]) + int(DP4[3])
                if 'delta' in variant_type:
                    delta_AF.append(round(alt_sum/DP, 4))
                    omicron21K_AF.append('')
                    omicron21L_AF.append('')
                    weights_delta.append(DP)
                elif 'omicron (21K)' in variant_type:
                    omicron21K_AF.append(round(alt_sum/DP, 4))
                    delta_AF.append('')
                    omicron21L_AF.append('')
                    weights_21K.append(DP)
                elif 'omicron (21L)' in variant_type:
                    omicron21L_AF.append(round(alt_sum/DP, 4))
                    delta_AF.append('')
                    omicron21K_AF.append('')
                    weights_21L.append(DP)
                else:
                    save = False

            if save:
                genomic_positions.append(position)
                mutations.append(', '.join(map(str, df_annotation['nonsynonymous_mutation'].tolist())))

    result = pd.DataFrame({'genomic_position': genomic_positions, 'nonsynonymous_mutation': mutations, 'delta_AF': delta_AF, 'omicron_AF (21K)': omicron21K_AF, 'omicron_AF (21L)': omicron21L_AF})

    delta_AF = list(filter(lambda x: x != "", delta_AF))
    omicron21K_AF = list(filter(lambda x: x != "", omicron21K_AF))
    omicron21L_AF = list(filter(lambda x: x != "", omicron21L_AF))
    
    avg = ['','','']  # weighted average [delta, 21K, 21L]

    if len(delta_AF) > min_variants and len(omicron21K_AF) > min_variants and len(omicron21L_AF) > min_variants:
        for (i, AF, w) in zip(range(3), [delta_AF, omicron21K_AF, omicron21L_AF], [weights_delta, weights_21K, weights_21L]):
            avg[i] = np.average(AF, weights=w)
        s = np.sum(avg)
        if s > 1:
            avg = avg/s
        avg = format_list(avg)
    
    result.loc[len(result.index)] = [input_file[:-4], 'delta weighted average: ' + avg[0], 'omicron (21K) weighted average: ' + avg[1], 'omicron (21L) weighted average: ' + avg[2], '']
    result.to_csv(input_file[:-4] + '.tsv', sep='\t', index=False)
    summary.append([input_file[:-4]] + avg)

df_summary = pd.DataFrame(summary, columns=['sample_id', 'delta_weighted_average', 'omicron_weighted_average (21K)', 'omicron_weighted_average (21L)'])
df_summary.to_csv(output_file, sep='\t', index=False)
