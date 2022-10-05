import pandas as pd
import sys

input_tsv = sys.argv[1]
metadata = sys.argv[2]
class_col = sys.argv[3]
if class_col == 'None':
    raise ValueError("Class column must be set, please specify class_col in config under differential_abundance")
subclass_col = None if sys.argv[4] == 'None' else sys.argv[4]
if subclass_col is not None:
    subject_col = None if sys.argv[5] == 'None' else sys.argv[5]
else:
    subject_col = None if sys.argv[4] == 'None' else sys.argv[4]
output_file = sys.argv[-1]

def process_row(x : str):
    """Convert taxa string into lefse-acceptable format"""
    for c in ['-', '[', ']', '(', ')', ' ']:
        x = x.replace(c, '_')

    # abbreviations for kingdom, phylum, class etc.
    abbr = ['k','p','c','o','f','g','s']

    # split input taxa on ; to get individual taxa ranks
    ranks = x.split(';')
    base = ranks[0]
    if '__' not in base:
        result = [abbr[0]+'__'+base]
    else:
        result = [base]

    # i denotes the rank index, i.e. kingdom=1, phylum=2 etc.
    for i, rank in enumerate(ranks[1:], start=2):
        if rank == '__':
            result.append(f'{base}_x__L{i}')
        elif '__' not in rank and rank != '':
            result.append(f'{abbr[i-1]}__{rank}')
            base = abbr[i-1]+'__'+rank
        elif rank == '':
            result.append(f'{base}_{abbr[i-1]}__L{i}')
        elif rank.split('__')[1] == '':
            result.append(f'{base}_{rank}L{i}')
        else:
            result.append(rank)
            base = rank

    # return string of ranks separated with |
    return '|'.join(result)

# load input taxa csv file
df = pd.read_csv(input_tsv, sep='\t')
zero = df.columns[0]
df[zero] = df[zero].apply(lambda x: process_row(x))
df.rename(columns = {zero:'sample-id'},inplace=True)
header = list(df.columns)[1:]

mf = pd.read_csv(metadata, sep='\t')
mf = mf.replace(' ', '_', regex=True)

mf = mf.loc[mf['sample-id'].isin(header)]
cols = mf.columns.to_list()

if subclass_col is None and subject_col is None:
    mf = mf[['sample-id',class_col]]
elif subclass_col is not None and subject_col is None:
    mf = mf[['sample-id',class_col, subclass_col]]
elif subclass_col is None and subject_col is not None:
    mf = mf[['sample-id',class_col, subject_col]]
else:
    mf = mf[['sample-id',class_col, subclass_col, subject_col]]

class_cols = [class_col]
subclass_cols = [subclass_col]
subject_cols = [subject_col]
for i in header:
    class_cols.append(mf.loc[mf['sample-id']==i][class_col].values[0])
    if subclass_col:
        subclass_cols.append(mf.loc[mf['sample-id']==i][subclass_col].values[0])
    if subject_col:
        subject_cols.append(mf.loc[mf['sample-id']==i][subject_col].values[0])

dct = {}
dct[class_col] = class_cols
num_cols = 1
if len(subclass_cols)>1:
    num_cols += 1
    dct[subclass_col] = subclass_cols
if len(subject_cols)>1:
    num_cols += 1
    dct[subject_col] = subject_cols
dct["empty"] = []

# store converted taxas into lefse-acceptable format of dataframe
new_df = pd.DataFrame.from_dict(dct,orient='index',columns = list(df.columns))
final_df = new_df.append(df, ignore_index=True)
final_df.iloc[num_cols] = list(df.columns)
final_df.to_csv(output_file, header=False, sep='\t',index = False)
