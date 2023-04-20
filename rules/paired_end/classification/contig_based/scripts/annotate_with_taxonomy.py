import sys
import pickle
import pandas as pd
from glob import glob
from Bio import SeqIO

input_blast = sys.argv[1]
input_taxes = sys.argv[2]
output_blast = sys.argv[3]

taxes = pickle.load(open(input_taxes, 'rb'))

def get_tax(taxid):
    if not taxid or taxid == 'N/A':
        return 'Unknown'
    taxid = int(taxid)
    return ';'.join(taxes.get(taxid, ''))

with open(input_blast) as in_blast, open(output_blast, 'w') as out_blast:
    for i, line in enumerate(in_blast):
        items = line.strip().split("\t")
        if i == 0:
            tax_col = items.index('staxid')
            to_insert = 'taxonomy'
        else:
            tax = get_tax(items[tax_col])
            to_insert = tax

        items = items[:tax_col] + [to_insert] + items[tax_col:]
        out_blast.write('\t'.join(items))
        out_blast.write('\n')
