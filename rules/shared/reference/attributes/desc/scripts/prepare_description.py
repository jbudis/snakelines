import sys
from Bio import SeqIO

input_reference = sys.argv[2]
output_dict = sys.argv[3]


with open(output_dict, 'w') as out:
    for seq in SeqIO.parse(input_reference, 'fasta'):
        desc = seq.description
        if desc.startswith(seq.id):
            desc = desc[len(seq.id) + 1:]
        out.write('{}\t{}\n'.format(seq.id, desc))

