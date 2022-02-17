import sys
from Bio import SeqIO

input_reference = sys.argv[1]
output_dict = sys.argv[2]


with open(output_dict, 'w') as out:
    for seq in SeqIO.parse(input_reference, 'fasta'):
        desc = seq.description
        if desc.startswith(seq.id):
            desc = desc[len(seq.id) + 1:]
        out.write('{}\t{}\n'.format(seq.id, desc))

