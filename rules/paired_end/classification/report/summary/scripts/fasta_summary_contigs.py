import os
import sys
import shutil
import zlib
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import SeqUtils


input_fasta = sys.argv[1]
params_fasta = sys.argv[2]
output_fasta_dir = sys.argv[3]
output_seqinfo = sys.argv[4]

os.makedirs(output_fasta_dir,exist_ok=True)
shutil.copyfile(input_fasta, params_fasta)

seqids, info_list = [], []
for i, seq in enumerate(SeqIO.parse(input_fasta, 'fasta')):

    # Store each sequence in the fasta to the separate file
    seq_path = os.path.join(output_fasta_dir, seq.id)
    with open(seq_path, 'w') as out:
        SeqIO.write(seq, out, 'fasta')

    seqids.append(seq.id)
    contig_info = {
        'Length':     len(seq),
        'GC content': str(np.round(SeqUtils.GC(seq.seq), 1)),
        'Compress ratio':  len(zlib.compress(bytes(str(seq.seq), 'utf-8'))) / len(seq.seq)
    }
    info_list.append(contig_info)

infos = pd.DataFrame(info_list, index=seqids)
infos.to_csv(output_seqinfo, sep='\t')
