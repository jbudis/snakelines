import sys
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

input_r1, input_r2, input_rm, params_n_ambiguous_bases = sys.argv[1:5]
output_rc = sys.argv[5]

with gzip.open(output_rc, 'wt') as out:

    with gzip.open(input_rm, 'rt') as rm:
        for read in SeqIO.parse(rm, 'fastq'):
            SeqIO.write(read, out, 'fastq')

    with gzip.open(input_r1, 'rt') as r1, gzip.open(input_r2, 'rt') as r2:
        for read1, read2 in zip(SeqIO.parse(r1, 'fastq'), SeqIO.parse(r2, 'fastq')):

            joined_seq = '{}{}{}'.format(read1.seq, 'N'*params_n_ambiguous_bases, read2.seq)
            joined_qual = {'phred_quality': read1.letter_annotations['phred_quality'] + \
                                                                                                [1]*params_n_ambiguous_bases + \
                                                                                            read2.letter_annotations['phred_quality']}
            joined_description = '{} {}|{}'.format(read1.id, read1.description.split(' ')[1], read2.description.split(' ')[1])
            joined = SeqRecord(id=read1.id, description=joined_description,
                                seq=Seq(joined_seq, IUPAC.IUPACAmbiguousDNA()), letter_annotations=joined_qual)

            SeqIO.write(joined, out, 'fastq')
