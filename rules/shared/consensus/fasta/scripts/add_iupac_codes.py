import sys
import vcf
from Bio import SeqIO
from Bio.Seq import Seq

"""
    Change the nucleotide to the particular IUPAC code at every position in which the SNP's
    alternative allele frequence lies in the [0.35, 0.5) interval.
"""


IUPAC_CODE = {
    'A': {'G':'R', 'T':'W', 'C':'M'},
    'C': {'T':'Y', 'G':'S', 'A':'M'},
    'G': {'A':'R', 'C':'S', 'T':'K'},
    'T': {'C':'Y', 'A':'W', 'G':'K'},
}

ALLELE_FREQ_LOW = 0.35
ALLELE_FREQ_HIGH = 0.65


if __name__ == '__main__':
    IN_FASTA_PATH = sys.argv[1]
    IN_VCF_PATH = sys.argv[2]
    OUT_FASTA_PATH = sys.argv[3]

    sequences = list(SeqIO.parse(open(IN_FASTA_PATH), 'fasta'))
    assert len(sequences) == 1
    fasta = list(sequences[0].seq)

    vcf_reader = vcf.Reader(filename=IN_VCF_PATH)
    for record in vcf_reader:
        if record.is_snp:
            assert len(record.alleles) == 2
            sample = record.samples[0]
            if ALLELE_FREQ_LOW <= sample['AD'][1]/(sample['AD'][0]+sample['AD'][1]) < ALLELE_FREQ_HIGH:
                fasta[record.POS - 1] = IUPAC_CODE[record.alleles[0]][str(record.alleles[1])]

    sequences[0].seq = Seq(''.join(fasta))
    with open(OUT_FASTA_PATH, 'w') as out_file:
        SeqIO.write(sequences, out_file, 'fasta')
