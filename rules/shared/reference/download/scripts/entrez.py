import sys
from Bio import Entrez
from Bio import SeqIO

params_genbank_ids = sys.argv[1]
params_email = sys.argv[2]
output_reference = sys.argv[3]
output_taxonomy = sys.argv[4]

with open(output_reference, 'w') as refs, open(output_taxonomy, 'w') as taxes:
    for genbank_id in params_genbank_ids:
        handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        SeqIO.write(record, refs, 'fasta')

        taxonomy = ';'.join(record.annotations['taxonomy'])
        if 'organism' in record.annotations:
            organism = record.annotations['organism']
            taxonomy = '{taxonomy};{organism}'.format(taxonomy=taxonomy, organism=organism)

        taxes.write('{}\t{}\n'.format(genbank_id, taxonomy))
