#!/usr/bin/python3

import sys
import pickle

input_file = sys.argv[1]    # spades/KMB1/joined/blast/contigs.fasta.kingdoms
output_file = sys.argv[2]   # spades/KMB1/joined/blast/contigs.fasta.annotated

taxonomies = pickle.load(open('/data/genome/metagenome/blast/nt/17-01-17/taxonomy/code_taxonomy.pickle', 'rb'))
node = ""
record = ""
output = open(output_file, 'w')

with open(input_file) as f:
    while True:

        line = f.readline()

        if line == "":
            break

        line = line.split("\t")

        if node != line[0]:
            output.write(record + "\n")
            record = line[0] + "\t"
            node = line[0]

        try:
            new_taxonomy = taxonomies[int(line[4][:-1])]
            record += ",".join(new_taxonomy) + ";"
        except KeyError:
            print("KeyError" + line[4])

        # print(record)

output.close()
