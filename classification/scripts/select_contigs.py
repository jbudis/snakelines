#!/usr/bin/python3

import sys

contigs_fasta = sys.argv[1]
viral_contigs = sys.argv[2]
output_file = sys.argv[3]

viral_list = open(viral_contigs).readlines()
for i,record in enumerate(viral_list):
    viral_list[i] = record.rstrip()

output = open(output_file, mode='w')
viral = False

#print(viral_list)

with open(contigs_fasta) as f:
    while True:
        line = f.readline()
        #print(line)

        if line == "":
            break
        
        if line[0] == ">":
            node = line[1:].split(" ")[0].rstrip()
            #print(node)
            if node in viral_list:
                viral = True
            else:
                viral = False
            
        if viral:
            output.write(line)
