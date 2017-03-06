#!/usr/bin/python3

import sys

input_file = sys.argv[1]
searched_expression = sys.argv[2]
output = open(sys.argv[3], 'w')

with open(input_file) as f:
    while True:
        
        line = f.readline()

        if line == "\n":
            continue
        if line == "":
            break       

        all_taxonomies = line.split('\t')[1].split(';')
        is_searched = False
        for taxonomy in all_taxonomies:
            taxonomy = taxonomy.split(',')
            for taxonomy_level in taxonomy:
                if taxonomy_level.lower().replace(" ", "_") == searched_expression:
                    is_searched = True
                    break
            if is_searched:
                break
        if is_searched:
            output.write(line.split('\t')[0] + "\n")

output.close()
