#!/usr/bin/python3

import sys
import pandas as pd


def create_working_accessories(input_file):
    """
    Creates working directory and file according to absolute path of input
    :param input_file:
    :return:
    """
    working_dir = list()
    working_file = list()
    no = input_file.count('/')
    
    for i in input_file:
        if i == '/':
            no -= 1

        if no != 0:
            working_dir.append(i)
            working_file.append(i)
                   
        if no == 0:
            working_file.append(i)
        
    working_dir = ''.join(working_dir)
    working_file = ''.join(working_file)
    
    return working_dir, working_file


def fillna_with_virusname(input_file):
    """
    Replaces 'Nan' values in first and second column with #VirusIdentifier values
    :param input_file:
    """
    load_krona_table = pd.read_csv(input_file, sep = '\t')

    # If table contains at least one 'NaN' value, condition statement is done
    if load_krona_table.isnull().values.any() == True:
        load_krona_table['VirusName'] = load_krona_table['#VirusIdentifier']
        load_krona_table['kingdom;phylum;class;order;family;genus;species'] = 6 * 'Unclassified;' + load_krona_table['#VirusIdentifier']
        load_krona_table.to_csv(output_file, na_rep = 'NaN', index = False, sep = '\t')
    else:
        load_krona_table.to_csv(output_file, index = False, sep = '\t')
        

# Creating of intial variables
input_file = sys.argv[1]  # As input serves first shell argument
working_dir, working_file = create_working_accessories(input_file)[0], create_working_accessories(input_file)[1]
output_file = list()
output_file += ['.'.join(input_file.split('.')[:-1]), '_checked.', input_file.split('.')[-1]]  # Creates new output name
output_file = ''.join(output_file)

# Calling main function of script
fillna_with_virusname(input_file)