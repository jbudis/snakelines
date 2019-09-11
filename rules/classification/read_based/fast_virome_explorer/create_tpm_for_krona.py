#!/usr/bin/python3

import sys
import pandas as pd

def create_working_accessories(input_file):
    """
    Creates working directory and file based on absolute path of input
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


def create_tpm_table(input_file, dictionary, abundance_file, krona_table):
    """
    Creates tpm table for krona Tools
    :param input_file:
    :param dictionary:
    :param abundance_file:
    :param krona_table:
    """
    load_input_file = pd.read_csv(input_file, header = None, sep='\t')
    for i in range(len(load_input_file[0])):
        if i != 0:
            dictionary[load_input_file[0][i]] = None

    # load abundace file and search into this file ncbi_id from input file, then save ncbi_id to dictionary
    load_abundance_file = pd.read_csv(abundance_file, header = None, sep='\t')
    for i in dictionary:
        if i != 0:
            matched_line = load_abundance_file.loc[load_abundance_file[0] == i] 
            dictionary[i] = float(matched_line[4])
    
    # duplicate input file, but with new name of file and load it
    load_input_file.to_csv(krona_table, sep = '\t')
    load_krona_table = pd.read_csv(krona_table, header = None, sep='\t')
    
    # replace values in original krona table with tpm values from abundance file
    for ncbi_id, index in zip(load_krona_table[1], range(len(load_krona_table[1]))):
        if ncbi_id[0] == '#' or ncbi_id[0] == '0':
            pass
        else:
            load_krona_table.loc[index, len(load_krona_table.columns)-1] = dictionary[ncbi_id]

    # remove first column, rename column's names and write this table to file
    load_krona_table = load_krona_table.drop(0,1)
    load_krona_table.rename(columns={1: load_krona_table[1][1], 2: load_krona_table[2][1], 3 :load_krona_table[3][1], 4: load_krona_table[4][1] }, inplace=True)
    load_krona_table = load_krona_table.drop([0, 1])

    if load_krona_table.isnull().values.any() == True:
        load_krona_table['VirusName'] = load_krona_table['#VirusIdentifier']
        load_krona_table['kingdom;phylum;class;order;family;genus;species'] = load_krona_table['#VirusIdentifier']

    load_krona_table.to_csv(krona_table, na_rep = 'NaN',  sep = '\t', index = False)


# creating of intial variables
input_file = sys.argv[1]
working_dir, working_file = create_working_accessories(input_file)[0], create_working_accessories(input_file)[1]
dictionary = dict()
abundance_file = '{}/abundance.tsv'.format(working_dir)
krona_table = '{}/FastViromeExplorer-final-sorted-abundance_checked_tpm.tsv'.format(working_dir)

# calling main function of script
create_tpm_table(input_file, dictionary, abundance_file, krona_table)