#!/usr/bin/python3

#import modules

import sys
import numpy as np
import pandas as pd
from glob import glob
import os
import csv


#function that creates working directory and file from input absolute path to file

def create_working_accessories(input_file):
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



#function below replace 'Nan' values in first and second column with #VirusIdentifier values

def fillna_with_virusname(input_file):


    load_krona_table = pd.read_csv(input_file, sep = '\t')                                                              #import krona table

    if load_krona_table.isnull().values.any() == True:                                                                  #If table contains at least one 'NaN' value, condition statement is done
        load_krona_table['VirusName'] = load_krona_table['#VirusIdentifier']                                             #values in column 'VirusName' are replaced with values of column '#VirusIdentifier'
        load_krona_table['kingdom;phylum;class;order;family;genus;species'] = 6 * 'Unclassified;' + load_krona_table['#VirusIdentifier']      #values in column 'kingdom;...' are replaced with values of column '#VirusIdentifier' as genus other taxonomic levels are filled as 'Unclassified'
        load_krona_table.to_csv(output_file, na_rep = 'NaN', index = False, sep = '\t')

    else:
        load_krona_table.to_csv(output_file, index = False, sep = '\t')
        

#creating of intial variables

input_file = sys.argv[1]                                                                                                #as input serves first shell argument
working_dir, working_file = create_working_accessories(input_file)[0], create_working_accessories(input_file)[1]
output_file = list()
output_file += ['.'.join(input_file.split('.')[:-1]), '-checked.', input_file.split('.')[-1]]                           #divide string into substring, dividing character is '.' and then take all substring except last from list and merge them together with string '-chcecked' and last substring from list
output_file = ''.join(output_file)

#calling main function of script

fillna_with_virusname(input_file)                                                                                         #function calling

