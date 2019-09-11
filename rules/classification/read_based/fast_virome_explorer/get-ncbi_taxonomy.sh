#!/bin/bash
#title           :get-ncbi_id.sh
#description     :This script assign taxonomy to ncbi_id of organism.
#author          :smolak
#date            :2019/08/20
#version         :non-version    
#usage           :bash get-ncbi_id.sh {ncbi_id}
#notes           :necesarry to have installed blastdbcmd.
#bash_version    :4.1.5(1)-release
#==============================================================================
NCBI_ID=$1
stdbuf -oL blastdbcmd -db /data/genome/metagenome/blast/refseq/19-01-07/refseq_genomic -entry $NCBI_ID -outfmt '%T' |
  while IFS= read -r LINE
  do
    cat /data/genome/taxonomy/lineages-2019-02-20.csv | egrep "^$LINE," | cut -d ',' -f 2-8 | sed 's/,/;/g'
  done
