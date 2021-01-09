#!/usr/bin/env Rscript

library(VirFinder)

fasta <- commandArgs(TRUE)[1]
out <- commandArgs(TRUE)[2]

predResult <- VF.pred(fasta)
# predResult$qvalue <- VF.qvalue(predResult$pvalue)
write.table(predResult, out, sep='\t', row.names=F, quote=F)