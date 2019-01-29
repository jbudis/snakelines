library(edgeR)
library(locfit)
library(statmod)

args = commandArgs(trailingOnly=TRUE)

count_file = args[1]
metadata_file = args[2]
output_norm_counts = args[3]
output_diff = args[4]
output_design = args[5]
group_by_attribute = args[6]
if (length(args) == 7) {
    batch_attribute = args[7]
}

counts <- read.delim(count_file, row.names=1, header = T, sep = '\t', check.names = FALSE)
metadata <- read.delim(metadata_file, row.names = 1, header = T, sep = '\t', comment.char = '#', check.names = FALSE)

groups = c()
for (sid in colnames(counts)) {
    groups = append(groups, metadata[sid, group_by_attribute])
}

if (exists("batch_attribute")) {
    litters = c()
    for (sid in colnames(counts)) {
        litters = append(litters, metadata[sid, batch_attribute])
    }
    design <- model.matrix(~litters+groups)
} else {
    design <- model.matrix(~groups)
}

print(typeof(design['litters']))

dge <- DGEList(counts=counts, group=groups)
dge <- calcNormFactors(dge)
dge <- estimateGLMCommonDisp(dge, design)
norm_counts = cpm(dge, normalized.lib.sizes = TRUE)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

sorted_results = topTags(lrt, n=nrow(counts))$table
diff = merge(sorted_results, norm_counts, by=0, all=TRUE)
diff <- diff[order(diff$FDR),]
diff[,2]=2^diff[,2]
colnames(diff)[2] = "FC"
up_down = c()
for(i in diff[,2]){
    if(i < 1){up_down = c(up_down,'down')}
    else{up_down = c(up_down,'up')}
}
FC = c()

i=1
for(f in diff[,2]){
    if(f<1){
        FC = c(FC,diff[i,2]^-1)
    }
    else{FC = c(FC,diff[i,2])}
    i=i+1
}
diff = cbind(diff, FC, up_down)
colnames(diff)[ncol(diff)-1] = 'FCexp^-1'

write.table(diff, file = output_diff, quote=FALSE, sep='\t', row.names=F)
write.table(cbind(sampleid = colnames(counts), design[, c('litters', 'groups')]), file = output_design, quote=FALSE, sep='\t', row.names=F)

# TODO name the row index as "name"
write.table(norm_counts, file = output_norm_counts, quote=FALSE, sep='\t', row.names=F)