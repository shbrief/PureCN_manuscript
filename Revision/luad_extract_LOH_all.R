# Extract LOH 

library(dplyr)
out_dir = "~/data2/PureCN_manuscript/Revision"

data_dir = "/nobackup/16tb_a/CNVworkflow_LUAD/purecn_output/PureCN"
tumor_dir = file.path(data_dir, "tumor_only")
paired_dir = file.path(data_dir, "matching_normal")

# a list of SampleID, which has '_genes.csv' output
sample_dir = list.files(paired_dir)
file_list = sapply(sample_dir, function(x) list.files(file.path(paired_dir, x)))
mut = sapply(file_list, function(x) x[grep("*_genes.csv$", x)])

# subset the 442 samples used in the paper
source("~/Documents/github/PureCN_manuscript/Figures/final/non_dup.R")
sub_ind = which(stringr::str_extract(names(mut), "TCGA.{11}") %in% luad_442)
mut = mut[sub_ind]

## extract LOH info
# C = Segment integer copy number
# M = Minor integer copy number (M + N = C, M ≤ N)
loh_list = list()
for (i in seq_along(mut)) {
    loh_tumor = read.csv(file.path(tumor_dir, names(mut[i]), mut[[i]]))
    loh_tumor_sub = loh_tumor[, c("gene.symbol","chr","start","end","C","M")]
    loh_list[[names(mut[i])]] = loh_tumor_sub
}

# save a table of LOH outputs
saveRDS(loh_list, file = file.path(out_dir, "luad_LOH_all.rds"))

