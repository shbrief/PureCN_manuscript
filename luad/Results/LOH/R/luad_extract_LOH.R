# Extract LOH of HLAA, HLAA/B/C

library(dplyr)
LOH_dir = "~/wallabe4_backup/github/PureCN_manuscript/luad/Results/LOH"
out_dir = file.path(LOH_dir, "extdata")

data_dir = "~/wallabe4_backup/data/CNVworkflow_LUAD/purecn_output/PureCN"
tumor_dir = file.path(data_dir, "tumor_only")
paired_dir = file.path(data_dir, "matching_normal")

# a list of SampleID, which has '_genes.csv' output
mut = list.files(paired_dir)
mut = mut[grep("*_genes.csv$", mut)]

## extract LOH info
# C = Segment integer copy number
# M = Minor integer copy number (M + N = C, M ≤ N)
loh_list = lapply(seq_along(mut), function(i){
  loh_tumor = read.csv(file.path(tumor_dir, mut[i]))
  data.frame(sampleID = mut[i],
             HLAA_C_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-A"),]$C,
             HLAA_M_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-A"),]$M,
             HLAA_LOH_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-A"),]$loh,
             HLAB_C_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-B"),]$C,
             HLAB_M_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-B"),]$M,
             HLAB_LOH_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-B"),]$loh,
             HLAC_C_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-C"),]$C,
             HLAC_M_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-C"),]$M,
             HLAC_LOH_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-C"),]$loh,
             TP53_C_t = loh_tumor[which(loh_tumor$gene.symbol == "TP53"),]$C,
             TP53_M_t = loh_tumor[which(loh_tumor$gene.symbol == "TP53"),]$M,
             TP53_LOH_t = loh_tumor[which(loh_tumor$gene.symbol == "TP53"),]$loh
  )
})

loh_rate = data.frame()
for (i in seq_along(loh_list)) {loh_rate = rbind(loh_rate, loh_list[[i]])}

# save a table of LOH outputs
write.table(loh_rate, file = file.path(out_dir, "luad_LOH.tsv"))


# # Gene Of Interest
# goi = c("TP53", "HLA-A", "HLA-B", "HLA-C")   
# 
# # extract loh info
# loh_list = lapply(seq_along(mut), function(i){
#   fname = mut[i]
#   loh_fname = file.path(fname, paste0(fname, "_genes.csv"))
#   loh_tumor = read.csv(file.path(tumor_dir, loh_fname)) %>% GRanges(.)
#   loh_tumor = loh_tumor[which(loh_tumor$gene.symbol %in% goi),] 
# }) 
# 
# loh_tumor_gr = GRangesList(loh_list)
# names(loh_tumor_gr) = mut
# rm(loh_list)
# 
# loh_list = lapply(seq_along(mut), function(i){
#   fname = mut[i]
#   loh_fname = file.path(mut[i], paste0(mut[i], "_genes.csv"))
#   loh_paired = read.csv(file.path(paired_dir, loh_fname)) %>% GRanges(.)
#   loh_paired = loh_paired[which(loh_paired$gene.symbol %in% goi),]
# }) 
# 
# loh_paired_gr = GRangesList(loh_list)
# names(loh_paired_gr) = mut
# rm(loh_list)
# 
# # rename
# names(loh_tumor_gr) = stringr::str_extract(names(loh_tumor_gr), "TCGA.{11}")
# names(loh_paired_gr) = stringr::str_extract(names(loh_paired_gr), "TCGA.{11}")
# 
# # save a table of LOH outputs
# saveRDS(loh_tumor_gr, file = file.path(out_dir, "luad_LOH_tumor.rds"))
# saveRDS(loh_paired_gr, file = file.path(out_dir, "luad_LOH_paired.rds"))

