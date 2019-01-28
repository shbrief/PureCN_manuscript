# Extract LOH of HLAA, HLAA/B/C

# catalog = "931070"
catalog = "S0293689"

library(dplyr)
LOH_dir = "~/Documents/github/PureCN_manuscript/Results/LOH"
out_dir = file.path(LOH_dir, "extdata")

data_dir = file.path("/data/16tb/CNVworkflow", catalog, "purecn_output", paste0(catalog, "_PureCN"))
tumor_dir = file.path(data_dir, "tumor_only")
paired_dir = file.path(data_dir, "matching_normal")

# a list of SampleID, which has '_mutation_burden.csv' output
mut = list.files(paired_dir)
mut = mut[file.exists(file.path(paired_dir, mut, paste0(mut, "_genes.csv")))]

# Gene Of Interest
goi = c("TP53", "HLA-A", "HLA-B", "HLA-C")   

# extract somatic.rate.ontarget
loh_list = lapply(seq_along(mut), function(i){
  fname = mut[i]
  loh_fname = file.path(fname, paste0(fname, "_genes.csv"))
  loh_tumor = read.csv(file.path(tumor_dir, loh_fname)) %>% GRanges(.)
  loh_tumor = loh_tumor[which(loh_tumor$gene.symbol %in% goi),] 
}) 

loh_tumor_gr = GRangesList(loh_list)
names(loh_tumor_gr) = mut
rm(loh_list)

loh_list = lapply(seq_along(mut), function(i){
  fname = mut[i]
  loh_fname = file.path(mut[i], paste0(mut[i], "_genes.csv"))
  loh_paired = read.csv(file.path(paired_dir, loh_fname)) %>% GRanges(.)
  loh_paired = loh_paired[which(loh_paired$gene.symbol %in% goi),]
}) 

loh_paired_gr = GRangesList(loh_list)
names(loh_paired_gr) = mut
rm(loh_list)

if (catalog == "931070") {
  sampleMap_931070 = read_csv("/data/16tb/CNVworkflow/931070/sampleMap_931070.csv")
  sampleMap_931070$submitter = stringr::str_extract(sampleMap_931070$barcode, "TCGA.{11}")
  for (i in seq_along(loh_tumor_gr)) {
    fname = paste0(names(loh_tumor_gr)[i], ".bam")
    names(loh_tumor_gr)[i] = sampleMap_931070[which(sampleMap_931070$filename == fname),]$submitter
    fname = paste0(names(loh_paired_gr)[i], ".bam")
    names(loh_paired_gr)[i] = sampleMap_931070[which(sampleMap_931070$filename == fname),]$submitter
  }
} else if (catalog == "S0293689") {
  sampleMap_S0293689 = read_csv("/data/16tb/CNVworkflow/S0293689/sampleMap_S0293689.csv")
  sampleMap_S0293689$submitter = stringr::str_extract(sampleMap_S0293689$barcode, "TCGA.{11}")
  for (i in seq_along(loh_paired_gr)) {
    fname = paste0(names(loh_tumor_gr)[i], ".bam")
    names(loh_tumor_gr)[i] = sampleMap_S0293689[which(sampleMap_S0293689$filename == fname),]$submitter
    fname = paste0(names(loh_paired_gr)[i], ".bam")
    names(loh_paired_gr)[i] = sampleMap_S0293689[which(sampleMap_S0293689$filename == fname),]$submitter
  }
}

# save a table of LOH outputs
saveRDS(loh_tumor_gr, file = file.path(out_dir, paste0(catalog, "_LOH_tumor.rds")))
saveRDS(loh_paired_gr, file = file.path(out_dir, paste0(catalog, "_LOH_paired.rds")))

# # extract somatic.rate.ontarget
# loh_list = lapply(seq_along(mut), function(i){
#   loh_fname = file.path(mut[i], paste0(mut[i], "_genes.csv"))
#   loh_tumor = read.csv(file.path(tumor_dir, loh_fname))
#   loh_paired = read.csv(file.path(paired_dir, loh_fname))
#   
#   data.frame(sampleID = mut[i], 
#             HLAA_LOH_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-A"),]$loh, 
#             HLAA_LOH_p = loh_paired[which(loh_tumor$gene.symbol == "HLA-A"),]$loh,
#             HLAB_LOH_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-B"),]$loh, 
#             HLAB_LOH_p = loh_paired[which(loh_tumor$gene.symbol == "HLA-B"),]$loh,
#             HLAC_LOH_t = loh_tumor[which(loh_tumor$gene.symbol == "HLA-C"),]$loh, 
#             HLAC_LOH_p = loh_paired[which(loh_tumor$gene.symbol == "HLA-C"),]$loh
#             )
# }) 
# 
# loh_rate = data.frame()
# for (i in seq_along(loh_list)) {loh_rate = rbind(loh_rate, loh_list[[i]])}
# 
# # save a table of LOH outputs
# write.table(loh_rate, file = file.path(out_dir, paste0(catalog, "_HLA_LOH.tsv")))