# Extract LOH of HLAA, HLAA/B/C

# catalog = "931070"
catalog = "S0293689"

library(dplyr)
library(readr)
library(GenomicRanges)
LOH_dir = "~/wallabe4_backup/github/PureCN_manuscript/Results/LOH"
out_dir = file.path(LOH_dir, "data")

data_dir = file.path("~/wallabe4_backup/data/CNVworkflow", catalog, "purecn_output", paste0(catalog, "_PureCN"))
tumor_dir = file.path(data_dir, "tumor_only")
paired_dir = file.path(data_dir, "matching_normal")

# a list of SampleID, which has '_genes.csv' output
mut = list.files(paired_dir)
mut = mut[grep("*_genes.csv$", mut)]

# Gene Of Interest
goi = c("TP53", "HLA-A", "HLA-B", "HLA-C")   

# extract somatic.rate.ontarget
loh_list = lapply(seq_along(mut), function(i){
  loh_tumor = read.csv(file.path(tumor_dir, mut[i])) %>% GRanges(.)
  loh_tumor = loh_tumor[which(loh_tumor$gene.symbol %in% goi),] 
}) 

loh_tumor_gr = GRangesList(loh_list)
names(loh_tumor_gr) = mut
rm(loh_list)

loh_list = lapply(seq_along(mut), function(i){
  loh_paired = read.csv(file.path(paired_dir, mut[i])) %>% GRanges(.)
  loh_paired = loh_paired[which(loh_paired$gene.symbol %in% goi),]
}) 

loh_paired_gr = GRangesList(loh_list)
names(loh_paired_gr) = mut
rm(loh_list)

if (catalog == "931070") {
  sampleMap_931070 = read.csv("~/wallabe4_backup/data/CNVworkflow/931070/sampleMap_931070.csv")[,-1]
  sampleMap_931070$submitter = stringr::str_extract(sampleMap_931070$barcode, "TCGA.{11}")
  for (i in seq_along(loh_tumor_gr)) {
    fname = gsub("_genes.csv$", ".bam", names(loh_tumor_gr[i]))
    names(loh_tumor_gr)[i] = sampleMap_931070[which(sampleMap_931070$filename == fname),]$submitter
    fname = gsub("_genes.csv$", ".bam", names(loh_paired_gr[i]))
    names(loh_paired_gr)[i] = sampleMap_931070[which(sampleMap_931070$filename == fname),]$submitter
  }
} else if (catalog == "S0293689") {
  sampleMap_S0293689 = read.csv("~/wallabe4_backup/data/CNVworkflow/S0293689/sampleMap_S0293689.csv")[,-1]
  sampleMap_S0293689$submitter = stringr::str_extract(sampleMap_S0293689$barcode, "TCGA.{11}")
  for (i in seq_along(loh_paired_gr)) {
      fname = gsub("_genes.csv$", ".bam", names(loh_tumor_gr[i]))
    names(loh_tumor_gr)[i] = sampleMap_S0293689[which(sampleMap_S0293689$filename == fname),]$submitter
    fname = gsub("_genes.csv$", ".bam", names(loh_paired_gr[i]))
    names(loh_paired_gr)[i] = sampleMap_S0293689[which(sampleMap_S0293689$filename == fname),]$submitter
  }
}

# save a table of LOH outputs (in GRangesList)
saveRDS(loh_tumor_gr, file = file.path(out_dir, paste0(catalog, "_LOH_tumor.rds")))
saveRDS(loh_paired_gr, file = file.path(out_dir, paste0(catalog, "_LOH_paired.rds")))

# extract LOH info
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
write.table(loh_rate, file = file.path(LOH_dir, "extdata", paste0(catalog, "_LOH_tumor_all.tsv")))
