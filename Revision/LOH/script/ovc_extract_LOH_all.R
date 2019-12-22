# Extract LOH 

library(dplyr)
out_dir = "~/data2/PureCN_manuscript/Revision"

# 931070
data_dir = "/nobackup/16tb_b/CNVworkflow/931070/purecn_output/931070_PureCN"
tumor_dir_931070 = file.path(data_dir, "tumor_only")
sample_dir = list.files(tumor_dir_931070)
file_list = sapply(sample_dir, function(x) list.files(file.path(tumor_dir_931070, x)))
mut_931070 = sapply(file_list, function(x) x[grep("*_genes.csv$", x)])

# 931070
data_dir = "/nobackup/16tb_b/CNVworkflow/S0293689/purecn_output/S0293689_PureCN"
tumor_dir_S0293689 = file.path(data_dir, "tumor_only")
sample_dir = list.files(tumor_dir_S0293689)
file_list = sapply(sample_dir, function(x) list.files(file.path(tumor_dir_S0293689, x)))
mut_S0293689 = sapply(file_list, function(x) x[grep("*_genes.csv$", x)])

mut = c(mut_931070, mut_S0293689)

# match filenames with sample names
mani = read.table("~/data2/PureCN_manuscript/Data/manifest/ovc_manifest.tsv")
x = as.data.frame(matrix(NA, nrow = length(mut), ncol=3))
colnames(x) = c("file_name", "barcode", "sample_name")
x$file_name = mut
for (i in 1:nrow(x)) {
    y = gsub("_genes.csv", ".bam", x$file_name[i])
    if (!identical(y, "character(0)")) {
        z = mani$barcode[which(mani$filename == y)]
        x$barcode[i] = as.character(z)
        x$sample_name[i] = stringr::str_extract(x$barcode[i], "TCGA.{11}")
    } else {
        next
    }
}
# saveRDS(x, "~/data2/PureCN_manuscript/Revision/LOH/data/ovc_name_convert.rds")

# subset the 236 samples used in the paper
source("~/data2/PureCN_manuscript/Figures/Final_Figures/non_dup.R")
ovc_236_ind = which(x$sample_name %in% ovc_236)
mut = x$file_name[ovc_236_ind]




## extract LOH info
# C = Segment integer copy number
# M = Minor integer copy number (M + N = C, M ≤ N)
loh_list = list()
for (i in seq_along(mut)) {
    file_path = file.path(tumor_dir_931070, names(mut[i]), mut[[i]])
    if (!file.exists(file_path)) {
        file_path = file.path(tumor_dir_S0293689, names(mut[i]), mut[[i]])
    }
    loh_tumor = read.csv(file_path)
    loh_tumor_sub = loh_tumor[, c("gene.symbol","chr","start","end","C","M")]
    loh_list[[names(mut[i])]] = loh_tumor_sub
}

# save a table of LOH outputs
saveRDS(loh_list, file = file.path(out_dir, "LOH/ovc_LOH_all.rds"))

