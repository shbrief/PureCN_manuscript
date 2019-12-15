# Extract SNV info from PureCN


##### PureCN output file list ##################################################
# a list of SampleID, which has '_variants.csv' output
sample_dir = list.files(data_dir)
file_list = sapply(sample_dir, function(x) list.files(file.path(data_dir, x)))
mut = sapply(file_list, function(x) x[grep("*_variants.csv$", x)])

# subset the 442 samples used in the paper
source("~/Documents/github/PureCN_manuscript/Figures/Final_Figures/non_dup.R")
sub_ind = which(stringr::str_extract(names(mut), "TCGA.{11}") %in% luad_442)
mut = mut[sub_ind]


##### SNV table in a dNdSCV input format #######################################
snv_list = list()
for (i in seq_along(mut)) {
    snv_table = read.csv(file.path(data_dir, names(mut[i]), mut[[i]]))
    snv_table_sub = snv_table[, c("Sampleid","chr","start","end","REF","ALT", "ML.SOMATIC")]
    snv_list[[names(mut[i])]] = snv_table_sub
}

saveRDS(snv_list, file = file.path(out_dir, paste0("LUAD_mutTable_", runtype, ".rds")))

# out_dir = "~/data2/PureCN_manuscript/Revision/dNdSCV"
# snv_list = readRDS(file.path(out_dir, paste0("luad_snv", runtype, ".rds")))
