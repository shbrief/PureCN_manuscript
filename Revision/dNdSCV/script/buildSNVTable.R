# Extract SNV info from PureCN
# - input: `i` (= runtype), `out_dir`, `data_dir`   
# - output: `snv_list` 

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
    # p <- predictSomatic output
    p = read.csv(file.path(data_dir, names(mut[i]), mut[[i]]))
    
    # remove variants in germline databases
    min.prior.somatic = 0.1 
    max.prior.somatic = 1
    p = p[p$prior.somatic >= min.prior.somatic & 
              p$prior.somatic <= max.prior.somatic,]
    # keep only somatic in targeted regions  
    p = p[p$ML.SOMATIC & p$on.target == 1, ]
    # remove flagged         
    p = p[!p$FLAGGED,]
    
    ##### subset only the columns required for dNdSCV
    snv_table_sub = p[, c("Sampleid","chr","start","end","REF","ALT")]
    snv_list[[names(mut[i])]] = snv_table_sub
}

saveRDS(snv_list, file = file.path(out_dir, paste0("LUAD_mutTable_", runtype, ".rds")))

# out_dir = "~/data2/PureCN_manuscript/Revision/dNdSCV"
# snv_list = readRDS(file.path(out_dir, paste0("luad_snv", runtype, ".rds")))
