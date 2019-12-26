# Extract SNV info from PureCN
# - input: `i` (= runtype), `out_dir`   
# - output: `snv_list` 

##### OVC 901070 ###############################################################
# PureCN output file list
data_dir = "/nobackup/16tb_b/CNVworkflow/931070/purecn_output/931070_PureCN"
sample_dir = file.path(data_dir, runtype)
sample_name = list.files(sample_dir)
file_list = sapply(sample_name, function(x) list.files(file.path(sample_dir, x)))
mut = sapply(file_list, function(x) x[grep("*_variants.csv$", x)])

# subset the 236 samples used in the paper
source("~/Documents/github/PureCN_manuscript/Figures/Final_Figures/non_dup.R")
sub_ind = which(stringr::str_extract(names(mut), "TCGA.{11}") %in% ovc_236)
mut = mut[sub_ind]

# SNV table in a dNdSCV input format
snv_list = list()
for (i in seq_along(mut)) {
    if (mut[i] != "character(0)") {
        # p <- predictSomatic output
        p = read.csv(file.path(data_dir, runtype, names(mut[i]), mut[[i]]))
        
        # remove variants in germline databases
        min.prior.somatic = 0.1 
        max.prior.somatic = 1
        p = p[p$prior.somatic >= min.prior.somatic & 
                  p$prior.somatic <= max.prior.somatic,]
        # keep only somatic in targeted regions  
        # p = p[p$ML.SOMATIC & p$on.target == 1, ]
        p = p[p$POSTERIOR.SOMATIC > 0.8 & p$on.target == 1, ]
        # remove flagged         
        p = p[!p$FLAGGED,]
        
        ##### subset only the columns required for dNdSCV
        snv_table_sub = p[, c("Sampleid","chr","start","end","REF","ALT")]
        snv_list[[names(mut[i])]] = snv_table_sub
    } else {
        next
    }
}

saveRDS(snv_list, file = file.path(out_dir, "data", paste0("OVC_931070_mutTable_", runtype, ".rds")))
rm(sample_dir, sample_name, file_list, mut, sub_ind, snv_list, snv_table_sub, p)

##### OVC S0293689 #############################################################
# PureCN output file list
data_dir = "/nobackup/16tb_b/CNVworkflow/S0293689/purecn_output/S0293689_PureCN"
sample_dir = file.path(data_dir, runtype)
sample_name = list.files(sample_dir)
file_list = sapply(sample_name, function(x) list.files(file.path(sample_dir, x)))
mut = sapply(file_list, function(x) x[grep("*_variants.csv$", x)])

# match filenames with sample names
mani = read.table("~/data2/PureCN_manuscript/Data/manifest/ovc_manifest.tsv")
x = as.data.frame(matrix(NA, nrow = length(mut), ncol=3))
colnames(x) = c("file_name", "barcode", "sample_name")
x$file_name = mut
for (i in 1:nrow(x)) {
    y = gsub("_variants.csv", ".bam", x$file_name[i])
    if (!identical(y, "character(0)")) {
        z = mani$barcode[which(mani$filename == y)]
        x$barcode[i] = as.character(z)
        x$sample_name[i] = stringr::str_extract(x$barcode[i], "TCGA.{11}")
    } else {
        next
    }
}

# subset the 236 samples used in the paper
sub_ind = which(x$sample_name %in% ovc_236)
mut = x$file_name[sub_ind]

# subset SNV info 
snv_list = list()
for (i in seq_along(mut)) {
    if (mut[i] != "character(0)") {
        # p <- predictSomatic output
        p = read.csv(file.path(data_dir, runtype, names(mut[i]), mut[[i]]))
        
        # remove variants in germline databases
        min.prior.somatic = 0.1 
        max.prior.somatic = 1
        p = p[p$prior.somatic >= min.prior.somatic & 
                  p$prior.somatic <= max.prior.somatic,]
        
        # keep only somatic in targeted regions  
        p = p[p$ML.SOMATIC & p$on.target == 1, ]
        # p = p[p$POSTERIOR.SOMATIC > 0.8 & p$on.target == 1, ]
        
        # remove flagged         
        p = p[!p$FLAGGED,]
        
        ##### subset only the columns required for dNdSCV
        snv_table_sub = p[, c("Sampleid","chr","start","end","REF","ALT")]
        snv_list[[names(mut[i])]] = snv_table_sub
    } else {
        next
    }
}

# changes the name
for (i in seq_along(snv_list)) {
    ind = which(x$file_name == paste0(names(snv_list)[i], "_variants.csv"))
    sampleId = x$sample_name[ind]
    names(snv_list)[i] = sampleId
    snv_list[[i]]$Sampleid = sampleId
}

saveRDS(snv_list, file = file.path(out_dir, "data", paste0("OVC_S0293689_mutTable_", runtype, ".rds")))
