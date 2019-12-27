# Extract SNV info from PureCN
# - input: `i` (= runtype), `out_dir`, `data_dir`   
# - output: `snv_list` 

##### PureCN output file list ##################################################
# a list of SampleID, which has '_variants.csv' output
sample_dir = list.files(data_dir)
file_list = sapply(sample_dir, function(x) list.files(file.path(data_dir, x)))

if (runtype == "tumor_only") {
    mut = sapply(file_list, function(x) x[grep("*.rds$", x)])
} else if (runtype == "matching_normal") {
    mut = sapply(file_list, function(x) x[grep("*_variants.csv$", x)])
}

# subset the 442 samples used in the paper
source("~/Documents/github/PureCN_manuscript/Figures/Final_Figures/non_dup.R")
sub_ind = which(stringr::str_extract(names(mut), "TCGA.{11}") %in% luad_442)
mut = mut[sub_ind]


##### SNV table in a dNdSCV input format #######################################
source("~/data2/PureCN_manuscript/Revision/dNdSCV/script/rescale_priors.R")
vcf_all = readVcf("~/data2/PureCN_manuscript/Revision/dNdSCV/data/CosmicCodingMuts.vcf.bgz", "hg38")
vcf = vcf_all[info(vcf_all)$CNT >= 20 & !info(vcf_all)$SNP]
seqlevelsStyle(vcf) = "UCSC"

snv_list = list()
for (i in seq_along(mut)) {
    
    if (runtype == "matching_normal") {
        p = read.csv(file.path(data_dir, names(mut[i]), mut[[i]]))
    } else if (runtype == "tumor_only") {
        # rescale prior
        x = readCurationFile(file.path(data_dir, names(mut[i]), mut[[i]]))
        p = rescale_overlaps(x, vcf)
        p_old = predictSomatic(x)
        p = data.frame(Sampleid = x$input$sampleid, p)
    }
    
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
    
    # # add hotspots back
    # p = dplyr::union(p, hotspots)
    
    ##### subset only the columns required for dNdSCV
    snv_table_sub = p[, c("Sampleid","chr","start","end","REF","ALT")]
    snv_list[[names(mut[i])]] = snv_table_sub
}

saveRDS(snv_list, file = file.path(out_dir, "data", paste0("LUAD_mutTable_", runtype, "_rescalePriors.rds")))

# out_dir = "~/data2/PureCN_manuscript/Revision/dNdSCV"
# snv_list = readRDS(file.path(out_dir, paste0("luad_snv", runtype, ".rds")))
