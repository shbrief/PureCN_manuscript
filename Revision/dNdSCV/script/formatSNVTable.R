# This script is to format SNV table created from `buildSNVTable.R` script. 
# Input : `snv_list` list object with 
# Output : 


##### Format mutation input ####################################################

for (i in seq_along(snv_list)) {

    # format for dNdSCV
    snv_list[[i]]$Sampleid = stringr::str_extract(snv_list[[i]]$Sampleid, "TCGA.{11}")
    snv_list[[i]] = snv_list[[i]][c("Sampleid", "chr", "start", "REF", "ALT")]
    colnames(snv_list[[i]]) = c("sampleID", "chr", "pos", "ref", "mut")
    
    # ref and mut as 'character' for paste
    snv_list[[i]]$ref = as.character(snv_list[[i]]$ref)
    snv_list[[i]]$mut = as.character(snv_list[[i]]$mut)
    
    # remove non-unique calls
    dup = duplicated(snv_list[[i]]$pos)
    snv_list[[i]] = snv_list[[i]][!dup,]
}


##### Function to merge MNPs ###################################################
# remove sample with less than 3 SNVs and manually check them
num_snv = sapply(snv_list, function(x) {nrow(x)})
snv_ind = which(num_snv <= 3)

if (length(snv_ind) != 0) {
    mannual_snv_check = snv_list[snv_ind]
    print("Manually check MNP in the list,`mannual_snv_check`")
    snv_list[snv_ind] = NULL
}

# merge MNPs
source('~/data2/PureCN_manuscript/Revision/dNdSCV/script/findMNPs.R')


##### Format MNPs (m = 2 cases) ################################################
res = lapply(snv_list, findMNPs)



# ##### How to check MNPs ########################################################
# i = 7   # random example 
# res[[i]][which(nchar(res[[i]]$ref) != 1),]   # subset of MNPs in i-th sample
# 
# ori = sapply(snv_list, function(x) {nrow(x)})
# merged = sapply(res, function(x) {nrow(x)})
# names(ori)[ori == merged]   # name of samples with no MNPs
