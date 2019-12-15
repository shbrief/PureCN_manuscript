# This script is to format SNV table created from `buildSNVTable.R` script. 
# Input : `snv_list` list object with 
# Output : 


##### Format mutation input ####################################################

for (i in seq_along(snv_list)) {
    
    # subset to only somatic mutation (`ML.SOMATIC = TRUE`)
    snv_list[[i]] = snv_list[[i]][which(snv_list[[i]]$ML.SOMATIC == TRUE),]
    
    # format for dNdSCV
    snv_list[[i]]$Sampleid = stringr::str_extract(snv_list[[i]]$Sampleid, "TCGA.{11}")
    snv_list[[i]] = snv_list[[i]][c("Sampleid", "chr", "start", "REF", "ALT")]
    colnames(snv_list[[i]]) = c("sampleID", "chr", "pos", "ref", "mut")
    
    # remove non-unique calls
    dup = duplicated(snv_list[[i]]$pos)
    snv_list[[i]] = snv_list[[i]][!dup,]
}


##### Function to merge MNPs ###################################################
findMNPs = function(x) {
    # `mnps` object will collect the index of the second nucleotide(s) of MNPs (m = 2)
    mnps = c()
    
    # merge any two, consecutive nucleotides changes
    for (i in 2:nrow(x)) {
        if ((x$pos[i] - x$pos[i-1]) == 1) {  
            x$ref[i-1] = paste0(x$ref[i-1], x$ref[i]) 
            x$mut[i-1] = paste0(x$mut[i-1], x$mut[i])
            mnps = c(mnps, x$pos[i]) 
        }
    }
    
    # remove the second nucleotide(s) of MNPs (m = 2)
    if (!is.null(mnps)) {
        x = x[-which(x$pos %in% mnps),]
    } else {x = x}
    
    return(x)
}


##### Format MNPs (m = 2 cases) ################################################
res = lapply(snv_list, findMNPs)
## an example of merged index from the first sample
# merged_ind = which(!snv_list[[1]]$pos %in% res[[1]]$pos)
# snv_list[[1]]$pos[sort(c(merged_ind-1, merged_ind))]


##### Longer MNPs (m != 2 cases) ###############################################
res2 = lapply(res, findMNPs)

merge1 = sapply(res, function(x) {nrow(x)})
merge2 = sapply(res2, function(x) {nrow(x)})

# If the second application of the function,`findMNPs`, merges rows again, there 
# is likely a MNP(s) longer than dinucleotides. Please manually check them. 
if (sum(merge1 != merge2) != 0) {
    longer_mnps = names(merge1)[merge1 != merge2]
    print("Double check the samples in the list, called 'longer_mnps'")
}

rm(merge1, merge2, res, res2)
