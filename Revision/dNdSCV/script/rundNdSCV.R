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

# A list of the second index of SNV next to each other
findMNPs = function(x) {
    mnps = c()
    for (i in 2:nrow(x)) {
        if ((x$pos[i] - x$pos[i-1]) == 1) {  
            mnps = c(mnps, i) 
        }
    }
    return(mnps)
}


test = sapply(snv_list, findMNPs)
test[[1]] -> test_ind
snv_list[[1]][sort(c(test_ind, test_ind-1)),]



diffMNPs = sapply(test, function(x) {diff(x)})
a = sapply(diffMNPs, function(x) {1 %in% x})
b = sapply(diffMNPs, function(x) {1 %in% x & 2 %in% x})

sub = diffMNPs[c(1,4,7)]
sapply(sub, function(x) {1 %in% x})
sapply(sub, function(x) {1 %in% x & 2 %in% x})
c = sapply(sub, function(x) {c(1,2) %in% x})


for (j in seq_along(diffMNPs)) {
    mnp_length = c()
    for (i in 1:6) {
        sTable = c(1:i) %in% diffMNPs[[j]]
        if (sum(sTable) == i) {
            mnp_length = c(mnp_length, i)
            return(mnp_length)
        }
    }
    k = max(mnp_length)
    print(paste(j, "th sample,", names(diffMNPs)[j], ", has the maximum length of MNP =", k))
}





minMNPs = sapply(test, function(x) {min(diff(x))})
triMNPs = sapply(test, function(x) {length(intersect(c(1:6), (diff(x)))) == 6})



triMNPs = sapply(test, function(x) {1 %in% (diff(x)) & 2 %in% (diff(x))})
sum(triMNPs)

##### Combine all mutations ####################################################
mut_all = as.data.frame(matrix(NA, ncol = 5))
colnames(mut_all) = colnames(snv_list[[1]])
for (i in seq_along(snv_list)) {mut_all = rbind(mut_all, snv_list[[i]])}
mut_all = mut_all[-1,]



##### Run dNdSCV ###############################################################
dndsout = dndscv(mut_all, refdb = file.path(wd, "GRCh38_refcds.rda"), cv = NULL)







# ##### Run dNdSCV #############################################################
# runDNDS = function(i) {
#     tryCatch({
#         sample_name = names(snv_list[i])
#         dndsout = dndscv(snv_list[[i]], 
#                          refdb = file.path(wd, "GRCh38_refcds.rda"),
#                          max_muts_per_gene_per_sample = 3,   # default
#                          max_coding_muts_per_sample = 3000,   # default
#                          cv = NULL)
#         write.csv(dndsout$sel_cv, file = file.path(dnds_out_dir, "sel_cv", paste0(sample_name, "_sel_cv.csv")))
#         saveRDS(dndsout, file = file.path(dnds_out_dir, paste0(sample_name, ".rds")))
#     }, error = function(e) {cat("ERROR with ", paste0("sample #", i, "of", sample_name), "\n")})
# }
# 
# sample_name = names(snv_list)
# sample_name = sample_name[!file.exists(file.path(dnds_out_dir, paste0(sample_name, ".rds")))]
# snv_list_sub = snv_list[sample_name]
# 
# bplapply(seq_along(snv_list_sub), runDNDS, BPPARAM = MulticoreParam(24))
