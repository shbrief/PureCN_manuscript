##### Check all the MNPs are merged properly ###################################
res = lapply(snv_list, findMNPs)

merge1 = sapply(snv_list, function(x) {nrow(x)})
cfmerge2 = sapply(res, function(x) {nrow(x)})

if (sum(merge1 != merge2) == 0) {
    print("MNPs are merged properly. Good to go!")
} else {
    check_mnps = names(merge1)[merge1 != merge2]
    print("Double check the samples in the list, called 'check_mnps'")    
    stop
}

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
