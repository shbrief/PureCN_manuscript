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
