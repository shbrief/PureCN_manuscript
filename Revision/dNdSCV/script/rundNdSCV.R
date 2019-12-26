##### Combine all samples ######################################################
mut_all = as.data.frame(matrix(NA, ncol = 5))
colnames(mut_all) = colnames(snv_list[[1]])
for (i in seq_along(snv_list)) {mut_all = rbind(mut_all, snv_list[[i]])}
mut_all = mut_all[-1,]


##### Run dNdSCV ###############################################################
dndsout = dndscv(mut_all, refdb = file.path(wd, "reference/GRCh38_refcds.rda"), cv = NULL)
