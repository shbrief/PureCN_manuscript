# Run dNdSCV 
# [input] snv_list : a list of data frames with 5 columns for dNdScv ("sampleID", "chr", "pos", "ref", "mut")
# [output] dndsout : a list of dNdSCV outputs

wd = "~/data2/PureCN_manuscript/Revision/dNdSCV"

for (i in seq_along(snv_list)) {
    tryCatch({
        sapmle_name = names(snv_list[i])
        dndsout = dndscv(snv_list[[i]], 
                         refdb = file.path(wd, "GRCh38_refcds.rda"),
                         # max_muts_per_gene_per_sample = 5,   # default is 3
                         # max_coding_muts_per_sample = 20000,   # default is 3000
                         cv = NULL)
        write.csv(dndsout, file = file.path(wd, "output/LUAD", sample_name, ".csv"))
    }, error = function(e) {cat("ERROR with ", paste0("sample #",i, "of", sample_name), "\n")})
}