abs_gl = readRDS("/home/sehyun/Documents/github/subtypeHeterogeneity/inst/extdata/ABSOLUTE_grangeslist.rds")
puri_ploi = read.csv("/home/sehyun/data2/PureCN_manuscript/Figures/Final_Tables/Table1_puri_ploi.csv")[,-1]

diploid_ind = which(puri_ploi$Ploidy_ABS < 2.1 & puri_ploi$Ploidy_ABS > 1.9)
diploid_samples = puri_ploi[diploid_ind,]$SampleId

for (sample in diploid_samples) {
    x = table(paste0(abs_gl[[sample]]$Modal_HSCN_1, "-",abs_gl[[sample]]$Modal_HSCN_2))
    print(x)
}





# CreateLOHmetadata = function(x) {
#     for (i in seq_along(x)) {
#         loh = c(x$Modal_HSCN_1[i], x$Modal_HSCN_2[i])
#         diff = max(loh) - min(loh)
#         x$LOH[i] = paste0(diff, "-", min(loh))
#     }
#     return(x)
# }