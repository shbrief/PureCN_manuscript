# Import ABSOLUTE result

# Check "/home/sehyun/Documents/github/PureCN_manuscript/luad/Results/ABSOLUTE/luad_ABSOLUTE_sampleMap.Rmd"
file.dir = "/home/sehyun/Documents/github/PureCN_manuscript/luad/Results/ABSOLUTE"
sampleMap = readRDS(file.path(file.dir, "ABSOLUTE_luad.rds"))
sampleMap = sampleMap[, c("SampleId", "purity", "ploidy", "Genome.doublings", "Subclonal.genome.fraction", "fullname")]
abs_puri_ploi = sampleMap
rm(file.dir, sampleMap)

abs_puri_ploi = abs_puri_ploi[,1:3]
names(abs_puri_ploi) = c("SampleId", "Purity_ABS", "Ploidy_ABS")
abs_puri_ploi = abs_puri_ploi[-which(is.na(abs_puri_ploi$Purity_ABS)), ]
