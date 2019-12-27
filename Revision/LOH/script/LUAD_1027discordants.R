# 1027 discordant LOH in LUAD
# "2-1" from PureCN and "2-2" from ABSOLUTE

# SampleId
loh_summary = readRDS("~/data2/PureCN_manuscript/Revision/LOH/data/LUAD_LOH_summary.rds")
discordant = sapply(loh_summary, function(x) {sum(x$loh.pcn == "2-1" & x$loh.abs == "2-2")})
discordant_sampleId = names(which(discordant != 0))
discordant_sampleId

# Purity-Ploidy summary
puri_ploi = read.csv("~/data2/PureCN_manuscript/Figures/Final_Tables/Table1_puri_ploi.csv")[,-1]
ind_pp = which(puri_ploi$SampleId %in% names(discordant[which(discordant != 0)]))
puri_ploi[ind_pp,]

# Manifest
luad_manifest_annot <- readRDS("~/data2/PureCN_manuscript/luad/luad_manifest_annot.rds")
ind_m = which(stringr::str_extract(luad_manifest_annot$barcode, "TCGA.{11}") %in% names(discordant[which(discordant != 0)]))
luad_manifest_annot[ind_m,]




# luad_LOH_all = readRDS("~/data2/PureCN_manuscript/Revision/LOH/data/luad_LOH_all.rds")
# names(luad_LOH_all) = stringr::str_extract(names(luad_LOH_all), "TCGA.{11}")
# ind = which(names(luad_LOH_all) %in% discordant_sampleId)
# luad_LOH_discordant = luad_LOH_all[ind]
