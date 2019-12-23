# Respond to Markus' email about the paper revision (11/08/19)

# File paths to the PureCN outputs in wallabe4
luad_dir = "/nobackup/16tb_a/CNVworkflow_LUAD/purecn_output/PureCN/tumor_only"
ov_931070_dir = "/nobackup/16tb_b/CNVworkflow/931070/purecn_output/931070_PureCN/tumor_only"
ov_S0293689_dir = "/nobackup/16tb_b/CNVworkflow/S0293689/purecn_output/S0293689_PureCN/tumor_only"

# Load manifest files for sampleid extraction
home_dir = "~/Documents/github/PureCN_manuscript"
mani_luad <- readRDS(file.path(home_dir, "luad/luad_manifest_annot.rds"))
mani_ov <- readRDS(file.path(home_dir, "Data/ovc_manifest/ovc_manifest_112618.rds"))

# Subset the manifest files with the samples used for paper
x = read.csv(file.path(home_dir, "Figures/Final_Tables/Table1_puri_ploi.csv"))[,-1]$SampleId   # 442 LUAD + 233 OV
luad = which(stringr::str_extract(mani_luad$barcode, "TCGA.{11}") %in% as.character(x))
mani_luad = mani_luad[luad,]   # 442 samples
ov = which(stringr::str_extract(mani_ov$barcode, "TCGA.{11}") %in% as.character(x))
mani_ov = mani_ov[ov,]   # 259 samples


##### Write csv file ########################################################### 
# LUAD
path_to_dir = luad_dir
sampleid = gsub(".bam$", "", mani_luad$filename)
write.csv(do.call(rbind, lapply(file.path(path_to_dir, sampleid, paste0(sampleid,".csv")), 
                                read.csv, as.is=TRUE)), file = "luad.csv")

# OV_931070
path_to_dir = ov_931070_dir
sampleid = gsub(".bam$", "", mani_ov[which(mani_ov$target_capture_kit_catalog_number == "931070"),]$filename)
write.csv(do.call(rbind, lapply(file.path(path_to_dir, sampleid, paste0(sampleid,".csv")), 
                                           read.csv, as.is=TRUE)), file = "ov_931070.csv")

# OV_S0293689
path_to_dir = ov_S0293689_dir
sampleid = gsub(".bam$", "", mani_ov[which(mani_ov$target_capture_kit_catalog_number == "S0293689"),]$filename)
write.csv(do.call(rbind, lapply(file.path(path_to_dir, sampleid, paste0(sampleid,".csv")), 
                                read.csv, as.is=TRUE)), file = "ov_S0293689.csv")


##### Subset discordant samples ################################################
# List of "low purity or of high heterogeneity" samples, defined as purify difference
# more than 0.2 and ploidy difference more than 0.5
luad = read.table("Revision/Curation/luad.csv", sep = ",", header = TRUE)[,-1] # 442 LUAD
puri_ploi = read.csv(file.path(home_dir, "Figures/Final_Tables/Table1_puri_ploi.csv"))[,-1] # 675 samples = 233 OV + 442 LUAD
puri_ploi = puri_ploi[which(puri_ploi$Sample == "LUAD"),]

ploi_diff = abs(puri_ploi$Ploidy_ABS - puri_ploi$Ploidy_PCN_t)
puri_diff = abs(puri_ploi$Purity_ABS - puri_ploi$Purity_PCN_t)
discordant_ind = which(ploi_diff > 0.5 | puri_diff > 0.2)   # 116 discordant LUAD samples

discordant_sampleId = puri_ploi$SampleId[discordant_ind] # select the SampleId of discordant samples
luad$Sampleid = stringr::str_extract(luad$Sampleid, "TCGA.{11}") %>% gsub("\\.", "-",.)  # change SampleId format
luad_discordant = luad[which(luad$Sampleid %in% discordant_sampleId),]
write.csv(luad_discordant, "~/data2/PureCN_manuscript/Revision/Curation/luad_discordant.csv")
