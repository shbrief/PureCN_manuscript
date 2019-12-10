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