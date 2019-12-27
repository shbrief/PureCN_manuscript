#!/usr/bin/env Rscript
wd = "~/data2/PureCN_manuscript"

# All TCGA-LUAD samples were processed by cat# 931070
catalog = "931070"
source(file.path(wd, "Methods/input_info.R"))

# Load manifest file
luad_dir = (file.path(wd, "luad"))
manifest = readRDS(file.path(luad_dir, "luad_manifest_annot.rds"))

# Downloaded bam files which have corresponding ABSOLUTE results
luad_manifest_download = readr::read_tsv("/nobackup/16tb_a/luad_manifest_download_1.tsv")
manifest = manifest[manifest$id %in% luad_manifest_download$id,]
manifest$file_dir = "/data/ovarian/luad"
rm(luad_manifest_download)

# Separate manifest into normal(n) and tumor(t) samples
manifest_n = manifest[manifest$sample_definition %in% 
                        c("Blood Derived Normal", "Solid Tissue Normal"),]
manifest_t = manifest[manifest$sample_definition %in% 
                        c("Primary Solid Tumor", "Recurrent Solid Tumor"),]

# subset manifest file with the same capture kit/ runtype
if (exists("run_type")) {source("~/Documents/github/PureCN_manuscript/Methods/run_type.R")}
