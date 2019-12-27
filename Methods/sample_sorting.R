#!/usr/bin/env Rscript

# Load manifest file
manifest = readRDS("~/data2/PureCN_manuscript/Data/manifest/ovc_manifest_112618.rds")

# Subset manifest with only one bedfile
bedfiles = manifest$bedfiles
summary = data.frame()

for (i in 1:length(bedfiles)){
  summary = rbind(summary, c(i, length(bedfiles[[i]])))
  colnames(summary) = c("i", "num_bedfiles")
}

manifest = manifest[which(summary$num_bedfiles == 1),]
  
# Separate manifest into normal(n) and tumor(t) samples
manifest_n = manifest[manifest$sample_definition %in% 
                        c("Blood Derived Normal", "Solid Tissue Normal"),]
manifest_t = manifest[manifest$sample_definition %in% 
                        c("Primary Solid Tumor", "Recurrent Solid Tumor"),]

if (exists("run_type")) {source("~/data2/PureCN_manuscript/Methods/run_type.R")}
