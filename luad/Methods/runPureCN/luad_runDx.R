#!/usr/bin/env Rscript

# input arguments ----------------------------------------------------------
catalog = "931070"
run_type = "matching_normal"    # c("tumor_only", "matching_normal")
# --------------------------------------------------------------------------

# Paths to the files
PURECN = system.file("extdata", package = "PureCN")
data_dir = "/home/sehyun/Documents/github/PureCN_manuscript/Data"
BLACKLIST = file.path(data_dir, "snp_blacklist/hg38_simpleRepeats.bed")

library(dplyr)
out_dir = file.path("/data/ovarian/CNVworkflow_LUAD/purecn_output")
rds_dir = file.path(out_dir, "PureCN", run_type)
rds_list = list.files(rds_dir) %>% as.data.frame()

# Create/Set the directory for Dx.R outputs
Dx_out_dir = file.path(out_dir, "Dx")
if (!file.exists(Dx_out_dir)) {dir.create(Dx_out_dir)}

Dx_out_dir = file.path(out_dir, "Dx", paste0(run_type, "_cds"))
if (!file.exists(Dx_out_dir)) {dir.create(Dx_out_dir)}

# Convert file names to a more readable version
luad_dir = "~/Documents/github/PureCN_manuscript/luad"
manifest = readRDS(file.path(luad_dir, "luad_manifest_annot.rds"))
manifest$fname = tools::file_path_sans_ext(manifest[,"filename"])

for (i in 1:nrow(rds_list)) {
  colnames(rds_list)[1] = "fname"
  rds_list$barcode[i] = manifest[which(manifest$fname == rds_list[i,1]),]$barcode
}

# Run Dx.R
allcalls_Dx <- lapply(1:nrow(rds_list), function(i){
  Dx.R = file.path(PURECN, "Dx.R")
  dir.create(file.path(Dx_out_dir, rds_list$barcode[i]))
  mycall_Dx <- paste("Rscript", Dx.R,
                     "--out", file.path(Dx_out_dir, rds_list$barcode[i], rds_list$barcode[i]),
                     "--rds", file.path(rds_dir, rds_list$fname[i], paste0(rds_list$fname[i], ".rds")),
                     "--callable", file.path(out_dir, "Dx/callableLoci/filtered_cds", paste0(rds_list$fname[i], "_callable_status_filtered_cds.bed")),
                     "--exclude", BLACKLIST,
                     "--signatures")
})

allcalls_Dx <- allcalls_Dx[!file.exists(file.path(Dx_out_dir, rds_list$barcode, paste0(rds_list$barcode, "_signatures.pdf")))]

# test the first 10 samples 
library(BiocParallel)
res <- bplapply(allcalls_Dx, system, BPPARAM = MulticoreParam(20))