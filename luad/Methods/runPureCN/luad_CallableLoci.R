#!/usr/bin/env Rscript

# input arguments ----------------------------------------------------------
catalog = "931070"
run_type = "tumor_only"    # c("tumor_only", "matching_normal")

methods_dir = "/home/sehyun/Documents/github/PureCN_manuscript/luad/Methods"
source(file.path(methods_dir, "luad_sample_sorting.R"))
# --------------------------------------------------------------------------

java7.exec = "/usr/local/share/jre1.7.0_79/bin/java"
gatk.jar = "/usr/local/share/gatk/3.6/GenomeAnalysisTK.jar"

PURECN = system.file("extdata", package = "PureCN")

ref.dir = "/home/sehyun/reference"
REFERENCE = file.path(ref.dir, "GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
Dx_OUT = "/data/ovarian/CNVworkflow_LUAD/purecn_output/Dx/callableLoci"

kit_subset = manifest_t[which(manifest_t$target_capture_kit_name == kit),]
bed_subset = kit_subset[which(kit_subset$bedfiles == bed),]

## allcalls_callableLoci
if (!file.exists(file.path(Dx_OUT, "callable_status"))){
  dir.create(file.path(Dx_OUT, "callable_status"))
}
out.dir_1 = file.path(Dx_OUT, "callable_status")

allcalls_callableLoci <- lapply(seq_along(bed_subset$id), function(i){
  uuid <- bed_subset[i, "id"]
  fname <- bed_subset[i, "filename"]
  data.dir <- bed_subset[i, "file_dir"]
  fullname <- file.path(data.dir, uuid, fname)
  pre_fname <- tools::file_path_sans_ext(bed_subset[i,"filename"])
  
  mycall_callableLoci <- paste("java -Xmx24g -jar", gatk.jar, 
                       "-T CallableLoci",
                       "-R", REFERENCE,
                       "-I", fullname,
                       "--summary", file.path(out.dir_1, paste0(pre_fname, "_table.txt")),
                       "-o", file.path(out.dir_1, paste0(pre_fname, "_callable_status.bed")),
                       "--minDepth 30")
})

pre_fname <- tools::file_path_sans_ext(bed_subset[,"filename"])
allcalls_callableLoci <- allcalls_callableLoci[!file.exists(file.path(out.dir_1, paste0(pre_fname, "_callable_status.bed")))]

# test the first 10 samples 
library(BiocParallel)
res <- bplapply(allcalls_callableLoci, system, BPPARAM = MulticoreParam(8))

## Filter only callable region for Dx.R run
if (!file.exists(file.path(Dx_OUT, "filtered"))){
  dir.create(file.path(Dx_OUT, "filtered"))
}
out.dir_2 = file.path(Dx_OUT, "filtered")

allcalls_filtered <- lapply(seq_along(bed_subset$id), function(i){
  pre_fname <- tools::file_path_sans_ext(bed_subset[i,"filename"])
  mycall_filtered <- paste("grep CALLABLE", 
                           file.path(out.dir_1, paste0(pre_fname, "_callable_status.bed")), ">",
                           file.path(out.dir_2, paste0(pre_fname, "_callable_status_filtered.bed")))
})

pre_fname <- tools::file_path_sans_ext(bed_subset[,"filename"])
allcalls_filtered <- allcalls_filtered[!file.exists(file.path(out.dir_2, paste0(pre_fname, "_callable_status_filtered.bed")))]

library(BiocParallel)
res <- bplapply(allcalls_filtered, system, BPPARAM = MulticoreParam(8))

## Restrict mutation burden calculation ot coding sequences
if (!file.exists(file.path(Dx_OUT, "filtered_cds"))){
  dir.create(file.path(Dx_OUT, "filtered_cds"))
}
out.dir_3 = file.path(Dx_OUT, "filtered_cds")

allcalls_cds <- lapply(seq_along(bed_subset$id), function(i){
  filter_cds.R <- file.path(PURECN, "FilterCallableLoci.R")
  pre_fname = tools::file_path_sans_ext(bed_subset[i,"filename"])
  
  mycall_cds <- paste("nice -n 20 Rscript", filter_cds.R, "--genome hg38",
                         "--infile", file.path(out.dir_2, paste0(pre_fname, "_callable_status_filtered.bed")), 
                         "--outfile", file.path(out.dir_3, paste0(pre_fname, "_callable_status_filtered_cds.bed")))
})

pre_fname <- tools::file_path_sans_ext(bed_subset[,"filename"])
allcalls_cds <- allcalls_cds[!file.exists(file.path(out.dir_3, paste0(pre_fname, "_callable_status_filtered_cds.bed")))]

library(BiocParallel)
res <- bplapply(allcalls_cds, system, BPPARAM = MulticoreParam(8))