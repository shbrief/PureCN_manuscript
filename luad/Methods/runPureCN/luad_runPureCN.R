#!/usr/bin/env Rscript

# input arguments ----------------------------------------------------------
run_type = "tumor"
# run_type = "normal"

luad_dir = "/home/sehyun/Documents/github/PureCN_manuscript/luad"
source(file.path(luad_dir, "Methods/luad_sample_sorting.R"))

# run_type = "tumor_only"
run_type = "matching_normal"
# --------------------------------------------------------------------------

PURECN=system.file("extdata", package = "PureCN")
OUT="/data/ovarian/CNVworkflow_LUAD"
MUTECT_OUT=file.path(OUT, "mutect_output")
PURECN_OUT=file.path(OUT, "purecn_output")
REF="/home/sehyun/reference"
BEDFILES=file.path(REF, "bedfiles")
BLACKLIST="/home/sehyun/Documents/github/PureCN_manuscript/Data/snp_blacklist/hg38_simpleRepeats.bed"

bed_subset$pre_fname = tools::file_path_sans_ext(bed_subset[,"filename"])  

# Lung cancer is much more heterogeneous than ovarian cancer. 
# Default of —maxnonclonal is 0.2, so I would try if 0.3 works better.
# I didn't put 'seed' before Nov.23 2018 runs
allcalls_purecn <- lapply(seq_along(bed_subset$id), function(i){
  PureCN.R <- file.path(PURECN, "PureCN.R")
  pre_fname <- bed_subset[i, "pre_fname"]
  bed_pre <- gsub(".bed$", "", bed)
  dir.create(file.path(PURECN_OUT, "PureCN", run_type, pre_fname))
  out_dir = file.path(PURECN_OUT, "PureCN", run_type, pre_fname)
  
  mycall_purecn <- paste("nice -n 20 Rscript", PureCN.R, 
                         "--out", out_dir,
                         "--tumor", file.path(PURECN_OUT, "tumor_cov", paste0(pre_fname, "_coverage_loess.txt")), 
                         "--sampleid", pre_fname,
                         "--vcf", file.path(MUTECT_OUT, paste0("stat_", run_type), paste0(pre_fname, "_matching_mutect.vcf")),
                         "--statsfile", file.path(MUTECT_OUT, paste0("stat_", run_type), paste0(pre_fname, "_matching_mutect_stat.txt")),
                         "--normaldb", file.path(PURECN_OUT, "normalDB/normalDB_hg38.rds"),
                         "--normal_panel", file.path(PURECN_OUT, "normalDB/mapping_bias_hg38.rds"),
                         "--intervals", file.path(BEDFILES, bed_pre, paste0(bed_pre, "_hg38_gcgene.txt")),
                         "--targetweightfile", file.path(PURECN_OUT, "normalDB/interval_weights_hg38.txt"),
                         "--snpblacklist", BLACKLIST,
                         "--genome hg38 --postoptimize --seed 123 --maxnonclonal 0.3")
})

# remove any run that has output file already
pre_fname <- tools::file_path_sans_ext(bed_subset[,"filename"])
allcalls_purecn <- allcalls_purecn[!file.exists(file.path(PURECN_OUT, "PureCN", run_type, pre_fname, paste0(pre_fname, "_chromosomes.pdf")))]

library(BiocParallel)
res <- bplapply(allcalls_purecn, system, BPPARAM = MulticoreParam(24))
