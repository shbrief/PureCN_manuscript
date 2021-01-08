#!/usr/bin/env Rscript

# input arguments ----------------------------------------------------------
run_type = "tumor_only"

catalog = "931070"
# catalog = "S0293689"

methods_dir = "/home/sehyun/Documents/github/PureCN_manuscript/Methods"
source(file.path(methods_dir, "input_info.R"))
source(file.path(methods_dir, "sample_sorting.R"))
# --------------------------------------------------------------------------

# Paths to the files
# PURECN = system.file("extdata", package = "PureCN")
data_dir = "/home/sehyun/Documents/github/PureCN_manuscript/Data"
mutect_dir = file.path("/data/16tb/CNVworkflow", catalog, "mutect_output/stat_matching_normal")

out_dir = file.path("/nobackup/16tb_b/CNVworkflow", catalog, "purecn_output")
rds_dir = file.path(out_dir, paste0(catalog, "_PureCN"), run_type)
rds_list = list.files(rds_dir)

# Create/Set the directory for Benchmark.R outputs
Benchmark_out_dir = file.path(out_dir, paste0(catalog, "_Benchmark"))
if (!file.exists(Benchmark_out_dir)) {dir.create(Benchmark_out_dir)}

Benchmark_out_dir = file.path(out_dir, paste0(catalog, "_Benchmark/tumor_only_SGZ"))
if (!file.exists(Benchmark_out_dir)) {dir.create(Benchmark_out_dir)}

# Run Benchmark.R
allcalls_Benchmark <- lapply(seq_along(rds_list), function(i){
  # Benchmark.R = file.path(PURECN, "BenchmarkTumorOnly.R")
  Benchmark.R = "/home/sehyun/Documents/github/SGZ/BenchmarkTumorOnlySGZ.R"
  mycall_Benchmark <- paste("Rscript", Benchmark.R,
                     "--out", file.path(Benchmark_out_dir, rds_list[i]),
                     "--rds", file.path(rds_dir, rds_list[i], paste0(rds_list[i], ".rds")),
                     "--vcf", file.path(mutect_dir, paste0(rds_list[i], "_matching_mutect.vcf")),
                     "--sgzdir /home/sehyun/Documents/github/SGZ/SGZ"
                     )
})

out_fname = paste0(rds_list, "_sgz.tumor_only_benchmark.csv")
allcalls_Benchmark <- allcalls_Benchmark[!file.exists(file.path(Benchmark_out_dir, out_fname))]

# test the first 10 samples 
library(BiocParallel)
res <- bplapply(allcalls_Benchmark, system, BPPARAM = MulticoreParam(20))