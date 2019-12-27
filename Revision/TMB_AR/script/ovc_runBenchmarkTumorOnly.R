#!/usr/bin/env Rscript

# input arguments ----------------------------------------------------------
# run_type = "tumor_only"
run_type = "matching_normal"

catalog = "931070"
# catalog = "S0293689"

methods_dir = "~/data2/PureCN_manuscript/Methods"
source(file.path(methods_dir, "input_info.R"))
source(file.path(methods_dir, "sample_sorting.R"))
# --------------------------------------------------------------------------

# Paths to the files
PURECN = system.file("extdata", package = "PureCN")
data_dir = "~/data2/PureCN_manuscript/Data"
mutect_dir = file.path("/nobackup/16tb_b/CNVworkflow", catalog, "mutect_output/stat_matching_normal")

out_dir = file.path("/nobackup/16tb_b/CNVworkflow", catalog, "purecn_output")
rds_dir = file.path(out_dir, paste0(catalog, "_PureCN"), run_type)
rds_list = list.files(rds_dir)

# Create/Set the directory for Benchmark.R outputs
# For revision, I created BenchmarkTumorOnly outputs' directory
Benchmark_out_dir = file.path(out_dir, paste0(catalog, "_BenchmarkTumorOnly"))
if (!file.exists(Benchmark_out_dir)) {dir.create(Benchmark_out_dir)}

Benchmark_out_dir = file.path(out_dir, paste0(catalog, "_BenchmarkTumorOnly"), run_type)
if (!file.exists(Benchmark_out_dir)) {dir.create(Benchmark_out_dir)}

# Run Benchmark.R
allcalls_Benchmark <- lapply(seq_along(rds_list), function(i){
    Benchmark.R = "~/data2/PureCN_manuscript/Revision/TMB_AR/script/BenchmarkTumorOnly.R"
    mycall_Benchmark <- paste("Rscript", Benchmark.R,
                     "--out", file.path(Benchmark_out_dir, rds_list[i]),
                     "--rds", file.path(rds_dir, rds_list[i], paste0(rds_list[i], ".rds")),
                     "--vcf", file.path(mutect_dir, paste0(rds_list[i], "_matching_mutect.vcf")),
                     "--callable", file.path(out_dir, paste0(catalog, "_Dx/callableLoci/filtered_cds/", rds_list[i], "_callable_status_filtered_cds.bed"))
                     )
})

out_fname = paste0(rds_list, "_tumor_only_benchmark.csv")
allcalls_Benchmark <- allcalls_Benchmark[!file.exists(file.path(Benchmark_out_dir, out_fname))]

# test the first 10 samples 
library(BiocParallel)
res <- bplapply(allcalls_Benchmark, system, BPPARAM = MulticoreParam(20))