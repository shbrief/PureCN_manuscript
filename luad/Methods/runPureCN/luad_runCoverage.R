# generate GC-normalized coverage files and build normalDB

# input arguments ----------------------------------------------------------
# run_type = "normal"
run_type = "tumor"
# --------------------------------------------------------------------------

script_dir = "/home/sehyun/Documents/github/PureCN_manuscript/luad/Methods"
source(file.path(script_dir, "luad_sample_sorting.R"))

PURECN = system.file("extdata", package = "PureCN")
BEDFILES = "/home/sehyun/reference/bedfiles"
OUT = "/data/ovarian/CNVworkflow_LUAD/purecn_output"
output = file.path(OUT, paste0(run_type, "_cov"))

# run Coverage.R on normal samples
allcalls_coverage = lapply(seq_along(bed_subset$id), function(i){
  uuid = bed_subset[i, "id"]
  fname = bed_subset[i, "filename"]
  data.dir = bed_subset[i, "file_dir"]
  fullname = file.path(data.dir, uuid, fname)
  bed_pre = gsub(".bed$", "", bed)
  intervals.txt = file.path(BEDFILES, bed_pre, paste0(bed_pre, "_hg38_gcgene.txt"))
  Coverage.R = file.path(PURECN, "Coverage.R")
  mycall = paste("nice -n 20 Rscript", Coverage.R,
                  "--outdir", output,
                  "--bam", fullname,
                  "--intervals", intervals.txt,
                  "--cores", 24)
})

# remove any run that has output file already
pre_fname = tools::file_path_sans_ext(bed_subset[,"filename"])
allcalls_coverage = allcalls_coverage[!file.exists(file.path(output, paste0(pre_fname, "_coverage_loess.txt")))]

library(BiocParallel)
res = bplapply(allcalls_coverage, system, BPPARAM = MulticoreParam(16))