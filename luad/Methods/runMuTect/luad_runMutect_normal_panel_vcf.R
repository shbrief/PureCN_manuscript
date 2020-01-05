## Generating the normal panel VCF (normal.panel.vcf.file) for mapping bias correction
## This script is based on Levi's "/shared/ovc_purecn/src/runMutect.R"

# input arguments ----------------------------------------------------------
run_type = "normal"
script_dir = "~/data2/PureCN_manuscript/luad/Methods/"
source(file.path(script_dir, "luad_sample_sorting.R"))
# --------------------------------------------------------------------------

mutect.dir = "/home/lwaldron/mutect"
ref.dir = "/home/sehyun/reference"
MUTECT_OUT = file.path("/nobackup/16tb_a/CNVworkflow_LUAD/mutect_output")
out.dir = file.path(MUTECT_OUT, "normal_panel")

java7.exec = "/usr/local/share/jre1.7.0_79/bin/java"
mutect.jar = "/usr/local/share/mutect/1.1.7/bin/mutect-1.1.7.jar"
gatk.jar = "/usr/local/share/gatk/3.6/GenomeAnalysisTK.jar"

REFERENCE = file.path(ref.dir, "GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
DBSNP_VCF = file.path(mutect.dir, "hg38bundle/Homo_sapiens_assembly38.dbsnp.vcf.gz")
COSMIC_VCF = file.path(mutect.dir, "CosmicCodingMuts.vcf.gz")

bed_subset$pre_fname = tools::file_path_sans_ext(bed_subset[, "filename"])

# run MuTect using normal samples under '--artifact_detection_mode'
allcalls_pon <- lapply(seq_along(bed_subset$id), function(i){
  uuid <- bed_subset[i, "id"]
  fname <- bed_subset[i, "filename"]
  data.dir <- bed_subset[i, "file_dir"]
  fullname <- file.path(data.dir, uuid, fname)
  submitter_id <- bed_subset[i, "submitter_id"]
  pre_fname <- bed_subset[i, "pre_fname"]
  mycall_pon <- paste("nice -n 20 ",
                      java7.exec,
                      "-Xmx6g -jar",
                      mutect.jar,
                      "--analysis_type MuTect",
                      "-R", REFERENCE,
                      "--artifact_detection_mode",
                      "--dbsnp", DBSNP_VCF,
                      "--cosmic", COSMIC_VCF,
                      "-dt None",
                      "-I:tumor", fullname,
                      "-o", file.path(out.dir, paste0(pre_fname, "_pon_stats.txt")),
                      "-vcf", file.path(out.dir, paste0(pre_fname, "_pon.vcf")))
})

pre_fname <- bed_subset$pre_fname
allcalls_pon <- allcalls_pon[!file.exists(file.path(out.dir, paste0(pre_fname, "_pon_stats.txt")))]

library(BiocParallel)
res <- bplapply(allcalls_pon, system, BPPARAM = MulticoreParam(24))

# Remove the empty none sample from the VCF
# use VCF in $OUT/${ID}_no_none.vcf for CombineVariants
allcalls_no_none <- lapply(seq_along(bed_subset$id), function(i){
  pre_fname <- bed_subset[i, "pre_fname"]
  mycall_no_none <- paste("nice -n 20 java -Xmx24g -jar",
                          gatk.jar,
                          "--analysis_type SelectVariants",
                          "-R", REFERENCE,
                          "-V", file.path(out.dir, paste0(pre_fname, "_pon.vcf")),
                          "-o", file.path(out.dir, paste0(pre_fname, "_pon_no_none.vcf")))
})

allcalls_no_none <- allcalls_no_none[!file.exists(file.path(out.dir, paste0(pre_fname, "_pon_no_none.vcf")))]

res <- bplapply(allcalls_no_none, system, BPPARAM = MulticoreParam(24))

## Merge normal VCFs into one
pon_size = 250
merged.file = file.path(MUTECT_OUT, paste0(catalog, ".normals.merged.min5_", pon_size, ".vcf"))
merge.cmd = paste("java -Xmx24g -jar", gatk.jar, "-T CombineVariants -nt 4 --minimumN 5 --genotypemergeoption UNSORTED -R", REFERENCE)
merge.cmd = paste(merge.cmd, paste("--variant", dir(out.dir, pattern="^.*_pon_no_none\\.vcf$", full.names = TRUE)[1: pon_size], collapse=" "))
merge.cmd = paste(merge.cmd, ">", merged.file)
system(merge.cmd)


