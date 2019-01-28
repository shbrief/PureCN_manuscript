## Generating vcf.file and stats.txt file for artifact filtering

# input arguments ----------------------------------------------------------
run_type = "tumor"
script_dir = "/home/sehyun/Documents/github/PureCN_manuscript/luad/Methods/"
source(file.path(script_dir, "luad_sample_sorting.R"))
# --------------------------------------------------------------------------

mutect.dir = "/home/lwaldron/mutect"
ref.dir = "/home/sehyun/reference"
MUTECT_OUT = file.path("/data/ovarian/CNVworkflow_LUAD/mutect_output")
out.dir = file.path(MUTECT_OUT, "stat_tumor_only")

java7.exec = "/usr/local/share/jre1.7.0_79/bin/java"
mutect.jar = "/usr/local/share/mutect/1.1.7/bin/mutect-1.1.7.jar"

REFERENCE = file.path(ref.dir, "GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
DBSNP_VCF = file.path(mutect.dir, "hg38bundle/Homo_sapiens_assembly38.dbsnp.vcf.gz")
COSMIC_VCF = file.path(mutect.dir, "CosmicCodingMuts.vcf.gz")

kit_subset = manifest_t[which(manifest_t$target_capture_kit_name == kit),]
bed_subset = kit_subset[which(kit_subset$bedfiles == bed),]
bed_subset$pre_fname = tools::file_path_sans_ext(bed_subset[, "filename"])

allcalls_stat <- lapply(seq_along(bed_subset$id), function(i){
  uuid <- bed_subset[i, "id"]
  fname <- bed_subset[i, "filename"]
  data.dir <- bed_subset[i, "file_dir"]
  fullname <- file.path(data.dir, uuid, fname)
  submitter_id <- bed_subset[i, "submitter_id"]
  pre_fname <- bed_subset[i, "pre_fname"]
  
  mycall_stat <- paste("nice -n 20 ",
                  java7.exec,
                  "-Xmx6g -jar",
                  mutect.jar,
                  "--analysis_type MuTect",
                  "-R", REFERENCE,
                  "--dbsnp", DBSNP_VCF,
                  "--cosmic", COSMIC_VCF,
                  "-I:tumor", fullname,
                  "-o", file.path(out.dir, paste0(pre_fname, "_mutect_stat.txt")),
                  "-vcf", file.path(out.dir, paste0(pre_fname, "_mutect.vcf")))
})

output_fname <- file.path(out.dir, paste0(bed_subset$pre_fname, "_mutect.vcf"))
allcalls_stat <- allcalls_stat[!file.exists(output_fname)]


library(BiocParallel)
res <- bplapply(allcalls_stat, system, BPPARAM = MulticoreParam(16))
