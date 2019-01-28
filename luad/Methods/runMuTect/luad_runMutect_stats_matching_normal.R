## Generating vcf.file and stats.txt file for artifact filtering

# input arguments ----------------------------------------------------------
run_type = "normal"

# All TCGA-LUAD samples were processed by cat# 931070
catalog = "931070"
source("~/Documents/github/PureCN_manuscript/Methods/input_info.R")

# Load manifest file
luad_dir = "~/Documents/github/PureCN_manuscript/luad"
manifest = readRDS(file.path(luad_dir, "luad_manifest_annot.rds"))

# Downloaded bam files which have corresponding ABSOLUTE results
luad_manifest_download = readr::read_tsv("/data/ovarian/luad_manifest_download_1.tsv")
manifest = manifest[manifest$id %in% luad_manifest_download$id,]
manifest$file_dir = "/data/ovarian/luad"
rm(luad_manifest_download)

# Separate manifest into normal(n) and tumor(t) samples
manifest_n = manifest[manifest$sample_definition %in% 
                        c("Blood Derived Normal", "Solid Tissue Normal"),]
manifest_t = manifest[manifest$sample_definition %in% 
                        c("Primary Solid Tumor", "Recurrent Solid Tumor"),]
# --------------------------------------------------------------------------

mutect.dir = "/home/lwaldron/mutect"
ref.dir = "/home/sehyun/reference"
MUTECT_OUT = file.path("/data/ovarian/CNVworkflow_LUAD/mutect_output")
out.dir = file.path(MUTECT_OUT, "stat_matching_normal")

java7.exec = "/usr/local/share/jre1.7.0_79/bin/java"
mutect.jar = "/usr/local/share/mutect/1.1.7/bin/mutect-1.1.7.jar"

REFERENCE = file.path(ref.dir, "GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
DBSNP_VCF = file.path(mutect.dir, "hg38bundle/Homo_sapiens_assembly38.dbsnp.vcf.gz")
COSMIC_VCF = file.path(mutect.dir, "CosmicCodingMuts.vcf.gz")

kit_subset_t = manifest_t[which(manifest_t$target_capture_kit_name == kit),]
bed_subset_t = kit_subset_t[which(kit_subset_t$bedfiles == bed),]

kit_subset_n = manifest_n[which(manifest_n$target_capture_kit_name == kit),]
bed_subset_n = kit_subset_n[which(kit_subset_n$bedfiles == bed),]

# process samples with only one tumor/normal pair ------------------------------------
library(dplyr)
dup_t = bed_subset_t %>% group_by(submitter_id) %>% filter(n()>1)   # 89 obs.
bed_subset_t = bed_subset_t[which(!bed_subset_t$barcode %in% dup_t$barcode),]   # 479 obs.
dup_n = bed_subset_n %>% group_by(submitter_id) %>% filter(n()>1)   # 144 obs.
bed_subset_n = bed_subset_n[which(!bed_subset_n$barcode %in% dup_n$barcode),]   # 450 obs.
# ------------------------------------------------------------------------------------

matching = which(bed_subset_t$submitter_id %in% bed_subset_n$submitter_id)
bed_subset = bed_subset_t[matching,]

allcalls_stat <- lapply(seq_along(bed_subset$id), function(i){
  uuid <- bed_subset[i, "id"]
  fname <- bed_subset[i, "filename"]
  data.dir <- bed_subset[i, "file_dir"]
  fullname <- file.path(data.dir, uuid, fname)
  submitter_id <- bed_subset[i, "submitter_id"]
  
  matching_i <- which(bed_subset_n$submitter_id == submitter_id)
  matching_uuid <- bed_subset_n[matching_i, "id"]
  matching_fname <- bed_subset_n[matching_i, "filename"]
  matching_data.dir <- bed_subset_n[matching_i, "file_dir"]
  matching_fullname <- file.path(matching_data.dir, matching_uuid, matching_fname)
  
  pre_fname <- tools::file_path_sans_ext(bed_subset[i,"filename"])
  
  mycall_stat <- paste("nice -n 20 ",
                       java7.exec,
                       "-Xmx6g -jar",
                       mutect.jar,
                       "--analysis_type MuTect",
                       "-R", REFERENCE,
                       "--dbsnp", DBSNP_VCF,
                       "--cosmic", COSMIC_VCF,
                       "-I:tumor", fullname,
                       "-I:normal", matching_fullname,
                       "-o", file.path(out.dir, paste0(pre_fname, "_matching_mutect_stat.txt")),
                       "-vcf", file.path(out.dir, paste0(pre_fname, "_matching_mutect.vcf")))
})

pre_fname <- tools::file_path_sans_ext(bed_subset[,"filename"])
allcalls_stat <- allcalls_stat[!file.exists(file.path(out.dir, paste0(pre_fname, "_matching_mutect.vcf")))]

library(BiocParallel)
res <- bplapply(allcalls_stat, system, BPPARAM = MulticoreParam(16))
