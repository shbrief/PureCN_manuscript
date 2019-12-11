

data_dir = "/nobackup/16tb_a/CNVworkflow_LUAD/purecn_output/PureCN"
tumor_dir = file.path(data_dir, "tumor_only")

sample_dir = list.files(tumor_dir)
file_list = sapply(sample_dir, function(x) list.files(file.path(tumor_dir, x)))
mut = sapply(file_list, function(x) x[grep("*.rds$", x)])

source("~/Documents/github/PureCN_manuscript/Figures/Final_Figures/non_dup.R")
sub_ind = which(stringr::str_extract(names(mut), "TCGA.{11}") %in% luad_442)
mut = mut[sub_ind]


file_path = file.path(tumor_dir, names(mut[1]), mut[[1]])
res = readRDS(file_path)

allcalls_mutationBurden = lapply(seq_along(mut), function(i) {
    pre_fname <- bed_subset[i, "pre_fname"]
    mycall_no_none <- paste("nice -n 20 java -Xmx24g -jar",
                            gatk.jar,
                            "--analysis_type SelectVariants",
                            "-R", REFERENCE,
                            "-V", file.path(out.dir, paste0(pre_fname, "_pon.vcf")),
                            "-o", file.path(out.dir, paste0(pre_fname, "_pon_no_none.vcf")))
})

allcalls_no_none <- allcalls_no_none[!file.exists(file.path(out.dir, paste0(pre_fname, "_pon_no_none.vcf")))]

res <- bplapply(allcalls_no_none || true, system, BPPARAM = MulticoreParam(24))