# Extract TMB 
# Dx.R was run with '--callable' argument using '{SAMPLE}_callable_status_filtered_cds.bed'

catalog = "931070"

TMB_dir = "~/Documents/github/PureCN_manuscript/luad/Results/TMB"
out_dir = file.path(TMB_dir, "extdata")

data_dir = "/nobackup/16tb_a/CNVworkflow_LUAD/purecn_output/Dx"
tumor_dir = file.path(data_dir, "tumor_only_cds")
paired_dir = file.path(data_dir, "matching_normal_cds")

# a list of SampleID, which has '_mutation_burden.csv' output
mut = list.files(paired_dir)
mut = mut[file.exists(file.path(paired_dir, mut, paste0(mut, "_mutation_burden.csv")))]

# extract somatic.rate.ontarget
mut_list = lapply(seq_along(mut), function(i){
  mut_fname = file.path(mut[i], paste0(mut[i], "_mutation_burden.csv"))
  mut_tumor = read.csv(file.path(tumor_dir, mut_fname))
  mut_paired = read.csv(file.path(paired_dir, mut_fname))
  
  data.frame(sampleID = mut[i], 
            TMB_tumor = mut_tumor$somatic.rate.ontarget, 
            TMB_paired = mut_paired$somatic.rate.ontarget)
}) 

mut_rate = data.frame()
for (i in seq_along(mut_list)) {mut_rate = rbind(mut_rate, mut_list[[i]])}

# save a table of TMB outputs
write.table(mut_rate, file = file.path(out_dir, paste0(catalog, "_TMB.tsv")))