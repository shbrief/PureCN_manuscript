# Import purity and ploidy calls from PureCN outputs

output.dir = file.path("/data/16tb/CNVworkflow", kit, "purecn_output", 
                       paste0(kit, "_PureCN"), purecn_mode)

samples = list.files(output.dir)
optima_csv = file.path(output.dir, samples, paste0(samples, ".csv"))
optima_csv = optima_csv[file.exists(optima_csv)]

res = sapply(optima_csv[seq_along(optima_csv)], puri_ploi, USE.NAMES = FALSE) %>% 
  cbind %>% t %>% as.data.frame
res$capture_kit = kit
res_all = res   # save the additional columns

res = res[, c("submitter_id", "Purity", "Ploidy")]
res$submitter_id = as.character(res$submitter_id)
res$Purity = as.numeric(res$Purity)
res$Ploidy = as.numeric(res$Ploidy)
names(res) = c("SampleId", paste0("Purity_", purecn_mode), paste0("Ploidy_", purecn_mode))