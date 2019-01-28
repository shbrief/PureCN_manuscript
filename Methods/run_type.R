# subset manifest file with the same capture kit/ runtype
if (run_type == "normal") {
  kit_subset_n = manifest_n[which(manifest_n$target_capture_kit_name == kit),]
  bed_subset = kit_subset_n[which(kit_subset_n$bedfiles == bed),]
} else if (run_type == "tumor") {
  kit_subset_t = manifest_t[which(manifest_t$target_capture_kit_name == kit),]
  bed_subset = kit_subset_t[which(kit_subset_t$bedfiles == bed),]
}