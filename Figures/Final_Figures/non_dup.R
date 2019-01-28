# 236 OVC samples
ovc_236 = read.csv("~/wallabe4_backup/github/PureCN_manuscript_backup/Results/duplication/ovc_236.csv")[,2]

# 442 LUAD samples
luad_442 = read.csv("~/wallabe4_backup/github/PureCN_manuscript_backup/luad/Results/duplication/luad_442.csv")[,2]
luad_442 = stringr::str_extract(luad_442, "TCGA.{11}")
