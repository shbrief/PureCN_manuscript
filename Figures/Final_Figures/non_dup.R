# 236 OVC samples
ovc_236 = read.csv("/data2/PureCN_manuscript/Results/duplication/ovc_236.csv")[,2]

# 442 LUAD samples
luad_442 = read.csv("/data2/PureCN_manuscript/luad/Results/duplication/luad_442.csv")[,2]
luad_442 = stringr::str_extract(luad_442, "TCGA.{11}")
