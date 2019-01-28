# script to select files to download
# 'orig_mani' and 'ovcnorm_mani' seem like a same file
orig_mani = read.table("/data/ovarian/ovcnormals/orig.ovcnormals_manifest.tsv", header=TRUE)
ovcnorm_mani = read.table("/data/ovarian/ovcnormals/ovcnormals_manifest.tsv", header=TRUE)
gdc_mani = read.table("/data/16tb/ovc_data/gdc_manifest_20171213_185003.txt", header=TRUE)

# subset 'gdc_mani' with the files that haven't downloaded yet. 
gdc_mani_download = gdc_mani[!(gdc_mani$md5 %in% ovcnorm_mani$md5),]
write.table(gdc_mani_download, file="/data/16tb/ovc_data/gdc_mani_download.tsv", quote=FALSE, sep="\t", row.names=FALSE)

# Shell script for downloading bam files
gdc-client download -m /data/16tb/ovc_data/gdc_mani_download.tsv -t /data/16tb/ovc_data/gdc-user-token.2017-12-13T17_15_08.632Z.txt
