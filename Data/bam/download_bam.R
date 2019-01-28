# check '/home/sehyun/reference/ovc_manifest/gdc_client_download.R'
# script to select files to download
ovcnorm_mani = read.table("/data/ovarian/ovcnormals/ovcnormals_manifest.tsv", header=TRUE)
gdc_mani = read.table("/data/16tb/ovc_data/gdc_manifest_20171213_185003.txt", header=TRUE)

# subset 'gdc_mani' with the files that haven't downloaded yet
gdc_mani_download = gdc_mani[!(gdc_mani$md5 %in% ovcnorm_mani$md5),]

# save new manifest files for bam files to be downloaded
data.dir = "/data/16tb/ovc_data"
write.table(gdc_mani_download, file = file.path(data.dir, "gdc_mani_download.tsv"), quote=FALSE, sep="\t", row.names=FALSE)

