#!/usr/bin/env Rscript
# liftOver BED file from hg19 to hg38

bedfiles.dir = "/home/sehyun/reference/bedfiles"
data.dir = "/shared/ovc_purecn/data"
hg38.dir = file.path(data.dir, "hg19ToHg38")
BED = "cancer_2000gene_shift170.targetIntervals"
  
# Change chromosome notation from '1, 2, 3, ...' to 'chr1, chr2, chr3, ...
org_bed = file.path(data.dir, paste0(BED, ".bed"))
new_bed = file.path(hg38.dir, paste0(BED, "_new.bed"))
system(paste("cat", org_bed, "| sed 's/^/chr/' >", new_bed))

## Check Meetup RMD file for updated code
library(rtracklayer)
ch = import.chain(file.path(bedfiles.dir, "hg19ToHg38.over.chain"))
bed_df = read.table(new_bed)
bed_df_ordered = bed_df[order(bed_df$V1, bed_df$V2),]
bed_gr = GRanges(seqnames = bed_df_ordered$V1, 
                 ranges = IRanges(start = bed_df_ordered$V2, end = bed_df_ordered$V3), 
                 strand = bed_df_ordered$V5)
bed_hg38 = liftOver(bed_gr, ch)
unlist_bed_hg38 = unlist(bed_hg38)

dir.create(file.path(bedfiles.dir, BED))
new_dir = file.path(bedfiles.dir, BED)
export.bed(unlist_bed_hg38, file.path(new_dir, paste0(BED,"_hg38.bed")))
