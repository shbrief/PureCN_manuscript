# hg19toHg38 liftOver ABSOLUTE_OV_grangeslist
# system.time() for one range
#    user  system elapsed 
#   3.928   0.020   3.947
# total # of range = 76,062 --> This will take over 3.5 days?

library(rtracklayer)
library(IRanges)
library(GenomicRanges)

abs = readRDS("/home/sehyun/Documents/github/subtypeHeterogeneity/inst/extdata/ABSOLUTE_OV_grangeslist.rds")
ch = import.chain("~/reference/bedfiles/hg19ToHg38.over.chain")

abs_hg38 = abs

for (i in seq_along(abs)) {
  for (j in seq_along(abs[[i]])) {
    abs_ij_hg38 = liftOver(abs[[i]][j], ch)
    abs_ij_hg38 = reduce(abs_ij_hg38)
    
    gr = abs_ij_hg38
    gr = gr[seqnames(gr) == as.character(seqnames(abs[[i]][j]))]
    
    start_hg38 = min(start(ranges(gr)))
    end_hg38 = max(end(ranges(gr)))
    
    ranges(abs_hg38[[i]][j]) = IRanges(start_hg38, end_hg38)
  }
}

saveRDS(abs_hg38, file = "/home/sehyun/Documents/github/PureCN_manuscript/Results/ABSOLUTE/ABSOLUTE_OV_grangeslist_hg38")