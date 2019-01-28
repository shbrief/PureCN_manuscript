#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))

option_list = list(
  make_option(c("--catalog"), action="store", type="character", default=NULL,
              help = "Capture Kit Catalog Number"),
  make_option(c("--runtype"), action="store", type="character", default=NULL, 
              help = "'tumor_only' or 'matching_normal'")
)

opt = parse_args(OptionParser(option_list=option_list))
run_type = opt$runtype

# --catalog argument
if (opt$catalog == "931070") {
  catalog = "931070"
  kit = "Custom V2 Exome Bait, 48 RXN X 16 tubes"
  bed = "whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed"
} else if (opt$catalog == "S0293689") {
  catalog = "S0293689"
  kit = "SureSelect Human All Exon 38 Mb v2"
  bed = "SeqCap_EZ_Exome_v3_capture.bed"
} 