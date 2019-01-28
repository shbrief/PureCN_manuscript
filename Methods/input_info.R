#!/usr/bin/env Rscript

# --catalog argument
if (catalog == "931070") {
  catalog = "931070"
  kit = "Custom V2 Exome Bait, 48 RXN X 16 tubes"
  bed = "whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed"
} else if (catalog == "S0293689") {
  catalog = "S0293689"
  kit = "SureSelect Human All Exon 38 Mb v2"
  bed = "SeqCap_EZ_Exome_v3_capture.bed"
} 