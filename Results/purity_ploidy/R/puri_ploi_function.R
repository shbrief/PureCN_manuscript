puri_ploi = function(x) {
  output = read.csv(x, stringsAsFactors = FALSE)
  if (kit=="S0293689") {
    manifest$filename = gsub(".bam", "", manifest$filename)
    manifest$SampleID = str_extract(manifest$barcode, "TCGA.{16}")
    output$Sampleid = manifest[which(manifest$filename == gsub("^X", "", output$Sampleid)),]$SampleID 
    # output$submitter_id = str_extract(output$Sampleid, "TCGA.{11}")
    # output$Comment = gsub(".*;", "", output$Comment)
  }
  output = output[c("Sampleid", "Purity", "Ploidy", "Flagged", "Comment")]
  output$Sampleid = str_extract(output$Sampleid, "TCGA.{16}")
  output$Sampleid = gsub("\\.", "-", output$Sampleid)
  output$submitter_id = str_extract(output$Sampleid, "TCGA.{11}")
  output$Comment = gsub(".*;", "", output$Comment)
  return(output)
}