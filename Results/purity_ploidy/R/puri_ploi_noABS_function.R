puri_ploi= function(x) {
  output = read.csv(x, stringsAsFactors = FALSE)
  if (kit=="S0293689") {
    manifest$filename = gsub(".bam", "", manifest$filename)
    manifest$SampleID = str_extract(manifest$barcode, "TCGA.\\d+.\\d+.\\d+")
    output$Sampleid = manifest[which(manifest$filename == gsub("^X", "", output$Sampleid)),]$barcode 
  }
  output = output[c("Sampleid", "Purity", "Ploidy")]
  output$Sampleid = str_extract(output$Sampleid, "TCGA.\\d+.\\d+.\\d+..\\d+.")
  output$Sampleid = gsub("\\.", "-", output$Sampleid)
  return(output)
}