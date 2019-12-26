# Build reference for dNdSCV
wd = "~/data2/PureCN_manuscript/Revision/dNdSCV"
path_cds_table = file.path(wd, "reference/mart_export.txt")
path_genome_fasta = "~/reference/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Build GRCh38 reference
buildref(cdsfile = path_cds_table, 
         genomefile = path_genome_fasta, 
         outfile = file.path(wd, "reference/GRCh38_refcds.rda"), 
         excludechrs = "MT")