#!/bin/bash

# Rscript -e 'system.file("extdata", package = "PureCN")'
# [1] "/usr/local/lib/R/R-3.5-bioc-devel/PureCN/extdata"
export PURECN="/usr/local/lib/R/R-3.5-bioc-devel/PureCN/extdata"

MUTECT_OUT="/data/ovarian/CNVworkflow_LUAD/mutect_output"
PURECN_OUT="/data/ovarian/CNVworkflow_LUAD/purecn_output"
PREF="/home/sehyun/reference/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set"
REFERENCE="${PREF}.fna"
BEDFILES="/home/sehyun/reference/bedfiles"

CODE_DIR="/home/sehyun/Documents/github/PureCN_manuscript/Methods"
PURECN_CODE="${CODE_DIR}/runPureCN"
MUTECT_CODE="${CODE_DIR}/runMuTect"

# Prepare input files using MeTect
# --runtype argument takes two options: 'normal' or 'tumor'
Rscript ${MUTECT_CODE}/runMutect_stats.R --catalog ${catalog} --runtype ${runtype}
Rscript ${MUTECT_CODE}/runMutect_stats_matching_normal.R --catalog ${catalog} --runtype ${runtype}
Rscript ${MUTECT_CODE}/runMutect_normal_panel_vcf.R --catalog ${catalog} --runtype ${runtype}

# Prepare interval file
# check '--offtargetwidth' argument for NormalDB.R
Rscript ${PURECN_CODE}/lifeOver.R

Rscript ${PURECN}/IntervalFile.R \
--infile "${BEDFILES}/${BED}/${BED}_hg38.bed" \
--fasta ${REFERENCE} \
--outfile "${BEDFILES}/${BED}/${BED}_gcgene.txt" \
--offtarget --genome hg38 \
--export "${BEDFILES}/${BED}/${BED}_optimized.bed" \
--mappability "${PREF}_76.bw" \
--force

# Coverage files
# Calculate and GC-normalize coverage from a list of BAM files (both normal and tumor)
#/usr/bin/Rscript /home/sehyun/script/ovc_subtype_clonality/runCoverage_normal.R
#/usr/bin/Rscript /home/sehyun/script/ovc_subtype_clonality/runCoverage_tumor.R
Rscript runCoverage.R --catalog ${CATALOG} --kit ${KIT} --bed ${BED} --runtype "normal"
Rscript runCoverage.R --catalog ${CATALOG} --kit ${KIT} --bed ${BED} --runtype "tumor"

# Build normalDB
ls -a ${PURECN_OUT}/normal_cov/*_coverage_loess.txt | cat > ${PURECN_OUT}/normalDB/LUAD_normalDB.list

# Make an index file for normal.panel.vcf file
bgzip -c ${MUTECT_OUT}/${CATALOG}.normals.merged.min5_250.vcf > ${MUTECT_OUT}/${CATALOG}.normals.merged.min5_250.vcf.gz
tabix -p vcf ${MUTECT_OUT}/${CATALOG}.normals.merged.min5_250.vcf.gz

# From already GC-normalized files
Rscript $PURECN/NormalDB.R --outdir ${PURECN_OUT}/normalDB \
--coveragefiles ${PURECN_OUT}/normalDB/LUAD_normalDB.list \
--normal_panel ${MUTECT_OUT}/${CATALOG}.normals.merged.min5_250.vcf.gz \
--genome hg38 --force

# Run PureCN
Rscript ${PURECN_CODE}/runPureCN.R

# Analysis PureCN outputs
Rscript ${PURECN_CODE}/CallableLoci.R
Rscript ${PURECN_CODE}/runDx.R --catalog ${CATALOG} --runtype ${runtype}
