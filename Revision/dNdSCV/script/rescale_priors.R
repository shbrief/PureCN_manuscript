# I didn’t provide cosmic.vcf.file for runAbsoluteCN. --> dbSNP contains known cancer
# mutations. In tumor-only, when you provide a COSMIC VCF or a Cosmic.CNT info field 
# in the VCF, it checks if a variant is in both COSMIC and dbSNP and handles it appropriately.
# If only in dbSNP, it’s marked as likely germline and is ignored.


# take the CNT filtered COSMIC VCF and re-sets the prior for all overlapping variants
rescale_priors <- function(x, id, prior.somatic) {
    likelihoods.old <- x$results[[id]]$SNV.posterior$likelihood
    posteriors.old <- x$results[[id]]$SNV.posterior$posteriors
    prior.somatic.old <- posteriors.old$prior.somatic
    idx_somatic <- grep("SOMATIC.M", colnames(likelihoods.old))
    idx_germline <- grep("GERMLINE.M", colnames(likelihoods.old))
    likelihoods.new <- likelihoods.old
    likelihoods.new[, idx_somatic] <- likelihoods.new[,idx_somatic] / prior.somatic.old * prior.somatic
    likelihoods.new[, idx_germline] <- likelihoods.new[,idx_germline] / (1-prior.somatic.old) * (1-prior.somatic)
    posteriors.new <- likelihoods.new/pmax(rowSums(likelihoods.new), .Machine$double.xmin)
    xx <- PureCN:::.extractMLSNVState(posteriors.new)
    posteriors <- posteriors.old
    posteriors[, colnames(posteriors.new)] <- posteriors.new
    posteriors[, colnames(xx)] <- xx
    p <- x$results[[id]]$purity
    posteriors$ML.AR <- (p * posteriors$ML.M + ifelse(posteriors$ML.SOMATIC,
                                                      0, 1) * (1 - p))/(p * posteriors$ML.C + 2 * (1 - p))
    posteriors$ML.AR[posteriors$ML.AR > 1] <- 1
    depth <- posteriors$depth
    ar <- posteriors$AR.ADJUSTED
    ar[!posteriors$ML.SOMATIC] <- NA
    m <- t(apply(cbind(ar, depth, posteriors$ML.C), 1, function(x) PureCN:::.calculate_ccf(vaf = x[1],
                                                                                           depth = x[2], purity = p, C = x[3])))
    posteriors$CELLFRACTION <- as.numeric(m[, 1])
    posteriors$CELLFRACTION.95.LOWER <- as.numeric(m[, 2])
    posteriors$CELLFRACTION.95.UPPER <- as.numeric(m[, 3])
    posteriors$prior.somatic <- prior.somatic
    posteriors
}   

rescale_overlaps <- function(x, vcf, new.prior = 0.5) {
    p <- GRanges(predictSomatic(x))
    idx <- overlapsAny(p, vcf)
    prior.somatic <- p$prior.somatic
    prior.somatic[idx] <- new.prior
    rescale_priors(x, 1, prior.somatic)
}