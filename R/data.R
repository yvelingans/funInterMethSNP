#' Example DNA methylation data
#'
#' A  dataset containing raw beta methylation values for CpG sites.
#' Values range between 0 and 1, extracted from the public dataset GSE73103.
#'
#' @format A data frame with CpG sites in rows and samples in columns.
"data_methylation_raw"

#' Example SNP genotype matrix
#'
#' SNP genotypes coded as 0, 1, 2 for individuals.
#'
#' @format A matrix with individuals in rows and SNPs in columns.
"mat_snp"

#' Example covariate matrix
#'
#' Contains potential confounders such as age and sex.
#'
#' @format A matrix with individuals in rows and covariates in columns.
"mat_covar"

#' Example SNP base pair positions
#'
#' Physical base-pair positions of the SNPs.
#'
#' @format A numeric vector with length equal to the number of SNPs.
"pos_snp"

#' Example methylation curves
#'
#' Estimated methylation curves for individuals over a genomic region.
#'
#' @format A matrix with individuals in rows and grid points in columns.
"mat_methyl"

#' Simulated continuous response
#'
#' A  continuous outcome variable.
#'
#' @format A numeric vector of length equal to the number of individuals.
"resp_continu"

#' Binary response (0/1)
#'
#' A  binary outcome variable (e.g., obesity status).
#'
#' @format A numeric vector of length equal to the number of individuals, with values 0 or 1.
"resp_binaire"

#' Simulated Poisson response (counts)
#'
#' A  count outcome variable generated from a Poisson model.
#'
#' @format A numeric vector of length equal to the number of individuals.
"resp_pois"
