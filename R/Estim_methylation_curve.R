#' Estimate methylation curves per individual using adaptive Gaussian kernel smoothing
#'
#' This function estimates smoothed methylation trajectories for each individual
#' from raw methylation values in \[0,1\], using adaptive kernel smoothing based on CpG positions.
#'
#' @param mat_methylation_raw Data frame of raw methylation data. CpGs in rows, individuals in columns,
#'                            with a column named 'Position' indicating genomic positions (in base pairs).
#' @param H Numeric. Minimum bandwidth for adaptive kernel smoothing (default = 1000 bp).
#' @param region_range Numeric vector of length 2. Genomic region of interest: c(start_position, end_position).
#' @param RGN Numeric vector. Grid of values in \[0,1\] used for curve estimation
#'            (default = seq(0,1,length.out = 1000)).
#'
#' @return A matrix of smoothed methylation curves with individuals in rows and grid points in columns.
#' @export
#' @examples
#' \dontrun{
#' # Simulate methylation data for 3 individuals and 500 CpGs
#' set.seed(123)
#' n_cpg <- 500
#' n_ind <- 3
#'
#' # Generate increasing positions within a fictional genomic region
#' positions <- sort(sample(11190000:11460000, n_cpg))
#'
#' # Generate beta values between 0 and 1 for each individual
#' mat_beta <- matrix(runif(n_cpg * n_ind, min = 0, max = 1), ncol = n_ind)
#'
#' # Combine beta values with Position column
#' mat_df <- as.data.frame(cbind(mat_beta, Position = positions))
#' colnames(mat_df) <- c(paste0("Ind", 1:n_ind), "Position")
#'
#' # Apply the smoothing function
#' smoothed <- estim_methyl_curve(
#'   mat_methylation_raw = mat_df,
#'   H = 1000,
#'   region_range = c(11190000, 11460000),
#'   RGN = seq(0, 1, length.out = 100)
#' )
#'}
#' @export
#' @importFrom stats dnorm
#'
estim_methyl_curve<-function(mat_methylation_raw,H=1000,
                             region_range=c(11190000,11460000),
                             RGN=seq(0,1,length.out=1000)){

  if (!is.data.frame(mat_methylation_raw)) {
    stop("'mat_methylation_raw' must be a data.frame with columns for individuals and a 'Position' column.")
  }


  if (!("Position" %in% colnames(mat_methylation_raw))) {
    stop("The matrix 'mat_methylation_raw' must include a column named 'Position'.")
  }

  if (!is.numeric(H) || H <= 0) stop("H must be a positive numeric value.")

  if (!is.numeric(region_range) || length(region_range) != 2 || region_range[1] >= region_range[2]) {
    stop("'region_range' must be a numeric vector of length 2 with region_range[1] < region_range[2].")
  }

  # exclude Position column
  samp_size<-ncol(mat_methylation_raw)- 1

  pos<-mat_methylation_raw$Position

  ## Estimate methylation at one position
  Estim.Methyl.elem <- function(t_val,theta_met,pos,H,region_range)
  {
    genomic_pos <- region_range[1] + t_val * (region_range[2] - region_range[1])
    dists <- abs(genomic_pos - pos)
    h <- max(sort(dists)[70], H)
    K <- dnorm((genomic_pos - pos)/h)
    Estim <- sum((theta_met)*K)/sum(K)
    return(Estim)
  }

  # For a vector of position
  Estim.Methyl <- function(theta_met,pos,H,region_range, RGN)
  {
    sapply(RGN,Estim.Methyl.elem,theta_met=theta_met,pos=pos,H=H,region_range = region_range)
  }

  # Matrix of smoothed methylation
  mat_methyl<-matrix(NA,nrow = samp_size, ncol = length(RGN))

  # Estimate methylation curves for all individuals
  for (i in 1:samp_size) {
    mat_methyl[i,]<-Estim.Methyl(theta_met = mat_methylation_raw[, i],
                                 pos = pos,
                                 H=H,
                                 region_range=region_range,
                                 RGN=RGN)
  }
  rownames(mat_methyl)<-rownames(mat_methylation_raw)
  return(mat_methyl)
}
