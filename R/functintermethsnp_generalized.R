#' Global Test for SNP–Methylation Interactions on Generalized Outcomes via Functional Regression
#'
#' This function fits a functional regression to assess the global interaction between SNP genotypes and DNA methylation profiles on non-Gaussian outcomes, such as binary or count responses.
#' A likelihood ratio test is performed to evaluate the significance of the interaction terms.
#' @param Y Response variable (vector), assumed from exponential family (e.g., binary, Poisson)
#' @param W Matrix of covariates (individuals x covariates)
#' @param G Matrix of SNPs (individuals x SNPs, values in 0/1/2)
#' @param mat_methyl Matrix of methylation curves (individuals x grid points)
#' @param pos_snp SNP base pair positions (numeric vector)
#' @param rho Numeric: kernel parameter controlling the weight decay
#' @param nb_fct Number of basis functions in the GAM
#' @param kernel_name One of "convex", "gaussian", or "linear"
#' @param region_range Vector of two elements: start and end positions in bp of CpGs in region of interest
#' @param family A family object describing the error distribution (e.g., binomial(), poisson())
#' @param RGN Numeric vector. Grid of values in \[0,1\] used for curve estimation
#'            (default = seq(0,1,length.out = 1000)).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{pvalue}: p-value from the likelihood ratio test comparing models with and without interaction terms,
#'   \item \code{teta}: vector of estimated coefficients from the full model,
#'   \item \code{lambda}: Estimated smoothing parameter.,
#'   \item \code{object_gam}: the fitted \code{mgcv::gam} object (full model).
#' }
#' @examples
#' \dontrun{
#' # Simuler des données
#' set.seed(123)
#' n <- 100  # nombre d'individus
#' p <- 3    # nombre de SNPs
#' m <- 500  # nombre de points CpG
#'
#' Y <- rbinom(n, size = 1, prob = 0.4)
#' W <- matrix(rnorm(n * 2), ncol = 2)
#' colnames(W) <- c("age", "sex")
#'
#' G <- matrix(sample(0:2, n * p, replace = TRUE), ncol = p)
#' colnames(G) <- c("rs101", "rs202", "rs303")
#'
#' mat_methyl <- matrix(runif(n * m), nrow = n)
#'
#' pos_snp_bp <- c(120000, 150000, 180000)
#' region_range <- c(100000, 200000)
#'
#' res <- fitMethSNPfunctInter_generalized(
#'   Y = Y,
#'   W = W,
#'   G = G,
#'   mat_methyl = mat_methyl,
#'   pos_snp = pos_snp_bp,
#'   rho = 10,
#'   nb_fct = 10,
#'   kernel_name = "gaussian",
#'   region_range = region_range,
#'   family=binomial()
#' )
#'
#' res$pvalue  #pvalue for the global interaction testing
#' }
#' @export
#' @importFrom mgcv gam  gamm
#' @importFrom stats as.formula binomial pchisq logLik coef



fitMethSNPfunctInter_generalized <- function(Y, W, G, mat_methyl, pos_snp,
                                         rho, nb_fct, kernel_name,region_range,
                                         RGN = seq(0, 1, length.out = 1000),
                                         family=binomial()){
  samp_size<-length(Y)

  # Check inputs
  if (!is.numeric(Y) || !is.vector(Y)) stop("Y must be a numeric vector.")

  if (!is.matrix(W)) stop("W must be a matrix.")

  if (!is.matrix(G)) stop("G must be a matrix.")

  if (!is.matrix(mat_methyl)) stop("mat_methyl must be a matrix return by the function estim_methyl_curve.")

  if (!is.numeric(pos_snp)) stop("pos_snp must be a numeric vector.")

  if (!is.numeric(rho) || length(rho) != 1) stop("rho must be a single numeric value.")

  if (!is.numeric(nb_fct) || length(nb_fct) != 1) stop("nb_fct must be a single numeric value.")

  if (!kernel_name %in% c("convex", "gaussian", "linear")) stop("kernel_name must be one of: 'convex', 'gaussian', 'linear'.")

  if (is.null(colnames(W))) colnames(W) <- paste0("W", 1:ncol(W))

  if (is.null(colnames(G))) colnames(G) <- paste0("SNP", 1:ncol(G))

  # Check dimensions consistency
  if (length(Y) != nrow(W)) stop("Number of rows in W must match length of Y.")

  if (length(Y) != nrow(G)) stop("Number of rows in G must match length of Y.")

  if (length(Y) != nrow(mat_methyl)) stop("Number of rows in mat_methyl must match length of Y.")

  if (ncol(G) != length(pos_snp)) stop("Length of pos_snp must match number of columns in G.")

  # Check that 'family' is a valid family object (e.g., binomial(), poisson(), Gamma(), etc.)
  if (inherits(family, "family")) {
    fam_name <- family$family
  } else {
    stop("The 'family' argument must be a valid family object, e.g., binomial(), poisson(), Gamma().")
  }

  # Reject gaussian family — should be used with fitMethSNPfunctInter_continu()
  if (fam_name == "gaussian") {
    stop("Error: 'gaussian' family is not supported in this function. Please use fitMethSNPfunctInter_continu() for continuous outcomes.")
  }
  ##-------------------------------
  ##--------- Kernel function -----

  kernel_convex<-function(t,pos_snp,rho){
    exp((-rho*abs(t-pos_snp)))
  }

  kernel_gaussian<-function(t,pos_snp,rho){
    exp((-(rho^2)*(t-pos_snp)^2))
  }

  Kernel_linear <- function(t, pos_snp, rho) {
    pmax(1 - rho * abs(t - pos_snp), 0)
  }


  # compute the integral of the product of  kernel_convex and methylation curve
  intg_kernelvex_methyl<-function(id,pos_snp,rho,mat_methyl,RGN=seq(0,1,length.out=1000)){
    intg_kernelvex_methyl_elm<-function(pos_snp,rho,id,mat_methyl,RGN){
      (sum(kernel_convex(t=RGN,pos_snp = pos_snp,rho = rho)*mat_methyl[id,]))/length(RGN)
    }
    sapply(id,intg_kernelvex_methyl_elm,pos_snp=pos_snp,rho=rho,mat_methyl=mat_methyl,RGN=RGN)
  }

  # compute the integral of the product of  kernel_gaussian and methylation curve

  intg_kernelgaus_methyl<-function(id,pos_snp,rho,mat_methyl,RGN=seq(0,1,length.out=1000)){
    intg_kernelgaus_methyl_elm<-function(pos_snp,rho,id,mat_methyl,RGN){
      (sum(kernel_gaussian(t=RGN,pos_snp = pos_snp,rho = rho)*mat_methyl[id,]))/length(RGN)
    }
    sapply(id,intg_kernelgaus_methyl_elm ,pos_snp=pos_snp,rho=rho,mat_methyl=mat_methyl,RGN=RGN)
  }

  # compute the integral of the product of  kernel_linear and methylation

  intg_kernellin_methyl<-function(id,pos_snp,rho,mat_methyl,RGN=seq(0,1,length.out=1000)){
    intg_kernellin_methyl_elm<-function(pos_snp,rho,id,mat_methyl,RGN){
      (sum(Kernel_linear(t=RGN,pos_snp = pos_snp,rho = rho)*mat_methyl[id,]))/length(RGN)
    }
    sapply(id,intg_kernellin_methyl_elm,pos_snp=pos_snp,rho=rho,mat_methyl=mat_methyl,RGN=RGN)
  }

  # compute and return matrix of integrals of the  product of kernel_function and methylation,
  # ind in rows, SNPs positions in columns

  mat_intg_kernel_methyl<-function(pos_snp,id,rho,mat_methyl,kernel_name,RGN=seq(0,1,length.out=1000)){
    if (kernel_name=="convex"){
      all_intg<- as.vector(sapply(pos_snp,intg_kernelvex_methyl,id=id,rho=rho,mat_methyl=mat_methyl,RGN=RGN))
    }else if (kernel_name=="gaussian"){
      all_intg<- as.vector(sapply(pos_snp,intg_kernelgaus_methyl,id=id,rho=rho,mat_methyl=mat_methyl,RGN=RGN))
    }else if (kernel_name=="linear"){
      all_intg<- as.vector(sapply(pos_snp,intg_kernellin_methyl,id=id,rho=rho,mat_methyl=mat_methyl,RGN=RGN))
    }else{
      stop("Unsupported kernel")
    }
    nbr_snp<-length(pos_snp)
    mat_intg<-matrix(all_intg,ncol = nbr_snp)
    return(mat_intg)
  }
  #--------------------------------------------------
  #----------------Model fitting---------------------

  ## number of total coefficients in the model
  my_nb_cof<-1+ncol(W)+2*ncol(G)+nb_fct

  # Index of the first coefficient associated with interaction terms
  start_eta<-ncol(W)+ncol(G)+2

  # Index of the last coefficient associated with interaction terms
  end_eta<-my_nb_cof-nb_fct

  # Number of SNPs
  nb_snp <-ncol(G)


  if (!is.numeric(region_range) || length(region_range) != 2) {
    stop("region_range must be a numeric vector of length 2.")
  }

  if (region_range[1] >= region_range[2]) {
    stop("The first element of region_range must be less than the second.")
  }

  # Normalized SNPs positions
  pos_snp_normalize<-(pos_snp- region_range[1])/(region_range[2]-region_range[1])

  # Compute the matrix of integrals between product of kernel function and methylation curve
  mat_Omegarho_model<-mat_intg_kernel_methyl(pos_snp = pos_snp_normalize,id=1:samp_size,rho = rho,
                                             mat_methyl = mat_methyl,kernel_name = kernel_name,RGN=RGN)
  # Matrix of interaction terms (element-wise product)
  mat_Krho_model<-mat_Omegarho_model*G

  # Matrix of methylation curve scaled from GAMM
  mat_meth_gam<-mat_methyl/(length(RGN))

  # Matrix of grid points (replicated for each individual)
  mat_RGN<-matrix(RGN,samp_size,length(RGN),byrow = TRUE)

  names_W <- colnames(W)
  names_G <- colnames(G)
  names_Krho <- paste0("inter_meth_", colnames(G))
  colnames(mat_Krho_model) <- names_Krho



  # Design matrix for gamm
  data_gam <- data.frame(
    Y = Y,
    W,
    G,
    mat_Krho_model,
    mat_meth_gam = I(mat_meth_gam),
    mat_RGN = I(mat_RGN)
  )

  formule_text <- paste(
    "Y ~",
    paste(c(names_W, names_G, names_Krho), collapse = " + "),
    "+ s(mat_RGN, by = mat_meth_gam, k = ", nb_fct, ", m = 2, bs = 'ps')"
  )
  formule <- as.formula(formule_text)

  # Model without interaction (null model)
  formule_H0 <- as.formula(paste(
    "Y ~",
    paste(c(names_W, names_G), collapse = " + "),
    "+ s(mat_RGN, by = mat_meth_gam, k =", nb_fct, ", m = 2, bs = 'ps')"
  ))
  # Fit null model (no SNP-methylation interaction)
  mod_H0 <- gam(formule_H0, data = data_gam, method = "REML", family = family)

  # Fit full model (with interaction terms
  mod_H1 <- gam(formule, data = data_gam, method = "REML", family = family)

  # Perform likelihood ratio test (LRT)
  ll_H0 <- logLik(mod_H0)
  ll_H1 <- logLik(mod_H1)
  lr_stat <- -2 * (ll_H0 - ll_H1)
  df <- length(coef(mod_H1)) - length(coef(mod_H0))
  p_value_lrt <- as.numeric(pchisq(lr_stat, df, lower.tail = FALSE))

  resultat <- list(
    pvalue = p_value_lrt,
    teta = mod_H1$coef,
    lambda=mod_H1$sp,
    object_gam = mod_H1
  )
  return(resultat)
}


