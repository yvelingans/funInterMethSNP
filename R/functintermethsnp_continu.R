#' Functional Interaction Model between SNPs and Methylation for Continuous Outcomes
#'
#' @description
#' This function fits a functional regression model to assess the global interaction between SNP genotypes and DNA methylation profiles on continuous outcomes. A set of kernel-weighted interaction terms is constructed, and a functional smooth term is estimated using penalized splines. A Wald-type test is performed to evaluate the global significance of the interaction effects.
#' @import mgcv
#' @param Y Numeric vector: continuous response variable
#' @param W Matrix of covariates (individuals x covariates)
#' @param G Matrix of SNPs (individuals x SNPs, values in 0/1/2)
#' @param mat_methyl Matrix of methylation curves (individuals x grid points),
#'                   as obtained using the `estim_methyl_curve()` function.
#' @param pos_snp Real SNP positions
#' @param rho Numeric: kernel parameter controlling the weight decay
#' @param nb_fct Number of basis functions
#' @param kernel_name One of "convex", "gaussian", or "linear"
#' @param region_range Vector of two elements: start and end positions of CpGs in region of interest
#' @param RGN Numeric vector. Grid of values in \[0,1\] used for curve estimation
#'            (default = seq(0,1,length.out = 1000)).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{pvalue}: Global p-value of the interaction test,
#'   \item \code{teta}: Vector of estimated model coefficients,
#'   \item \code{sig2}: Estimated residual variance,
#'   \item \code{lambda}: Estimated smoothing parameter,
#'   \item \code{object_gam}: Object of class \code{gam} (from \pkg{mgcv}).
#'     Can be used with \code{summary()} to obtain individual tests,
#'     and \code{plot()} to visualize the estimated curve \eqn{\hat{\delta}(t)}.
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n <- 100  # nombre d'individus
#' p <- 3    # nombre de SNPs
#' m <- 500  # nombre de points CpG
#'
#' Y <- rnorm(n)
#' W <- matrix(rnorm(n * 2), ncol = 2)
#' colnames(W) <- c("age", "sex")
#'
#' G <- matrix(sample(0:2, n * p, replace = TRUE), ncol = p)
#' colnames(G) <- c("rs101", "rs202", "rs303")
#'
#' # Example of methylation curve estimation
#' mat_methyl <- matrix(runif(n * m), nrow = n)  # should be obtained beforehand.
#'
#' pos_snp_bp <- c(120000, 150000, 180000)
#' region_range <- c(100000, 200000)
#'
#' res <- fitMethSNPfunctInter_continu(
#'   Y = Y,
#'   W = W,
#'   G = G,
#'   mat_methyl = mat_methyl,
#'   pos_snp = pos_snp_bp,
#'   rho = 10,
#'   nb_fct = 10,
#'   kernel_name = "gaussian",
#'   region_range = region_range
#' )
#'
#' res$pvalue  # Pvalue for the globalinteraction test
#' }
#' @export
#' @importFrom mgcv gamm
#' @importFrom stats as.formula pf dnorm



fitMethSNPfunctInter_continu <- function(Y, W, G, mat_methyl, pos_snp,
                                      rho, nb_fct, kernel_name,region_range,
                                      RGN = seq(0, 1, length.out = 1000)){
  samp_size<-length(Y)

  # Check inputs
  if (!is.numeric(Y) || !is.vector(Y)) stop("Y must be a numeric vector.")

  if (!is.matrix(W)) stop("W must be a matrix.")

  if (!is.matrix(G)) stop("G must be a matrix.")

  if (!is.matrix(mat_methyl)) stop("mat_methyl must be a matrix return by fonction estim_methyl_curve.")

  if (!is.numeric(pos_snp)) stop("pos_snp must be a numeric vector.")

  if (!is.numeric(rho) || length(rho) != 1) stop("rho must be a single numeric value.")

  if (!is.numeric(nb_fct) || length(nb_fct) != 1) stop("nb_fct must be a single numeric value.")

  if (!kernel_name %in% c("convex", "gaussian", "linear")) stop("kernel_name must be one of: 'convex', 'gaussian', 'linear'.")

  # Check dimensions consistency
  if (length(Y) != nrow(W)) stop("Number of rows in W must match length of Y.")

  if (length(Y) != nrow(G)) stop("Number of rows in G must match length of Y.")

  if (length(Y) != nrow(mat_methyl)) stop("Number of rows in mat_methyl must match length of Y.")

  if (ncol(G) != length(pos_snp)) stop("Length of pos_snp must match number of columns in G.")

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

    # GAMM model  fitting using mgcv::gamm
    mod <- gamm(formule, data = data_gam, method = "REML")



    # Gam object extraction
    res_gam<-mod$gam
    object_gam<- res_gam

    # Variance component
    my_sigma2_int<- res_gam$sig2

    # Estimated smoothing parameter
    lambda_optim<-res_gam$sp

    # Vector of all coefficients
    teta_int<-as.vector(res_gam$coefficients)

    # Variance-covariance matrix of coefficients
    var_teta_int_vp<- res_gam$Vp

    # Extract estimated interaction coefficients
    eta_hat<-teta_int[start_eta:end_eta]

    #Variance-covariance matrix of interaction coefficients
    var_vp_etahat<- var_teta_int_vp[start_eta:end_eta,start_eta:end_eta]

    ## Compute Wald-type statistic
    Fisher_vp<- (t(eta_hat)%*%(solve(var_vp_etahat))%*%eta_hat)/nb_snp

    # Degree of freedom for denominator
    degden<-samp_size- my_nb_cof

    # Computation associated p-value for testing overall interaction
    my_pvalue_Fisher_vp<- pf(Fisher_vp,nb_snp,degden,lower.tail = FALSE)

    ## Output list
    resultat<-list(pvalue=my_pvalue_Fisher_vp,
                   teta=res_gam$coef,
                   lambda=lambda_optim,
                   sig2=my_sigma2_int,
                   object_gam=object_gam)
    return(resultat)
}

