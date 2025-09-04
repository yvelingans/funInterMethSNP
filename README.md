# funInterMethSNP: Detecting SNP–Methylation Interaction Effects on Phenotypes

Understanding how **DNA methylation** and **genetic variants (SNPs)** jointly influence phenotypes is essential to uncover the molecular basis of complex diseases.  

**funInterMethSNP** provides a functional regression framework to **detect interaction effects** between SNP genotypes (coded as 0/1/2) and DNA methylation profiles across various phenotypes.  
A key step in these models is the reconstruction of DNA methylation data as functional variables.  

The package implements models for **continuous**, **binary**, and **count outcomes**, and includes tools to reconstruct methylation curves and perform global hypothesis testing.  

## Main features

- Models DNA methylation as a **functional predictor**.  
- Tests **global SNP–CpG interaction effects** rather than individual SNP–CpG pairs.  
- Accounts for **spatial correlation** among neighboring CpG sites.  
- Provides flexible kernel weighting with parameter `rho` to emphasize **localized effects**.  
- Supports outcomes from the **exponential family** (Gaussian, binomial, Poisson, etc.).  

## Installation

Install the development version from GitHub :

```r
# install.packages("devtools")
devtools::install_github("yvelingans/funInterMethSNP", build_vignettes = TRUE)
```


## Example usage

```r
library(funInterMethSNP)

# Load example data (derived from GSE73103)
data(data_methylation_raw)
data(mat_snp)
data(mat_covar)
data(resp_continu)
data(pos_snp)

# Reconstruct methylation curves
mat_methyl <- estim_methyl_curve(
  mat_methylation_raw = data_methylation_raw,
  region_range = c(min(data_methylation_raw$Position),
                   max(data_methylation_raw$Position))
)

# Fit the continuous outcome model
fit_cont <- fitMethSNPfunctInter_continu(
  Y = resp_continu,
  W = mat_covar,
  G = mat_snp,
  mat_methyl = mat_methyl,
  pos_snp = pos_snp,
  rho = 10,
  nb_fct = 10,
  kernel_name = "convex",
  region_range = c(min(data_methylation_raw$Position),
                   max(data_methylation_raw$Position)),
  RGN = seq(0, 1, length.out = ncol(mat_methyl))
)

# Inspect results
fit_cont$pvalue
fit_cont$teta 
summary(fit_cont$object_gam)
plot(fit_cont$object_gam)
```



## Documentation
For a complete tutorial, see the vignette :

```r
browseVignettes("funInterMethSNP")
```

## Citation
If you use **funInterMethSNP** in your research, please cite:  
Gansou, Y. et al. (2025). *A functional approach to testing the overall effect of interaction between DNA methylation and SNPs*. Submitted.



