README
================
Michael B. Sohn
7/12/2023

## OPTIMEM: Optimal Normalization Method for Sparse Compositional Data

The optimem function determines non-absolutely, differentially abundant
(ADA) taxa that are used to remove the compositional effects using a
sequential removal and random amalgamation procedure. For detailed
information about the arguments, please see the documentation for
*optimem()*.

Note: OPTIMEM is computationally intensive. It may take several hours,
depending on the number of taxa and the proportion of removed taxa at
each removal step (eta). As a default, eta is set at 0.1. It has to be
reduced if the computed MSS does not decrease as the sequential removal
step increases. We also recommend using large numbers of the random
selection and amalgamation steps (e.g., n.b=500, n.r=1000). When eta is
set at a very small value (\< 1/(number of taxa)), one taxon will be
removed at each removal step and n.b will be ignored.

### Installation

install.packages(“devtools”)

devtools::install_github(“mbsohn/optimem”)

### Example: Determine non-ADA taxa

``` r
library(optimem)
# Simulate a dataset using a negative binomial model
set.seed(2023)
p <- 100; n.non.da <- sample(40:90, 1)
non.da.mu <- sample(1:100, n.non.da, rep=TRUE)
da.mu1 <- sample(1:100, p-n.non.da, rep=TRUE)
da.mu2 <- sample(1:100, p-n.non.da, rep=TRUE)
mu1 <- c(da.mu1, non.da.mu); mu2 <- c(da.mu2, non.da.mu)
sz <- 1; n1 <- 50; n2 <- 50
dat1 <- t(replicate(n1, MASS::rnegbin(length(mu1), mu1*sample(1:10, 1), sz)))
dat2 <- t(replicate(n2, MASS::rnegbin(length(mu1), mu2*sample(1:10, 1), sz)))
M <- proportions(rbind(dat1, dat2), margin=1)
n.sample <- nrow(M); n.taxa <- ncol(M)
y <- c(rep(1, n1), rep(2, n2))
colnames(M) <- paste0("T", 1:n.taxa)
true.da.taxa <- colnames(M)[which(mu1 != mu2)]
# Run OPTIMEM
rslt <- optimem(M, y)
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

### Differential abundance analysis using non-ADA taxa as a reference

``` r
det.nonADA <- rslt$nonADAtaxa_min
if(rslt$min_MSS < rslt$lower_limit_null){
     MDA <- sweep(M[,!(colnames(M) %in% det.nonADA)], 1, rowSums(M[, det.nonADA]), "/")
} else{
     MDA <- sweep(M, 1, rowSums(M[, det.nonADA]), "/")
}
ADA.test <- apply(MDA, 2, function(x) wilcox.test(x~y)$p.value)
ADA.p <- p.adjust(ADA.test, method="BH")
names(which(ADA.p < 0.05))
```

    ##  [1] "T1"  "T2"  "T5"  "T9"  "T11" "T13" "T16" "T18" "T20" "T25" "T26" "T29"
    ## [13] "T34" "T35" "T36" "T37" "T39" "T40" "T41" "T42" "T44" "T45"

To account for covariates, a probabilistic index model, which can be
seen as the rank-equivalent of the general linear model, can be used.
