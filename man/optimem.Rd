\name{optimem}
\alias{optimem}
\title{
Optimal Normalization Method for Sparse Compositional Data
}
\description{
Determine non-absolutely, differentially abundant (ADA) taxa that are used to remove the compositional effects using a sequential removal and random amalgamation procedure
}
\usage{
optimem(M, y, eta=0.1, alpha=0.05, n.k=NULL, n.b=200, n.r=300,
        min.nonADA=0.1, n.perm=20, k.sel.plot=TRUE)
}
\arguments{
  \item{M}{a matrix of compositional data}
  \item{y}{a vector of categorical outcomes}
  \item{eta}{a proportion of taxa removed at each removal step}
  \item{alpha}{a significance level used to determine the lower and upper bounds of the null distribution of MSS}
  \item{n.k}{a number of the removal step}
  \item{n.b}{a number of random selection of taxa at a removal step. It is inactive when eta*ncol(M) < 1}
  \item{n.r}{a number of random amalgamation}
  \item{min.nonADA}{a minimum number of non-ADA taxa. OPTIMEM stops when it is smaller than the number of the remaining taxa}
  \item{n.perm}{a number of permutations to generate null distributions of MSS}
  \item{k.sel.plot}{a logical value for generating a MSS vs k plot as an output}

  Note: OPTIMEM is computationally intensive. It may take several hours, depending on the number of taxa and the proportion of taxa removed at each removal step (eta). As a default, eta is set at 0.1. It has to be reduced if the computed MSS does not decrease as the sequential removal step increases.
}
\value{
  \item{nonADAtaxa_min}{a vector of nonADA taxa at the minimum of MSS}
  \item{min_MSS}{the minimum MSS}
  \item{mss_taxa}{MSS and corresponding nonADA taxa at each removal step}
  \item{mean_null}{a mean of the empirical distribution of MSS}
  \item{upper_limit_null}{an upper bound of the empirical distribution of MSS}
  \item{lower_limit_null}{a lower bound of the empirical distribution of MSS}
}
\references{
  Sohn, M.B., Monaco, C., and Gill, S.R. (2023). \emph{Optimal Normalization Method for Sparse Compositional Data} (Submitted)
}
\author{
  Michael B. Sohn

  Maintainer: Michael B. Sohn <michael_sohn@urmc.rochester.edu>
}
\examples{
\dontrun{
# Simulate a dataset
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

# Differential abundance analysis using the Wilcoxon signed-rank test
det.nonADA <- rslt$nonADAtaxa_min
if(rslt$min_MSS < rslt$lower_limit_null){
     MDA <- sweep(M[,!(colnames(M) \%in\% det.nonADA)], 1, rowSums(M[, det.nonADA]), "/")
} else{
     MDA <- sweep(M, 1, rowSums(M[, det.nonADA]), "/")
}
DA.test <- apply(MDA, 2, function(x) wilcox.test(x~y)$p.value)
DA.p <- p.adjust(DA.test, method="BH")
names(which(DA.p < 0.05))
}
}
