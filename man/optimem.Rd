\name{optimem}
\alias{optimem}
\title{
Optimal Normalization Method for Sparse Compositional Data
}
\description{
Determine non-differentially abundant (non-DA) features (e.g., taxa) that are used to remove the compositional effects using a sequential removal and random amalgamation procedure
}
\usage{
optimem(M, y, eta=0.1, alpha=0.05, n.k=NULL, n.b=200, n.r=300,
        min.nonDA=0.1, n.perm=20, k.sel.plot=TRUE, early.exit=TRUE, k.min.prop=0.1)
}
\arguments{
  \item{M}{a matrix of compositional data}
  \item{y}{a vector of categorical outcomes}
  \item{eta}{a proportion of taxa removed at each removal step}
  \item{alpha}{a significance level used to determine the lower and upper bounds of the null distribution of MSS}
  \item{n.k}{a number of the removal step}
  \item{n.b}{a number of random selection of feature at a removal step. It is inactive when eta*ncol(M) < 1}
  \item{n.r}{a number of random amalgamation}
  \item{min.nonDA}{a minimum number of non-DA feature. OPTIMEM stops when it is smaller than the number of the remaining taxa}
  \item{n.perm}{a number of permutations to generate null distributions of MSS}
  \item{k.sel.plot}{a logical value for generating a MSS vs k plot as an output}
  \item{early.exit}{a logical value for early exit when the first 20\% of features are within a prespecified \% CI of the null distribution or when MSS is smaller than the lower limit of the null distribution and then larger than the upper limit of the null}
  \item{k.min.prop}{a minimum proportion of features required to continue the removal step}

  Note: OPTIMEM is computationally intensive. It may take a few hours, depending on the number of features (e.g., taxa) and the proportion of removed features at each removal step (eta). As a default, eta is set at a very small value (< 1/(number of features)), so one feature is removed at a time. It has to be increased if the number of features is greater than several hundreds. In this case, we recommend using large numbers of the random selection and amalgamation steps (e.g., n.b=500, n.r=1000). Note n.b will be ignored when eta is set at a very small value.
}
\value{
  \item{nonDAtaxa_min}{a vector of nonDA taxa at the minimum of MSS}
  \item{min_MSS}{the minimum MSS}
  \item{mss_taxa}{MSS and corresponding nonDA taxa at each removal step}
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

# Differential abundance analysis using the Wilcoxon signed-rank test
det.nonDA <- rslt$nonDAtaxa_min
if(rslt$min_MSS < rslt$lower_limit_null){
     MDA <- sweep(M[,!(colnames(M) %in% det.nonDA)], 1, rowSums(M[, det.nonDA]), "/")
} else{
     MDA <- sweep(M, 1, rowSums(M[, det.nonDA]), "/")
}
DA.test <- apply(MDA, 2, function(x) wilcox.test(x~y)$p.value)
DA.p <- p.adjust(DA.test, method="BH")
names(which(DA.p < 0.05))
}
}
