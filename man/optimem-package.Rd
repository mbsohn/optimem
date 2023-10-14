\name{optimem-package}
\alias{optimem-package}
\docType{package}
\title{
Optimal Normalization Method for Sparse Compositional Data
}
\description{
Determine non-differentially abundant (non-DA) features (e.g., taxa) that are used to remove the compositional effects using a sequential removal and random amalgamation procedure
}
\author{
Michael B. Sohn

Maintainer: Michael B. Sohn <michael_sohn@urmc.rochester.edu>
}
\references{
Sohn, M.B., Monaco, C., and Gill, S.R. \emph{Optimal Normalization Method for Sparse Compositional Data} (Submitted)
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

# Plot MSS vs k
plot(rslt)

# Print OPTIMEM outputs
print(rslt)

# Differential abundance analysis using the Wilcoxon signed-rank test
det.nonDA <- rslt$nonDAtaxa_min
if(rslt$min_MSS < rslt$lower_limit_null){
     MDA <- sweep(M[,!(colnames(M) \%in\% det.nonDA)], 1, rowSums(M[, det.nonDA]), "/")
} else{
     MDA <- sweep(M, 1, rowSums(M[, det.nonDA]), "/")
}
DA.test <- apply(MDA, 2, function(x) wilcox.test(x~y)$p.value)
DA.p <- p.adjust(DA.test, method="BH")
names(which(DA.p < 0.05))
}
}
