## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(seqgendiff)
library(limma)
set.seed(31) ## for reproducibility

## ------------------------------------------------------------------------
nsamp <- 100
ngene <- 10000
mat <- matrix(stats::rpois(n = nsamp * ngene, lambda = 100), 
              nrow = ngene,
              ncol = nsamp)

## ------------------------------------------------------------------------
submat <- select_counts(mat = mat, nsamp = 6, ngene = 1000)

## ------------------------------------------------------------------------
head(rownames(submat))
head(colnames(submat))

## ------------------------------------------------------------------------
scaling_factor <- seq_len(ncol(submat))
scaling_factor
thout <- thin_lib(mat = submat, thinlog2 = scaling_factor)

## ------------------------------------------------------------------------
## Empirical thinning
colSums(thout$mat) / colSums(submat)

## Specified thinning
2 ^ -scaling_factor

## ------------------------------------------------------------------------
thout <- thin_all(mat = submat, thinlog2 = 1)

## ------------------------------------------------------------------------
sum(thout$mat) / sum(submat)

## ------------------------------------------------------------------------
designmat <- cbind(rep(c(0, 1), each = ncol(submat) / 2), 
                   rep(c(0, 1), length.out = ncol(submat)))
designmat

coefmat <- matrix(stats::rnorm(ncol(designmat) * nrow(submat)),
                  ncol = ncol(designmat),
                  nrow = nrow(submat))
head(coefmat)

## ------------------------------------------------------------------------
thout <- thin_diff(mat          = submat, 
                   design_fixed = designmat, 
                   coef_fixed   = coefmat)

## ------------------------------------------------------------------------
new_design <- cbind(thout$design_obs, thout$designmat)
vout <- limma::voom(counts = thout$mat, design = new_design)
lout <- limma::lmFit(vout)
coefhat <- coef(lout)[, -1, drop = FALSE]

## ------------------------------------------------------------------------
oldpar <- par(mar = c(2.5, 2.5, 1, 0) + 0.1,
              mgp = c(1.5, 0.5, 0))
plot(x    = coefmat[, 1], 
     y    = coefhat[, 1], 
     xlab = "True Coefficient", 
     ylab = "Estimated Coefficient",
     main = "First Variable",
     pch  = 16)
abline(a   = 0,
       b   = 1,
       lty = 2, 
       col = 2,
       lwd = 2)

plot(x    = coefmat[, 2], 
     y    = coefhat[, 2], 
     xlab = "True Coefficient", 
     ylab = "Estimated Coefficient",
     main = "Second Variable",
     pch  = 16)
abline(a   = 0,
       b   = 1,
       lty = 2, 
       col = 2,
       lwd = 2)
par(oldpar)

## ------------------------------------------------------------------------
target_cor <- matrix(c(0.9, 0), nrow = 2)
target_cor
thout_cor <- thin_diff(mat         = submat,
                       design_perm = designmat, 
                       coef_perm   = coefmat, 
                       target_cor  = target_cor)

## ------------------------------------------------------------------------
cor(thout_cor$designmat, thout_cor$sv)

## ------------------------------------------------------------------------
eout <- effective_cor(design_perm = designmat, 
                      sv          = thout_cor$sv, 
                      target_cor  = target_cor, 
                      iternum     = 50)
eout

## ------------------------------------------------------------------------
thout <- thin_2group(mat           = submat, 
                     prop_null     = 0.9, 
                     signal_fun    = stats::rgamma,
                     signal_params = list(shape = 1, rate = 1))

## ------------------------------------------------------------------------
new_design <- cbind(thout$design_obs, thout$designmat)
new_design
vout <- limma::voom(counts = thout$mat, design = new_design)
lout <- limma::lmFit(vout)
coefhat <- coef(lout)[, 2, drop = FALSE]

## ------------------------------------------------------------------------
oldpar <- par(mar = c(2.5, 2.5, 1, 0) + 0.1,
              mgp = c(1.5, 0.5, 0))
plot(x    = thout$coefmat, 
     y    = coefhat, 
     xlab = "True Coefficient", 
     ylab = "Estimated Coefficient",
     main = "First Variable",
     pch  = 16)
abline(a   = 0,
       b   = 1,
       lty = 2, 
       col = 2,
       lwd = 2)

