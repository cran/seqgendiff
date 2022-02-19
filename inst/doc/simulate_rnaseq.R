## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  echo = requireNamespace("airway", quietly = TRUE),
  eval = requireNamespace("airway", quietly = TRUE),
  comment = "#>",
  fig.align = "center",
  fig.height = 7,
  fig.width = 7
)

## ---- eval = FALSE, echo = !requireNamespace("airway", quietly = TRUE)--------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("airway")

## -----------------------------------------------------------------------------
set.seed(1)

## ---- eval = FALSE------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("airway")

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(airway)
data("airway")

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(seqgendiff)
library(SummarizedExperiment)
library(DESeq2)
library(limma)
library(sva)
library(qvalue)

## -----------------------------------------------------------------------------
coldat <- colData(airway)[, c("cell", "dex")]
true_sv <- model.matrix(~cell + dex, data = coldat)[, -1]
true_sv

## -----------------------------------------------------------------------------
allzero <- rowSums(assay(airway)) < 10^-6
airway <- airway[!allzero, ]

## -----------------------------------------------------------------------------
thout <- thin_2group(mat = assay(airway), 
                     prop_null = 0.9, 
                     signal_fun = stats::rnorm,
                     signal_params = list(mean = 0, sd = 0.8))

## -----------------------------------------------------------------------------
X <- cbind(thout$design_obs, thout$designmat)
Y <- log2(thout$mat + 0.5)
n_sv <- num.sv(dat = Y, mod = X)
svout <- sva(dat = Y, mod = X, n.sv = n_sv)
vout <- voom(counts = thout$mat, design = cbind(X, svout$sv))
lout <- lmFit(vout)
eout <- eBayes(lout)
qout <- qvalue(p = eout$p.value[, 2])
bhat <- eout$coefficients[, 2]

## -----------------------------------------------------------------------------
plot(thout$coefmat, 
     bhat, 
     xlab = "True Coefficients", 
     ylab = "Estimated Coefficients")
abline(0, 1, col = 2, lwd = 2)

## -----------------------------------------------------------------------------
is_null_gene <- abs(thout$coefmat) < 10^-6
boxplot(qout$qvalues ~ is_null_gene,
        xlab = "Null Gene",
        ylab = "q-value")

## ---- message=FALSE-----------------------------------------------------------
thout$design_obs <- cbind(thout$design_obs, svout$sv)
dds <- ThinDataToDESeqDataSet(obj = thout)
colData(dds)
design(dds)

## -----------------------------------------------------------------------------
dds <- DESeq(object = dds)
pval_dds <- rowData(dds)[, "WaldPvalue_P1"]
qval_dds <- qvalue(p = pval_dds)

## -----------------------------------------------------------------------------
coefmat <- rowData(dds)[, c("true_P1", "P1")]
plot(coefmat$true_P1, 
     coefmat$P1,
     xlab = "True Coefficients", 
     ylab = "Estimated Coefficients")
abline(0, 1, col = 2, lwd = 2)

## -----------------------------------------------------------------------------
boxplot(qval_dds$qvalues ~ is_null_gene,
        xlab = "Null Gene",
        ylab = "q-value")

## -----------------------------------------------------------------------------
mse_limma <- mean((bhat - thout$coefmat) ^ 2)
mse_deseq2 <- mean((coefmat$true_P1 - coefmat$P1) ^ 2, na.rm = TRUE)
mse_deseq2
mse_limma

## -----------------------------------------------------------------------------
mean(is_null_gene[qval_dds$qvalues < 0.1], na.rm = TRUE)
mean(is_null_gene[qout$qvalues < 0.1], na.rm = TRUE)

## -----------------------------------------------------------------------------
sum(qval_dds$qvalues < 0.1, na.rm = TRUE)
sum(qout$qvalues < 0.1, na.rm = TRUE)

## -----------------------------------------------------------------------------
boxplot(svout$sv[, 2] ~ true_sv[, 4], 
        xlab = "Dexamethasone Treatment", 
        ylab = "2nd Surrogate Variable")
points(jitter(true_sv[, 4] + 1), svout$sv[, 2])

