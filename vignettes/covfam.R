## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(CovBat)
  library(ComBatFamily)
})

# generate toy dataset
set.seed(8888)
n <- 20
p <- 5
bat <- as.factor(c(rep("a", n/2), rep("b", n/2)))
q <- 2
covar <- matrix(rnorm(n*q), n, q)
colnames(covar) <- paste0("x", 1:q)
data <- data.frame(matrix(rnorm(n*p), n, p))

cf <- covfam(data, bat, covar, lm, formula = y ~ x1 + x2)
c <- covbat(t(data), bat, covar)

max(cf$dat.covbat - t(c$dat.covbat))

## -----------------------------------------------------------------------------
# no covariates specified
cf <- covfam(data, bat)
c <- covbat(t(data), bat)

max(cf$dat.covbat - t(c$dat.covbat))

## -----------------------------------------------------------------------------
suppressPackageStartupMessages(
  library(mgcv)
)
cg <- covfam(data, bat, covar, gam, formula = y ~ s(x1) + x2)

## -----------------------------------------------------------------------------
# no covariates specified
cf <- covfam(data, bat, covar, formula = y ~ x1 + x2, score.model = lm,
             score.args = list(formula = y ~ x1 + x2))

