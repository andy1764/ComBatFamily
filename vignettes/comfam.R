## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(neuroCombat)
  library(ComBatFamily)
})

# generate toy dataset
set.seed(8888)
n <- 20
p <- 100
bat <- as.factor(c(rep("a", n/2), rep("b", n/2)))
q <- 2
covar <- matrix(rnorm(n*q), n, q)
colnames(covar) <- paste0("x", 1:q)
data <- data.frame(matrix(rnorm(n*p), n, p))

cf <- comfam(data, covar, bat, lm, formula = y ~ x1 + x2)
c <- neuroCombat(t(data), bat, covar, verbose = FALSE)

max(cf$dat.combat - t(c$dat.combat))

## -----------------------------------------------------------------------------
suppressPackageStartupMessages(
  library(mgcv)
)
cg <- comfam(data, covar, bat, gam, formula = y ~ s(x1) + x2)

## -----------------------------------------------------------------------------
# devtools::install_github("jcbeer/longCombat")
suppressPackageStartupMessages({
  library(longCombat)
  library(lme4)
})

# generate toy dataset
n <- 20
p <- 10
r <- 2 # repeats
bat <- as.factor(rep(c(rep("a", n/2), rep("b", n/2)), r))
q <- 2
covar <- matrix(rnorm(n*r*q), n*r, q)
colnames(covar) <- paste0("x", 1:q)
covar <- cbind(covar, ID = rep(1:n, r))
data <- data.frame(matrix(rnorm(n*r*p), n*r, p))

cf <- comfam(data, covar, bat, lmer, y ~ x1 + x2 + (1 | ID))
c <- longCombat("ID", "x1", "bat", 1:p, "x1 + x2", "(1 | ID)",
                data = cbind(data, covar, bat), method = "MSR", verbose = FALSE)

max(cf$dat.combat - c$data_combat[,-(1:3)])

## -----------------------------------------------------------------------------
suppressPackageStartupMessages(
  library(gamlss)
)

n <- 20
p <- 10
bat <- as.factor(c(rep("a", n/2), rep("b", n/2)))
q <- 2
covar <- matrix(rnorm(n*q), n, q)
colnames(covar) <- paste0("x", 1:q)
data <- data.frame(matrix(rnorm(n*p), n, p))

cgl <- comfam(data, covar, bat, gamlss, y ~ x1 + x2,
              sigma.formula = ~ x1 + x2, 
              control = gamlss.control(trace = FALSE))

