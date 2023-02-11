# test script for converting ComBat to a sensible format

# make toy data
n <- 20
p <- 100
bat <- as.factor(c(rep("a", n/2), rep("b", n/2)))

q <- 2
covar <- matrix(rnorm(n*q), n, q)

data <- matrix(rnorm(n*p), n, p)

# original ComBat jank
dat <- t(dat)

batch <- bat
n.batch <- nlevels(batch)
batches <- lapply(levels(batch), function(x)which(batch==x))
n.batches <- sapply(batches, length)
n.array <- sum(n.batches)

design <- cbind(model.matrix(~ -1 + bat), mod)

B.hat1 <- solve(crossprod(design))
B.hat1 <- tcrossprod(B.hat1, design)
B.hat <- tcrossprod(B.hat1, dat)

grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
stand.mean <- crossprod(grand.mean, t(rep(1,n.array)))

tmp <- design;tmp[,c(1:n.batch)] <- 0
stand.mean <- stand.mean+t(tmp%*%B.hat)



# what the heck is happening with stand.mean
predx1 <- design
predx1[,1:2] <- 0.5
predx1[,-(1:2)] <- 0

grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
stand.mean <- crossprod(grand.mean, t(rep(1,n.array)))

predx1 %*% B.hat - t(stand.mean)

predx2 <- design
predx2[,1:2] <- 0
predx2 %*% B.hat - tmp%*%B.hat

predx <- design
predx[,1:2] <- 0.5
predx %*% B.hat - t(stand.mean) - tmp%*%B.hat

# okay come on
fits <- lapply(1:nrow(dat), function(x) {
  y <- dat[x,]
  lm(y ~ design - 1)
})

px <- fits[[1]]$model
px$design[,1:2] <- 0.5

new.stand.mean <- sapply(fits, function(f) {
  predict(f, px)
})

new.stand.mean - predx %*% B.hat


# better but jank
batm <- model.matrix(~ -1 + bat)
batmod <- cbind(batm, mod)

fits <- lapply(1:nrow(dat), function(x) {
  y <- dat[x,]
  lm(y ~ batmod - 1)
})

predx <- batmod
predx[,1:2] <- 0.5
new.stand.mean <- t(sapply(fits, function(f) {
  predict(f, data.frame(predx))
}))

predx <- batmod
predx[,1:2] <- 0.5
predx[,-(1:2)] <- 0
new.grand.mean <- t(sapply(fits, function(f) {
  predict(f, data.frame(predx))
}))

fitted <- t(sapply(fits, predict))

var.pooled <- apply(dat-fitted, 1, var)*(n-1)/n

new.stand.mean - stand.mean
var.pooled - ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)


##### ComBat Family vs OG ComBat ####
n <- 20
p <- 100
bat <- as.factor(c(rep("a", n/2), rep("b", n/2)))
q <- 2
covar <- matrix(rnorm(n*q), n, q)
colnames(covar) <- paste0("x", 1:q)
data <- data.frame(matrix(rnorm(n*p), n, p))

cf <- combat.fam(data, covar, bat, lm, formula = y ~ x1 + x2)

library(neuroCombat)
c <- neuroCombat(t(data), bat, covar, eb = TRUE)

cf$estimates$stand.mean - t(c$estimates$stand.mean)
cf$estimates$gamma.hat - c$estimates$gamma.hat
cf$estimates$delta.hat - c$estimates$delta.hat

cf$dat.combat - t(c$dat.combat)


##### ComBat Family GAMLSS ####
n <- 20
p <- 100
bat <- as.factor(c(rep("a", n/2), rep("b", n/2)))
q <- 2
covar <- matrix(rnorm(n*q), n, q)
colnames(covar) <- paste0("x", 1:q)
data <- data.frame(matrix(rnorm(n*p), n, p))

library(gamlss)
cf <- combat.fam(data, covar, bat, gamlss, y ~ x1 + x2,
                 sigma.formula = ~ x1 + x2)



##### ComBat Family Longitudinal? ####
library(lme4)

n <- 20
p <- 100
r <- 2 # repeats
bat <- as.factor(rep(c(rep("a", n/2), rep("b", n/2)), r))
q <- 2
covar <- matrix(rnorm(n*r*q), n*r, q)
colnames(covar) <- paste0("x", 1:q)
covar <- cbind(covar, ID = rep(1:n, r))
data <- data.frame(matrix(rnorm(n*r*p), n*r, p))

cf <- combat.fam(data, covar, bat, lmer, y ~ x1 + x2 + (1 | ID))

# technically works but prob not correct


