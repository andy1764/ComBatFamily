#' Utility functions
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom matrixStats rowMeans2 rowVars rowSums2 rowSds
#' @importFrom stats var model.matrix prcomp
#'

# This is a copy of the original code from the standard version of the
# sva package that can be found at
# https://bioconductor.org/packages/release/bioc/html/sva.html
# The original and present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license.

# Author: Jean-Philippe Fortin, fortin946@gmail.com
# This is a modification of the ComBat function code from the sva package that
# can be found at https://bioconductor.org/packages/release/bioc/html/sva.html
# The original and present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license.
# Code optimization improved by Richard Beare
# Modified by Andrew Chen for covbat.R
# Added functionality to use only training data as input and to have
# residualized observations as output
combat_modded <- function(dat, batch, mod=NULL, train = NULL, resid = FALSE,
                          eb = TRUE, parametric = TRUE, mean.only = FALSE,
                          verbose = TRUE)
{
  dat <- as.matrix(dat)

  .checkConstantRows <- function(dat){
    sds <- rowSds(dat)
    ns <- sum(sds==0)
    if (ns>0){
      message <- paste0(ns, " rows (features) were found to be constant across samples. Please remove these rows before running ComBat.")
      stop(message)
    }
  }
  .checkConstantRows(dat)
  if (eb){
    if (verbose) cat("[combat] Performing ComBat with empirical Bayes\n")
  } else {
    if (verbose) cat("[combat] Performing ComBat without empirical Bayes (L/S model)\n")
  }
  # make batch a factor and make a set of indicators for batch
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1+batch)
  if (verbose) cat("[combat] Found",nlevels(batch),'batches\n')

  # A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- lapply(levels(batch), function(x)which(batch==x))
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  #combine batch variable and covariates
  design <- cbind(batchmod,mod)
  # check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[,!check])

  # Number of covariates or covariate levels
  if (verbose) cat("[combat] Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')

  # Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    if(ncol(design)==(n.batch+1)){
      stop("[combat] The covariate is confounded with batch. Remove the covariate and rerun ComBat.")
    }
    if(ncol(design)>(n.batch+1)){
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.')
      } else {
        stop("At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.")
      }
    }
  }

  ## Standardize Data across features
  if (verbose) cat('[combat] Standardizing Data across features\n')

  # Estimate coefficients using training set if specified, otherwise use full data
  if (!is.null(train)) {
    design_tr <- design[train,]

    B.hat1 <- solve(crossprod(design_tr))
    B.hat1 <- tcrossprod(B.hat1, design_tr)
    B.hat <- tcrossprod(B.hat1, dat[,train])
  } else {
    B.hat1 <- solve(crossprod(design))
    B.hat1 <- tcrossprod(B.hat1, design)
    B.hat <- tcrossprod(B.hat1, dat)
  }

  # Standardization Model
  grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
  stand.mean <- crossprod(grand.mean, t(rep(1,n.array)))

  if(!is.null(design)){
    tmp <- design;tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp%*%B.hat)
  }
  s.data <- (dat-stand.mean)/(tcrossprod(sqrt(var.pooled), rep(1,n.array)))

  ## Get regression batch effect parameters
  if (eb){
    if (verbose) cat("[combat] Fitting L/S model and finding priors\n")
  } else {
    if (verbose) cat("[combat] Fitting L/S model\n")
  }
  batch.design <- design[,1:n.batch]
  gamma.hat <- tcrossprod(solve(crossprod(batch.design, batch.design)), batch.design)
  gamma.hat <- tcrossprod(gamma.hat, s.data)
  delta.hat <- NULL
  for (i in batches){
    delta.hat <- rbind(delta.hat,rowVars(s.data[,i], na.rm=TRUE)) # fixed error
  }

  # Empirical Bayes correction:
  gamma.star <- delta.star <- NULL
  gamma.bar <- t2 <- a.prior <- b.prior <- NULL
  if (eb){
    ##Find Priors
    #gamma.bar <- apply(gamma.hat, 1, mean)
    #t2 <- apply(gamma.hat, 1, var)
    gamma.bar <- rowMeans(gamma.hat)
    t2 <- rowVars(gamma.hat)
    a.prior <- apriorMat(delta.hat)
    b.prior <- bpriorMat(delta.hat)

    ##Find EB batch adjustments
    if (parametric){
      if (verbose) cat("[combat] Finding parametric adjustments\n")
      for (i in 1:n.batch){
        temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
        gamma.star <- rbind(gamma.star,temp[1,])
        delta.star <- rbind(delta.star,temp[2,])
      }
    } else {
      if (verbose) cat("[combat] Finding non-parametric adjustments\n")
      for (i in 1:n.batch){
        temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),gamma.hat[i,], delta.hat[i,])
        gamma.star <- rbind(gamma.star,temp[1,])
        delta.star <- rbind(delta.star,temp[2,])
      }
    }

  }

  if (mean.only) {
    delta.star <- array(1, dim = dim(delta.star))
  }

  ### Normalize the Data ###
  if (verbose) cat("[combat] Adjusting the Data\n")
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    if (eb){
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/tcrossprod(sqrt(delta.star[j,]), rep(1,n.batches[j]))
    } else {
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.hat))/tcrossprod(sqrt(delta.hat[j,]), rep(1,n.batches[j]))
    }
    j <- j+1
  }

  if (resid == FALSE) {
    bayesdata <- (bayesdata*(tcrossprod(sqrt(var.pooled), rep(1,n.array))))+stand.mean
  } else {
    bayesdata <- bayesdata*(tcrossprod(sqrt(var.pooled), rep(1,n.array)))
  }

  return(list(dat.combat=bayesdata,
              s.data=s.data,
              gamma.hat=gamma.hat, delta.hat=delta.hat,
              gamma.star=gamma.star, delta.star=delta.star,
              gamma.bar=gamma.bar, t2=t2, a.prior=a.prior, b.prior=b.prior, batch=batch, mod=mod,
              stand.mean=stand.mean, stand.sd=sqrt(var.pooled)[,1])
  )
}

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat){
  m=mean(gamma.hat)
  s2=var(gamma.hat)
  (2*s2+m^2)/s2
}
bprior <- function(gamma.hat){
  m=mean(gamma.hat)
  s2=var(gamma.hat)
  (m*s2+m^3)/s2
}
postmean <- function(g.hat,g.bar,n,d.star,t2){
  (t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)
}
postvar <- function(sum2,n,a,b){
  (.5*sum2+b)/(n/2+a-1)
}

apriorMat <- function(gamma.hat) {
  m <- rowMeans2(gamma.hat)
  s2 <- rowVars(gamma.hat)
  return((2*s2+m^2)/s2)
}
bpriorMat <- function(gamma.hat) {
  m <- rowMeans2(gamma.hat)
  s2 <- rowVars(gamma.hat)
  return((m*s2+m^3)/s2)
}
# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
  #n <- apply(!is.na(sdat),1,sum)
  n <- rowSums(!is.na(sdat))
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  ones <- rep(1,ncol(sdat))

  while(change>conv){
    g.new  <- postmean(g.hat,g.bar,n,d.old,t2)
    #sum2   <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
    #sum2   <- apply((sdat-tcrossprod(g.new, rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
    sum2 <- rowSums2((sdat-tcrossprod(g.new, ones))^2, na.rm=TRUE)
    d.new  <- postvar(sum2,n,a,b)
    change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  #cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}


int.eprior <- function(sdat, g.hat, d.hat){
  g.star <- d.star <- NULL
  r <- nrow(sdat)
  for(i in 1:r){
    g <- g.hat[-i]
    d <- d.hat[-i]
    x <- sdat[i,!is.na(sdat[i,])]
    n <- length(x)
    j <- numeric(n)+1
    dat <- matrix(as.numeric(x), length(g), n, byrow=TRUE)
    resid2 <- (dat-g)^2
    sum2 <- resid2 %*% j
    LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
    LH[LH=="NaN"]=0
    g.star <- c(g.star, sum(g*LH)/sum(LH))
    d.star <- c(d.star, sum(d*LH)/sum(LH))
    ## if(i%%1000==0){cat(i,'\n')}
  }
  adjust <- rbind(g.star,d.star)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}



# Create indices for matched samples across sites, for one covariate.
#batch: site, study or scanner covariate
#x: continuous covariate
#xmin: minimum value for x to be considered
#xmax: maximum value for x to be considered
#step=1: step for the grid; must be a positive integer
createMatchingIndices <- function(x, batch, xmin=NULL, xmax=NULL, step=1){

  stopifnot(length(x)==length(batch))
  batches <- unique(batch)
  n.batches <- length(batches)
  x_per_batch <- split(x, f=batch)[batches]
  if (is.null(xmin)) xmin <- min(x)
  if (is.null(xmax)) xmax <- max(x)
  grid <- seq(xmin,xmax,step)
  n.bins <- length(grid)-1

  # Creating count matrix:
  counts <- matrix(0, n.bins, n.batches)
  for (i in 1:n.bins){
    counts[i,] <- unlist(lapply(x_per_batch, function(temp){
      sum(temp >= grid[i] & temp < grid[i+1])
    }))
  }
  mins <- unlist(apply(counts,1,min)) #Minimal count
  indices <- c()

  # Creating indices:
  for (i in 1:n.bins){
    for (j in 1:n.batches){
      min <- mins[i]
      if (min!=0){
        cand <- which(x >= grid[i] & x < grid[i+1] & batch==batches[j])

        if (length(cand) !=1){
          cand <- sample(cand,min) # Sampling at random
        }


        indices <- c(indices, cand)
      }
    }
  }
  return(indices)
}
