#' ComBat Family Harmonization
#'
#' Implementation of the ComBat Family of harmonization methods allowing for
#' flexible covariate modeling and alternative estimators for site effect
#' adjustment. Support for modeling of both location and scale via GAMLSS and
#' longitudinal harmonization via mixed effects models.
#'
#' @param data \emph{n x p} data frame or matrix of observations where
#'   \emph{p} is the number of features and \emph{n} is the number of subjects.
#' @param bat Factor indicating batch (often equivalent to site or scanner)
#' @param covar Data frame or matrix of covariates supplied to `model`
#' @param model Model function. ComBat Family supports any models that take
#'   arguments `formula` and `data`, but are limited to models fitting with
#'   identity link (e.g. `family = gaussian(link = "identity")`). This includes
#'   \link[stats]{lm}, \link[mgcv]{gam}, \link[gamlss]{gamlss},
#'   \link[quantreg]{rq}, \link[lme4]{lmer}, and more
#' @param formula Formula for `model` starting with `y ~` where `y` represents
#'   each feature
#' @param eb If \code{TRUE}, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization
#' @param robust.LS If \code{TRUE}, uses robust location and scale estimators
#'   for error variance and site effect parameters. Currently uses median and
#'   biweight midvariance
#' @param ref.batch Reference batch, must take value in `levels(bat)`
#' @param ... Additional arguments to `model`
#'
#' @return `comfam` returns a list containing the following components:
#' \item{dat.combat}{Harmonized data as a matrix with same dimensions as `data`}
#' \item{batch.info}{Batch information, including reference batch if specified}
#' \item{fits}{List of model fits from regression step, outputs of `model` for each feature}
#' \item{estimates}{List of estimates from standardization and batch effect correction}
#' @export
#'
#' @seealso
#' \link[ComBatFamily]{plot.comfam} for assessing regression fit via
#' diagnostic plots associated with `model`
#' \link[ComBatFamily]{predict.comfam} for applying ComBat parameters for
#' harmonization of new observations
#'
#' @examples
#' comfam(iris[,1:2], iris$Species)
#' comfam(iris[,1:2], iris$Species, iris[3:4], lm, y ~ Petal.Length + Petal.Width)
comfam <- function(data, bat, covar = NULL, model = lm, formula = NULL,
                   eb = TRUE, robust.LS = FALSE, ref.batch = NULL, ...) {
  if (hasArg(family)) {
    if (!(list(...)$family$family %in% c("gaussian", "Normal"))) {
      warning("Families other than Gaussian are supported but experimental, output dataset will not necessarily be in the original space.")

      warning("EB step will still assume Gaussian errors.")
    }
  }

  # Data details and formatting
  n <- nrow(data)
  p <- ncol(data)

  bat <- droplevels(bat)
  batch <- model.matrix(~ -1 + bat)
  batches <- lapply(levels(bat), function(x) which(bat == x))
  n_batches <- sapply(batches, length)

  # Specify robust location/scale estimators
  if (robust.LS) {
    loc <- median
    scl <- .biweight_midvar
  } else {
    loc <- mean
    scl <- var
  }

  #### Fit specified models ####
  if (is.null(covar)) {
    mod <- data.frame(I(batch))
  } else {
    mod <- data.frame(covar, I(batch))
  }

  if (is.null(covar) | is.null(formula)) {
    formula <- y ~ 1
  }

  fits <- apply(data, 2, function(y) {
    dat <- data.frame(y = y, mod)

    # include batch in formula to target pooled mean/variance
    bat_formula <- update(formula, ~ . + batch + -1)
    do.call(model, list(formula = bat_formula, data = dat, ...))
  })

  #### Standardize the data ####
  # Model matrix for obtaining pooled mean
  pmod <- mod
  pmod$batch[] <- matrix(n_batches/n, n, nlevels(bat), byrow = TRUE)

  # Reference batch
  if (!is.null(ref.batch)) {
    pmod$batch[] <- 0
    pmod$batch[,which(levels(bat) == ref.batch)] <- 1
  }

  stand_mean <- sapply(fits, predict, newdata = pmod, type = "response")
  resid_mean <- sapply(fits, predict, newdata = mod, type = "response")

  var_pooled <- apply(data - resid_mean, 2, scl)*(n-1)/n

  if (hasArg(sigma.formula)) {
    sd_mat <- sapply(fits, predict, newdata = pmod, what = "sigma",
                     type = "response")
  } else {
    sd_mat <- matrix(sqrt(var_pooled), n, p, byrow = TRUE)
  }

  data_stand <- (data-stand_mean)/sd_mat

  #### Obtain location and scale adjustments ####
  gamma_hat <- do.call(rbind, by(data_stand, bat, function(x) apply(x, 2, loc)))
  delta_hat <- do.call(rbind, by(data_stand, bat, function(x) apply(x, 2, scl)))

  # Empirical Bayes adjustments
  if (eb) {
    gamma_star <- NULL
    delta_star <- NULL

    for (i in 1:nlevels(bat)) {
      n_b <- n_batches[i]

      # method of moments estimates
      g_bar <- mean(gamma_hat[i,])
      g_var <- var(gamma_hat[i,])

      d_bar <- mean(delta_hat[i,])
      d_var <- var(delta_hat[i,])

      d_a <- (2 * d_var + d_bar^2)/d_var
      d_b <- (d_bar * d_var + d_bar^3)/d_var

      # adjust within batch
      bdat <- data_stand[batches[[i]],]
      g_orig <- gamma_hat[i,]
      g_old  <- gamma_hat[i,]
      d_old  <- delta_hat[i,]

      change_old <- 1
      change <- 1
      count  <- 0
      while(change > 10e-5){
        g_new <- (n_b*g_var*g_orig + d_old*g_bar)/(n_b*g_var + d_old)

        if (robust.LS) {
          sum2 <- (n_b-1) * sapply(1:p, function(v) {
            .biweight_midvar(bdat[,v], g_new[v])})
        } else {
          sum2   <- colSums(sweep(bdat, 2, g_new)^2)
        }

        d_new <- (sum2/2 + d_b)/(n_b/2 + d_a - 1)

        change <- max(abs(g_new - g_old)/g_old, abs(d_new - d_old)/d_old)

        if (count > 30) {
          if (change > change_old) {
            warning("Empirical Bayes step failed to converge after 30 iterations,
    	            using estimate before change between iterations increases.")
            break
          }
        }

        g_old <- g_new
        d_old <- d_new

        change_old <- change
        count <- count+1
      }

      gamma_star <- rbind(gamma_star, g_new)
      delta_star <- rbind(delta_star, d_new)
    }

    rownames(gamma_star) <- rownames(gamma_hat)
    rownames(delta_star) <- rownames(delta_hat)
  } else {
    gamma_star <- gamma_hat
    delta_star <- delta_hat
  }

  #### Harmonize the data ####
  # Remove batch effects
  data_nb <- data_stand
  for (i in 1:nlevels(bat)) {
    data_nb[batches[[i]],] <- sweep(data_nb[batches[[i]],], 2,
                                    gamma_star[i,], "-")
    data_nb[batches[[i]],] <- sweep(data_nb[batches[[i]],], 2,
                                    sqrt(delta_star[i,]), "/")
  }

  # Reintroduce covariate effects
  data_combat <- data_nb*sd_mat + stand_mean

  estimates <-  list(
    stand.mean = stand_mean,
    stand.sd = sd_mat,
    var.pooled = var_pooled,
    gamma.hat = gamma_hat,
    delta.hat = delta_hat,
    gamma.star = gamma_star,
    delta.star = delta_star
  )

  batch_info <- list(
    batch = bat,
    batch.mod = pmod$batch,
    ref.batch = ref.batch
  )

  out <- list(dat.combat = data_combat, batch.info = batch_info,
              fits = fits, estimates = estimates)
  class(out) <- "comfam"
  out
}

#' Apply Harmonization to New Data
#'
#' Using parameters estimated via `comfam`, apply harmonization on new data.
#' `predict.comfam` will estimate new batch adjustments if new batches are
#' specified. For batches with existing estimates, the estimates from `object`
#' are used. Harmonization targets are the same as `object` (e.g. `ref.batch`
#' from `object` if specified).
#'
#' @param object Object of class `comfam`, typically output of
#'   \link[ComBatFamily]{comfam}
#' @param newdata \emph{n x p} data frame or matrix of new observations where
#'   \emph{p} is the number of features and \emph{n} is the number of subjects.
#'   The features must match the original `data` used in `object`
#' @param newbat Factor indicating new batch (often equivalent to site or scanner)
#' @param newcovar Data frame or matrix of new covariates supplied to `model`.
#'   Must contain all variables specified in the original `formula` used in
#'   `object`.
#' @param eb If \code{TRUE}, uses ComBat model with empirical Bayes for new batches
#' @param robust.LS If \code{TRUE}, uses robust location and scale estimators
#'   for new batch effect estimates Currently uses median and biweight
#'   midvariance
#'
#' @return `comfam` returns a list containing the following components:
#' \item{dat.combat}{New harmonized data as a matrix with same dimensions as `newdata`}
#' \item{batch.info}{New batch information, including reference batch if specified}
#' \item{fits}{List of model fits from regression step, forwarded from `object`}
#' \item{estimates}{List of estimates from standardization and batch effect correction, including new batches if relevant}
#' @export
#'
#' @examples
#' com_out <- comfam(iris[1:75,1:2], iris$Species[1:75])
#'
#' # out-of-sample with new batch
#' out_pred <- predict(com_out, iris[76:150,1:2], iris$Species[76:150])
#'
#' # in-sample
#' in_pred <- predict(com_out, iris[1:25,1:2], iris$Species[1:25])
#' max(in_pred$dat.combat - com_out$dat.combat[1:25,])
predict.comfam <- function(object, newdata, newbat, newcovar = NULL,
                           robust.LS = FALSE, eb = TRUE) {
  data <- newdata
  n <- nrow(data)
  p <- ncol(data)

  # Specify robust location/scale estimators
  if (robust.LS) {
    loc <- median
    scl <- .biweight_midvar
  } else {
    loc <- mean
    scl <- var
  }

  bat <- object$batch.info$batch
  batch_mod <- object$batch.info$batch.mod
  stand_mean <- object$estimates$stand.mean
  var_pooled <- object$estimates$var.pooled
  gamma_hat <- object$estimates$gamma.hat
  delta_hat <- object$estimates$delta.hat
  gamma_star <- object$estimates$gamma.star
  delta_star <- object$estimates$delta.star
  fits <- object$fits

  #### Match new batches to old batches ####
  bat <- droplevels(bat)
  newbat <- droplevels(newbat)
  bat_levels <- union(levels(bat), levels(newbat))
  newbat_app <- factor(newbat, bat_levels)
  batches <- lapply(levels(newbat_app), function(x) which(newbat_app == x))

  # new batches to estimate/adjust
  newbat_est <- which(bat_levels %in% setdiff(levels(newbat), levels(bat)))
  newbat_adj <- which(bat_levels %in% levels(newbat))

  #### Standardize the data ####
  # resize batch_mod
  batch <- matrix(batch_mod[1,], n, nlevels(bat), byrow = TRUE)

  if (is.null(newcovar)) {
    pmod <- data.frame(I(batch))
  } else {
    pmod <- data.frame(newcovar, I(batch))
  }

  stand_mean <- sapply(fits, predict, newdata = pmod, type = "response")
  if (hasArg(sigma.formula)) {
    sd_mat <- sapply(fits, predict, newdata = pmod, what = "sigma",
                     type = "response")
  } else {
    sd_mat <- matrix(sqrt(var_pooled), n, p, byrow = TRUE)
  }

  data_stand <- (data - stand_mean)/sd_mat

  #### Obtain location and scale adjustments ####
  # get naive estimates for new batches
  for (i in newbat_est) {
    gamma_hat <- rbind(gamma_hat, apply(data_stand[batches[[i]],], 2, loc))
    delta_hat <- rbind(delta_hat, apply(data_stand[batches[[i]],], 2, scl))
  }

  rownames(gamma_hat) <- rownames(delta_hat) <- bat_levels

  # Empirical Bayes adjustments for new batches
  if (eb) {
    for (i in newbat_est) {
      n_b <- length(batches[[i]])

      # method of moments estimates
      g_bar <- mean(gamma_hat[i,])
      g_var <- var(gamma_hat[i,])
      d_bar <- mean(delta_hat[i,])
      d_var <- var(delta_hat[i,])
      d_a <- (2 * d_var + d_bar^2)/d_var
      d_b <- (d_bar * d_var + d_bar^3)/d_var

      # adjust within batch
      bdat <- data_stand[batches[[i]],]
      g_orig <- gamma_hat[i,]
      g_old  <- gamma_hat[i,]
      d_old  <- delta_hat[i,]

      change_old <- 1
      change <- 1
      count  <- 0
      while(change > 10e-5){
        g_new <- (n_b*g_var*g_orig + d_old*g_bar)/(n_b*g_var + d_old)

        if (robust.LS) {
          sum2 <- (n_b-1) * sapply(1:p, function(v) {
            .biweight_midvar(bdat[,v], g_new[v])})
        } else {
          sum2   <- colSums(sweep(bdat, 2, g_new)^2)
        }

        d_new <- (sum2/2 + d_b)/(n_b/2 + d_a - 1)

        change <- max(abs(g_new - g_old)/g_old, abs(d_new - d_old)/d_old)

        if (count > 30) {
          if (change > change_old) {
            warning("Empirical Bayes step failed to converge after 30 iterations,
    	            using estimate before change between iterations increases.")
            break
          }
        }

        g_old <- g_new
        d_old <- d_new

        change_old <- change
        count <- count+1
      }

      gamma_star <- rbind(gamma_star, g_new)
      delta_star <- rbind(delta_star, d_new)
    }

    rownames(gamma_star) <- rownames(gamma_hat)
    rownames(delta_star) <- rownames(delta_hat)
  } else {
    gamma_star <- gamma_hat
    delta_star <- delta_hat
  }

  #### Harmonize the data ####
  # Remove batch effects
  data_nb <- data_stand
  for (i in newbat_adj) {
    data_nb[batches[[i]],] <- sweep(data_nb[batches[[i]],], 2,
                                    gamma_star[i,], "-")
    data_nb[batches[[i]],] <- sweep(data_nb[batches[[i]],], 2,
                                    sqrt(delta_star[i,]), "/")
  }

  # Reintroduce covariate effects
  data_combat <- data_nb*sd_mat + stand_mean

  estimates <-  list(
    stand.mean = stand_mean,
    stand.sd = sd_mat,
    var.pooled = var_pooled,
    gamma.hat = gamma_hat,
    delta.hat = delta_hat,
    gamma.star = gamma_star,
    delta.star = delta_star
  )

  batch_info <- list(
    batch = newbat_app,
    batch.mod = batch,
    ref.batch = object$batch.info$ref.batch
  )

  out <- list(dat.combat = data_combat, batch.info = batch_info,
              fits = fits, estimates = estimates)
  class(out) <- "comfam"
  out
}

#' Plot Diagnostics for `comfam`
#'
#' Diagnostic plots for original model fits in `comfam`, leverages S3 plot
#' methods for `model` (e.g. \link[stats]{plot.lm})
#'
#' @param object Object of class `comfam`, typically output of
#'   \link[ComBatFamily]{comfam}
#' @param feature Feature to diagnose, either index or variable name
#'
#' @return
#' @export
#'
#' @examples
#' com_out <- comfam(iris[,1:2], iris$Species)
#' plot(com_out)
plot.comfam <- function(object, feature) {
  plot(object$fits[[feature]])
}

.biweight_midvar <- function(data, center=NULL, norm.unbiased = TRUE) {
  if (is.null(center)) {
    center <- median(data)
  }

  mad <- median(abs(data - center))
  d <- data - center
  c <- ifelse(norm.unbiased, 9/qnorm(0.75), 9)
  u <- d/(c*mad)

  n <- length(data)
  indic <- abs(u) < 1

  num <- sum(indic * d^2 * (1 - u^2)^4)
  dem <- sum(indic * (1 - u^2) * (1 - 5*u^2))^2

  n * num/dem
}
