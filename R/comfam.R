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
#' @param debug Whether to output model fits and intermediate data frames
#' @param ... Additional arguments to `model`
#'
#' @return `comfam` returns a list containing the following components:
#' \item{dat.combat}{Harmonized data as a matrix with same dimensions as `data`}
#' \item{estimates}{List of estimates from standardization and batch effect correction}
#' @export
#'
#' @examples
#' comfam(iris[,1:2], iris$Species)
#' comfam(iris, iris$Species, iris[3:4], lm, y ~ Petal.Length + Petal.Width)
comfam <- function(data, bat, covar = NULL, model = lm, formula = NULL, eb = TRUE,
                   robust.LS = FALSE, debug = FALSE, ...) {
  # Data details and formatting
  n <- nrow(data)
  p <- ncol(data)
  fn <- colnames(data)

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

  stand_mean <- sapply(fits, predict, newdata = pmod)
  resid_mean <- sapply(fits, predict, newdata = mod)

  var_pooled <- apply(data - resid_mean, 2, scl)*(n-1)/n

  if (hasArg(sigma.formula)) {
    sd_mat <- sapply(fits, predict, newdata = pmod, what = "sigma")
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
      g_old  <- gamma_hat[i,]
      d_old  <- delta_hat[i,]
      change_old <- 1
      change <- 1
      count  <- 0
      while(change > 10e-5){
        g_new <- (n_b*g_var*gamma_hat[i,] + d_old*g_bar)/(n_b*g_var + d_old)

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

  out <- list(dat.combat = data_combat, estimates = estimates)
  if (debug) {
    out$debug <- list(fits = fits, dat.standardized = data_stand)
  }
  out
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
