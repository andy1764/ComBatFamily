
combat.fam <- function(data, covar, bat, model, formula, eb = TRUE,
                       robust.LS = FALSE, debug = TRUE, ...) {
  # Data details and formatting
  n <- nrow(data)
  p <- ncol(data)
  fn <- colnames(data)

  bat <- droplevels(bat)
  batch <- model.matrix(~ -1 + bat)

  # Specify robust location/scale estimators
  if (robust.LS) {
    loc <- median
    scl <- .biweight_midvar
  } else {
    loc <- mean
    scl <- var
  }

  #### Fit specified models ####
  fits <- apply(data, 2, function(y) {
    dat <- data.frame(y = y, covar, I(batch))

    # include batch in formula to target pooled mean/variance
    bat_formula <- update(formula, ~ . + batch + -1)
    do.call(model, list(bat_formula, data = dat, ...))
  })

  #### Standardize the data ####

  # Model matrix for obtaining pooled mean
  mod <- data.frame(covar, I(batch))
  pmod <- mod
  pmod$batch[] <- 1/nlevels(bat)

  stand_mean <- sapply(fits, predict, pmod)
  resid_mean <- sapply(fits, predict, mod)

  var_pooled <- apply(data-resid_mean, 2, scl)*(n-1)/n
  sd_mat <- matrix(sqrt(var_pooled), n, p, byrow = TRUE)

  # TODO: Variance modeling via GAMLSS

  data_stand <- (data-stand_mean)/sd_mat

  #### Obtain location and scale adjustments ####
  gamma_hat <- do.call(rbind, by(data_stand, bat, function(x) apply(x, 2, loc)))
  delta_hat <- do.call(rbind, by(data_stand, bat, function(x) apply(x, 2, scl)))

  # TODO: Empirical Bayes adjustments
  if (eb) {
    gamma_star <- gamma_hat
    delta_star <- delta_hat
  } else {
    gamma_star <- gamma_hat
    delta_star <- delta_hat
  }

  #### Harmonize the data ####

  # Remove batch effects
  batches <- lapply(levels(bat), function(x) which(bat == x))

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
    var.pooled = var_pooled,
    gamma.hat = gamma_hat,
    delta.hat = delta_hat
  )

  debug_out <- NULL
  if (debug) {
    debug_out <- list(
      fits = fits,
      dat.standardized = data_stand
    )
  }

  return(
    list(
      dat.combat = data_combat,
      estimates = estimates,
      debug = debug_out
    )
  )
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
