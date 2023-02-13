#' CovBat Family Harmonization
#' Implementation of the CovBat Family of harmonization methods allowing for
#' removal of multivariate batch effects, flexible covariate modeling and
#' alternative estimators for site effect adjustment. Support for modeling of
#' both location and scale via GAMLSS. Additional support for modeling of
#' covariate effects in score location and scale.
#'
#' @param data \emph{n x p} data frame or matrix of observations where
#'   \emph{p} is the number of features and \emph{n} is the number of subjects.
#' @param covar Data frame or matrix of covariates supplied to `model`
#' @param bat Factor indicating batch (often equivalent to site or scanner)
#' @param model Model function. ComBat Family supports any models that take
#'   arguments `formula` and `data`, but are limited to models fitting with
#'   identity link (e.g. `family = gaussian(link = "identity")`). This includes
#'   \link[stats]{lm}, \link[mgcv]{gam}, \link[gamlss]{gamlss},
#'   \link[quantreg]{rq}, \link[lme4]{lmer}, and more
#' @param formula Formula for `model`, format is dependent on choice of model.
#' @param score.model Model for scores, defaults to NULL for fitting basic
#'   location and scale model without covariates on the scores
#' @param score.args List of arguments for score model, requires `formula`
#' @param eb If \code{TRUE}, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization.
#' @param robust.LS If \code{TRUE}, uses robust location and scale estimators
#'   for error variance and site effect parameters. Currently uses median and
#'   biweight midvariance
#' @param percent.var Numeric. The number of harmonized principal component
#'    scores is selected to explain this proportion of the variance.
#' @param n.pc Optional numeric. If specified, this number of principal
#'    component scores is harmonized. Overrides \code{percent.var}.
#' @param std.var If \code{TRUE}, scales variances to be equal to 1 before PCA.
#' @param debug Whether to output model fits and intermediate data frames
#' @param ... Additional arguments to `model`
#'
#' @return
#' @export
#'
#' @examples
covfam <- function(data, covar, bat, model, formula,
                   score.model = NULL, score.args = NULL, eb = TRUE,
                   robust.LS = FALSE, percent.var = 0.95, n.pc = NULL,
                   std.var = TRUE, debug = FALSE,
                   ...)
{
  # Use ComBat to remove mean/variance effects
  com_out <- combat.fam(data, covar, bat, model, formula, eb,
                        robust.LS, debug = TRUE, ...)
  x_com <- com_out$debug$data.resid

  # PC on ComBat-adjusted data
  x_pc <- prcomp(x_com, center = TRUE, scale. = std.var)

  # Subset scores based on percent of variance explained
  if (!is.null(n.pc)) {
    npc <- n.pc
  } else {
    npc <- which(cumsum(x_pc$sdev^2/sum(x_pc$sdev^2)) > percent.var)[1]
  }
  scores <- x_pc$x[,1:npc]

  # ComBat without covariates to remove site effect in score mean/variance
  # If score.model specified, fits that model instead
  if (is.null(score.model)) {
    scores_com <- combat.fam(scores, covar, bat, lm, y ~ 1, eb = FALSE,
                             robust.LS = FALSE)
  } else {
    scores_com <- do.call(combat.fam,
                          c(list(scores, covar, bat, model = score.model,
                                 eb = FALSE, robust.LS = FALSE, debug = TRUE),
                            score.args))
  }

  full_scores <- x_pc$x
  full_scores[,1:npc] <- scores_com$dat.combat

  # Project scores back into observation space
  if (std.var) {
    data_covbat <- full_scores %*% t(x_pc$rotation) *
      matrix(x_pc$scale, dim(x_com)[1], dim(x_com)[2], byrow = TRUE) +
      matrix(x_pc$center, dim(x_com)[1], dim(x_com)[2], byrow = TRUE)
  } else {
    data_covbat <- full_scores %*% t(x_pc$rotation) +
      matrix(x_pc$center, dim(x_com)[1], dim(x_com)[2], byrow = TRUE)
  }

  data_covbat <- data_covbat + com_out$estimates$stand.mean

  debug_out <- NULL
  if (debug) {
    debug_out <- list(
      scores.combat = scores_com
    )
  }

  return(
    list(
      dat.covbat = data_covbat,
      combat.out = com_out,
      debug = debug_out
    )
  )
}
