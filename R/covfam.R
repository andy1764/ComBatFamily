#' CovBat Family Harmonization
#'
#' Implementation of the CovBat Family of harmonization methods allowing for
#' removal of multivariate batch effects, flexible covariate modeling and
#' alternative estimators for site effect adjustment. Support for modeling of
#' both location and scale via GAMLSS. Additional support for modeling of
#' covariate effects in score location and scale.
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
#' @param score.model Model for scores, defaults to NULL for fitting basic
#'   location and scale model without covariates on the scores
#' @param score.args List of arguments for score model
#' @param eb If \code{TRUE}, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization.
#' @param robust.LS If \code{TRUE}, uses robust location and scale estimators
#'   for error variance and site effect parameters. Uses median and
#'   biweight midvariance
#' @param ref.batch Reference batch, must take value in `levels(bat)`
#' @param percent.var Numeric. The number of harmonized principal component
#'    scores is selected to explain this proportion of the variance
#' @param n.pc Optional numeric. If specified, this number of principal
#'    component scores is harmonized. Overrides \code{percent.var}
#' @param std.var If \code{TRUE}, scales variances to be equal to 1 before PCA.
#' @param ... Additional arguments to `model`
#'
#' @return `covfam` returns a list containing the following components:
#' \item{dat.covbat}{Harmonized data as a matrix with same dimensions as `data`}
#' \item{batch.info}{Batch information, including reference batch if specified}
#' \item{combat.out}{List output of \link[ComBatFamily]{comfam} from the ComBat step}
#' \item{pc.output}{Output of `prcomp` from PCA step}
#' \item{n.pc}{Numeric, number of PCs harmonized}
#' \item{scores.com}{List output of \link[ComBatFamily]{comfam} from the CovBat step}
#'
#' @import stats
#' @export
#'
#' @examples
#' covfam(iris[,1:2], iris$Species)
#' covfam(iris[,1:2], iris$Species, iris[3:4], lm, y ~ Petal.Length + Petal.Width)
covfam <- function(data, bat, covar = NULL, model = lm, formula = NULL,
                   score.model = NULL, score.args = NULL, eb = TRUE,
                   robust.LS = FALSE, ref.batch = NULL, percent.var = 0.95,
                   n.pc = NULL, std.var = TRUE, ...)
{
  n <- nrow(data)
  p <- ncol(data)

  #### Remove mean/variance effects ####
  com_out <- comfam(data, bat, covar, model, formula, eb, robust.LS, ref.batch,
                    ...)
  com_res <- com_out$dat.combat - com_out$estimates$stand.mean

  #### Adjust for multivariate batch effects via PCA ####
  d_pc <- prcomp(com_res, center = TRUE, scale. = std.var)

  # Only adjust PCs specified via percent.var or npc
  if (!is.null(n.pc)) {
    npc <- n.pc
  } else {
    npc <- which(cumsum(d_pc$sdev^2/sum(d_pc$sdev^2)) > percent.var)[1]
  }
  scores <- d_pc$x[,1:npc]

  # ComBat without covariates to remove site effect in score mean/variance
  # If score.model specified, fits that model instead
  if (is.null(score.model)) {
    scores_com <- comfam(scores, bat, eb = FALSE, ref.batch = ref.batch)
  } else {
    scores_com <- do.call(comfam, c(list(scores, bat, covar,
                                         model = score.model, eb = FALSE,
                                         ref.batch = ref.batch), score.args))
  }
  full_scores <- d_pc$x
  full_scores[,1:npc] <- scores_com$dat.combat

  #### Project scores back into observation space ####
  if (std.var) {
    data_covbat <- full_scores %*% t(d_pc$rotation) *
      matrix(d_pc$scale, n, p, byrow = TRUE) +
      matrix(d_pc$center, n, p, byrow = TRUE)
  } else {
    data_covbat <- full_scores %*% t(d_pc$rotation) +
      matrix(d_pc$center, n, p, byrow = TRUE)
  }

  # Reintroduce covariate effects
  data_covbat <- data_covbat + com_out$estimates$stand.mean

  batch_info <- list(
    batch = bat,
    levels = levels(bat)
  )
  batch_info$ref.batch <- ref.batch

  out <- list(dat.covbat = data_covbat, batch.info = batch_info,
              combat.out = com_out, pc.output = d_pc, n.pc = npc,
              scores.combat = scores_com)
  class(out) <- c("covfam")
  out
}
