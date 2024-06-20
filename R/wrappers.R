#' ComBat Harmonization
#'
#' Implementation of ComBat (Johnson et al., 2007) using \link[ComBatFamily]{comfam} with `model` as \link[stats]{lm}.
#'
#' @param data \emph{n x p} data frame or matrix of observations where
#'   \emph{p} is the number of features and \emph{n} is the number of subjects
#' @param bat Factor indicating batch (often equivalent to site or scanner)
#' @param covar Data frame or matrix of covariates supplied to \link[stats]{lm}
#' @param formula Formula for \link[stats]{lm} starting with `y ~` where `y` represents
#'   each feature
#' @param eb If \code{TRUE}, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization
#' @param robust.LS If \code{TRUE}, uses robust location and scale estimators
#'   for error variance and site effect parameters. Currently uses median and
#'   biweight midvariance
#' @param ref.batch Reference batch, must take value in `levels(bat)`
#' @param ... Additional arguments passed to \link[stats]{lm}
#'
#' @return `combat` returns a list containing the following components:
#' \item{dat.combat}{Harmonized data as a matrix with same dimensions as `data`}
#' \item{batch.info}{Batch information, including reference batch if specified}
#' \item{fits}{List of model fits from regression step, outputs of \link[stats]{lm} for each feature}
#' \item{estimates}{List of estimates from standardization and batch effect correction}
#'
#' @import stats
#' @importFrom methods hasArg
#' @export
#'
#' @seealso
#' \link[ComBatFamily]{plot.comfam} for assessing regression fit via
#' diagnostic plots associated with \link[stats]{lm}
#'
#' \link[ComBatFamily]{predict.comfam} for applying ComBat parameters for
#' harmonization of new observations
#'
#' @examples
#' combat(iris[,1:2], iris$Species)
#' combat(iris[,1:2], iris$Species, iris[3:4], y ~ Petal.Length + Petal.Width)
#'
#' @references
#' Fortin, J.-P., Cullen, N., Sheline, Y. I., Taylor, W. D., Aselcioglu, I., Cook, P. A., Adams, P., Cooper, C., Fava, M., McGrath, P. J., McInnis, M., Phillips, M. L., Trivedi, M. H., Weissman, M. M., & Shinohara, R. T. (2018). Harmonization of cortical thickness measurements across scanners and sites. *NeuroImage*, 167, 104–120. https://doi.org/10.1016/j.neuroimage.2017.11.024
#'
#' Fortin, J.-P., Parker, D., Tunç, B., Watanabe, T., Elliott, M. A., Ruparel, K., Roalf, D. R., Satterthwaite, T. D., Gur, R. C., Gur, R. E., Schultz, R. T., Verma, R., & Shinohara, R. T. (2017). Harmonization of multi-site diffusion tensor imaging data. *NeuroImage*, 161, 149–170. https://doi.org/10.1016/j.neuroimage.2017.08.047
#'
#' Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8(1), 118–127. https://doi.org/10.1093/biostatistics/kxj037
combat <- function(data, bat, covar = NULL, formula = NULL,
                   eb = TRUE, robust.LS = FALSE, ref.batch = NULL, ...) {
  comfam(data, bat, covar, lm, formula, eb, robust.LS, ref.batch, ...)
}

#' ComBat-GAM harmonization
#'
#' Implementation of ComBat-GAM (Pomponio et al., 2020) using \link[ComBatFamily]{comfam} with `model` as \link[mgcv]{gam}.
#'
#' @param data \emph{n x p} data frame or matrix of observations where
#'   \emph{p} is the number of features and \emph{n} is the number of subjects
#' @param bat Factor indicating batch (often equivalent to site or scanner)
#' @param covar Data frame or matrix of covariates supplied to \link[mgcv]{gam}
#' @param formula Formula for \link[mgcv]{gam} starting with `y ~` where `y` represents
#'   each feature
#' @param eb If \code{TRUE}, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization
#' @param robust.LS If \code{TRUE}, uses robust location and scale estimators
#'   for error variance and site effect parameters. Currently uses median and
#'   biweight midvariance
#' @param ref.batch Reference batch, must take value in `levels(bat)`
#' @param ... Additional arguments passed to \link[mgcv]{gam}
#'
#' @return `combat` returns a list containing the following components:
#' \item{dat.combat}{Harmonized data as a matrix with same dimensions as `data`}
#' \item{batch.info}{Batch information, including reference batch if specified}
#' \item{fits}{List of model fits from regression step, outputs of \link[mgcv]{gam} for each feature}
#' \item{estimates}{List of estimates from standardization and batch effect correction}
#'
#' @import stats mgcv
#' @importFrom methods hasArg
#' @export
#'
#' @seealso
#' \link[ComBatFamily]{plot.comfam} for assessing regression fit via
#' diagnostic plots associated with \link[mgcv]{gam}
#'
#' \link[ComBatFamily]{predict.comfam} for applying ComBat parameters for
#' harmonization of new observations
#'
#' @examples
#' combat_gam(iris[,1:2], iris$Species)
#' combat_gam(iris[,1:2], iris$Species, iris[3:4], y ~ s(Petal.Length) + Petal.Width)
#'
#' @references
#' Pomponio, R., Erus, G., Habes, M., Doshi, J., Srinivasan, D., Mamourian, E., Bashyam, V., Nasrallah, I. M., Satterthwaite, T. D., Fan, Y., Launer, L. J., Masters, C. L., Maruff, P., Zhuo, C., Völzke, H., Johnson, S. C., Fripp, J., Koutsouleris, N., Wolf, D. H., … Shou, H., Davatzikos, C. (2020). Harmonization of large MRI datasets for the analysis of brain imaging patterns throughout the lifespan. *NeuroImage*, 208, 116450. https://doi.org/10.1016/j.neuroimage.2019.116450
combat_gam <- function(data, bat, covar = NULL, formula = NULL,
                       eb = TRUE, robust.LS = FALSE, ref.batch = NULL, ...) {
  comfam(data, bat, covar, gam, formula, eb, robust.LS, ref.batch, ...)
}

#' Longitudinal ComBat harmonization
#'
#' Implementation of Longitudinal ComBat (Beer et al., 2020) using \link[ComBatFamily]{comfam} with `model` as \link[lme4]{lmer}.
#' Currently, this implementation is lacking the REML option for variance estimates.
#' Please use the longComBat package (https://github.com/jcbeer/longCombat) if this is needed.
#'
#' @param data \emph{n x p} data frame or matrix of observations where
#'   \emph{p} is the number of features and \emph{n} is the number of subjects.
#' @param bat Factor indicating batch (often equivalent to site or scanner)
#' @param covar Data frame or matrix of covariates supplied to \link[lme4]{lmer}
#' @param formula Formula for \link[lme4]{lmer} starting with `y ~` where `y` represents
#'   each feature
#' @param eb If \code{TRUE}, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization
#' @param robust.LS If \code{TRUE}, uses robust location and scale estimators
#'   for error variance and site effect parameters. Currently uses median and
#'   biweight midvariance
#' @param ref.batch Reference batch, must take value in `levels(bat)`
#' @param ... Additional arguments passed to \link[lme4]{lmer}
#'
#' @return `combat` returns a list containing the following components:
#' \item{dat.combat}{Harmonized data as a matrix with same dimensions as `data`}
#' \item{batch.info}{Batch information, including reference batch if specified}
#' \item{fits}{List of model fits from regression step, outputs of \link[lme4]{lmer} for each feature}
#' \item{estimates}{List of estimates from standardization and batch effect correction}
#'
#' @import stats lme4
#' @importFrom methods hasArg
#' @export
#'
#' @seealso
#' \link[ComBatFamily]{plot.comfam} for assessing regression fit via
#' diagnostic plots associated with \link[lme4]{lmer}
#'
#' \link[ComBatFamily]{predict.comfam} for applying ComBat parameters for
#' harmonization of new observations
#'
#' @examples
#' long_combat(iris[,1:2], iris$Species)
#' long_combat(iris[,1:2], iris$Species, iris[3:4], y ~ Petal.Length + Petal.Width)
#'
#' @references
#' Beer, J. C., Tustison, N. J., Cook, P. A., Davatzikos, C., Sheline, Y. I., Shinohara, R. T., & Linn, K. A. (2020). Longitudinal ComBat: A method for harmonizing longitudinal multi-scanner imaging data. *NeuroImage*, 220, 117129. https://doi.org/10.1016/j.neuroimage.2020.117129
long_combat <- function(data, bat, covar = NULL, formula = NULL,
                        eb = TRUE, robust.LS = FALSE, ref.batch = NULL, ...) {
  comfam(data, bat, covar, lmer, formula, eb, robust.LS, ref.batch, ...)
}

#' ComBatLS: Location- and scale-preserving harmonization
#'
#' Implementation of ComBatLS (Gardner et al.) using \link[ComBatFamily]{comfam} with `model` as \link[gamlss]{gamlss}
#' and `family = NO()`. Use `sigma.formula` to specify the scale-preserving model.
#'
#' @param data \emph{n x p} data frame or matrix of observations where
#'   \emph{p} is the number of features and \emph{n} is the number of subjects.
#' @param bat Factor indicating batch (often equivalent to site or scanner)
#' @param covar Data frame or matrix of covariates supplied to \link[gamlss]{gamlss}
#' @param formula Formula for \link[gamlss]{gamlss} starting with `y ~` where `y` represents
#'   each feature
#' @param sigma.formula Formula for variance modeling, formatted following \link[gamlss]{gamlss}
#' @param eb If \code{TRUE}, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization
#' @param robust.LS If \code{TRUE}, uses robust location and scale estimators
#'   for error variance and site effect parameters. Currently uses median and
#'   biweight midvariance
#' @param ref.batch Reference batch, must take value in `levels(bat)`
#' @param ... Additional arguments passed to \link[gamlss]{gamlss}
#'
#' @return `combat` returns a list containing the following components:
#' \item{dat.combat}{Harmonized data as a matrix with same dimensions as `data`}
#' \item{batch.info}{Batch information, including reference batch if specified}
#' \item{fits}{List of model fits from regression step, outputs of \link[gamlss]{gamlss} for each feature}
#' \item{estimates}{List of estimates from standardization and batch effect correction}
#'
#' @import stats gamlss
#' @importFrom methods hasArg
#' @export
#'
#' @seealso
#' \link[ComBatFamily]{plot.comfam} for assessing regression fit via
#' diagnostic plots associated with \link[gamlss]{gamlss}
#'
#' \link[ComBatFamily]{predict.comfam} for applying ComBat parameters for
#' harmonization of new observations
#'
#' @examples
#' combatls(iris[,1:2], iris$Species)
#' combatls(iris[,1:2], iris$Species, iris[3:4], y ~ Petal.Length + Petal.Width,
#'   ~ Petal.Length)
#'
#' @references
#' Gardner et al. to be posted on biorxiv
combatls <- function(data, bat, covar = NULL, formula = NULL,
                     sigma.formula = ~ 1, eb = TRUE, robust.LS = FALSE,
                     ref.batch = NULL, ...) {
  comfam(data, bat, covar, gamlss, formula, eb, robust.LS, ref.batch,
         sigma.formula = sigma.formula, family = NO(), ...)
}
