% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{long_combat}
\alias{long_combat}
\title{Longitudinal ComBat harmonization}
\usage{
long_combat(
  data,
  bat,
  covar = NULL,
  formula = NULL,
  eb = TRUE,
  robust.LS = FALSE,
  ref.batch = NULL,
  ...
)
}
\arguments{
\item{data}{\emph{n x p} data frame or matrix of observations where
\emph{p} is the number of features and \emph{n} is the number of subjects.}

\item{bat}{Factor indicating batch (often equivalent to site or scanner)}

\item{covar}{Data frame or matrix of covariates supplied to \link[lme4]{lmer}}

\item{formula}{Formula for \link[lme4]{lmer} starting with \verb{y ~} where \code{y} represents
each feature}

\item{eb}{If \code{TRUE}, uses ComBat model with empirical Bayes for mean
and variance harmonization}

\item{robust.LS}{If \code{TRUE}, uses robust location and scale estimators
for error variance and site effect parameters. Currently uses median and
biweight midvariance}

\item{ref.batch}{Reference batch, must take value in \code{levels(bat)}}

\item{...}{Additional arguments passed to \link[lme4]{lmer}}
}
\value{
\code{combat} returns a list containing the following components:
\item{dat.combat}{Harmonized data as a matrix with same dimensions as \code{data}}
\item{batch.info}{Batch information, including reference batch if specified}
\item{fits}{List of model fits from regression step, outputs of \link[lme4]{lmer} for each feature}
\item{estimates}{List of estimates from standardization and batch effect correction}
}
\description{
Implementation of Longitudinal ComBat (Beer et al., 2020) using \link[ComBatFamily]{comfam} with \code{model} as \link[lme4]{lmer}.
Currently, this implementation is lacking the REML option for variance estimates.
Please use the longComBat package (https://github.com/jcbeer/longCombat) if this is needed.
}
\examples{
long_combat(iris[,1:2], iris$Species)
long_combat(iris[,1:2], iris$Species, iris[3:4], y ~ Petal.Length + Petal.Width)

}
\references{
Beer, J. C., Tustison, N. J., Cook, P. A., Davatzikos, C., Sheline, Y. I., Shinohara, R. T., & Linn, K. A. (2020). Longitudinal ComBat: A method for harmonizing longitudinal multi-scanner imaging data. \emph{NeuroImage}, 220, 117129. https://doi.org/10.1016/j.neuroimage.2020.117129
}
\seealso{
\link[ComBatFamily]{plot.comfam} for assessing regression fit via
diagnostic plots associated with \link[lme4]{lmer}

\link[ComBatFamily]{predict.comfam} for applying ComBat parameters for
harmonization of new observations
}
