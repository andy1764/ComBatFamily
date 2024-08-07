% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{combatls}
\alias{combatls}
\title{ComBatLS: Location- and scale-preserving harmonization}
\usage{
combatls(
  data,
  bat,
  covar = NULL,
  formula = NULL,
  sigma.formula = ~1,
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

\item{covar}{Data frame or matrix of covariates supplied to \link[gamlss]{gamlss}}

\item{formula}{Formula for \link[gamlss]{gamlss} starting with \verb{y ~} where \code{y} represents
each feature}

\item{sigma.formula}{Formula for variance modeling, formatted following \link[gamlss]{gamlss}}

\item{eb}{If \code{TRUE}, uses ComBat model with empirical Bayes for mean
and variance harmonization}

\item{robust.LS}{If \code{TRUE}, uses robust location and scale estimators
for error variance and site effect parameters. Currently uses median and
biweight midvariance}

\item{ref.batch}{Reference batch, must take value in \code{levels(bat)}}

\item{...}{Additional arguments passed to \link[gamlss]{gamlss}}
}
\value{
\code{combat} returns a list containing the following components:
\item{dat.combat}{Harmonized data as a matrix with same dimensions as \code{data}}
\item{batch.info}{Batch information, including reference batch if specified}
\item{fits}{List of model fits from regression step, outputs of \link[gamlss]{gamlss} for each feature}
\item{estimates}{List of estimates from standardization and batch effect correction}
}
\description{
Implementation of ComBatLS (Gardner et al.) using \link[ComBatFamily]{comfam} with \code{model} as \link[gamlss]{gamlss}
and \code{family = NO()}. Use \code{sigma.formula} to specify the scale-preserving model.
}
\examples{
combatls(iris[,1:2], iris$Species)
combatls(iris[,1:2], iris$Species, iris[3:4], y ~ Petal.Length + Petal.Width,
  ~ Petal.Length)

}
\references{
Gardner, M., Shinohara, R. T., Bethlehem, R. A. I., Romero-Garcia, R., Warrier, V., Dorfschmidt, L., Shanmugan, S., Seidlitz, J., Alexander-Bloch, A., & Chen, A. A. (2024). ComBatLS: A location- and scale-preserving method for multi-site image harmonization. bioRxiv, 2024.06.21.599875. https://doi.org/10.1101/2024.06.21.599875
}
\seealso{
\link[ComBatFamily]{plot.comfam} for assessing regression fit via
diagnostic plots associated with \link[gamlss]{gamlss}

\link[ComBatFamily]{predict.comfam} for applying ComBat parameters for
harmonization of new observations
}
