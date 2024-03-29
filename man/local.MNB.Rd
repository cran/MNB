% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Local.R
\name{local.MNB}
\alias{local.MNB}
\title{Local influence}
\usage{
local.MNB(star, formula, dataSet, schemes, cova, plot = TRUE)
}
\arguments{
\item{star}{Initial values for the parameters to be optimized over.}

\item{formula}{The structure matrix of covariates of dimension n x p (in models that include an intercept x
should contain a column of ones).}

\item{dataSet}{data}

\item{schemes}{Perturbation scheme. Possible values: "cases" for Case weights perturbation on ith subject
or cluster, "cases.obs" for Case weights perturbation on jth measurement
taken on the ith subject or cluster, "cova.pertu" for explanatory variable perturbation, "dispersion" for
dispersion parameter perturbation}

\item{cova}{Indicator which column from dataset (continuous covariate) must be perturbation.}

\item{plot}{TRUE or FALSE. Indicates if a graph should be plotted.}
}
\value{
L and graphics
}
\description{
It performes influence analysis by a local influence approach by Cook (1986). It is considering
three perturbation schemes: Case weights, explanatory variable and dispersion parameter perturbation. Another
procedure which considering is the total local curvature corresponding to the ith element approach by
Lesaffre and Verbeke (1998).
}
\details{
The function returns a list (L) with the eigenvector associated with the maximum curvature, the total
local influence and the index plot.
}
\examples{

\donttest{

data(seizures)
head(seizures)

star <-list(phi=1, beta0=1, beta1=1, beta2=1, beta3=1)

local.MNB(formula=Y ~ trt + period + trt:period + offset(log(weeks)),star=star,dataSet=seizures,
schemes="weight",plot=FALSE)

local.MNB(formula=Y ~ trt + period + trt:period + offset(log(weeks)),star=star,dataSet=seizures,
schemes="weight.obs",plot=FALSE)

local.MNB(formula=Y ~ trt + period + trt:period + offset(log(weeks)),star=star,dataSet=seizures,
schemes="dispersion",plot=FALSE)

}
}
\references{
\itemize{
\item Cook, R. D. (1986). Assessment of local influence (with discussion). Journal of
the Royal Statistical Society B, 48, 133-169.
\item Lesaffre E. and Verbeke G. (1998). Local influence in linear mixed models. Biometrics, 54, 570-582.
\item Fabio, L. C., Villegas, C., Carrasco, J. M. F., and de Castro, M. (2021). D
Diagnostic tools for a multivariate negative binomial model for fitting correlated data with
overdispersion. Communications in Statistics - Theory and Methods.
https://doi.org/10.1080/03610926.2021.1939380.
}
}
\author{
Jalmar M F Carrasco <carrascojalmar@gmail.com>,
Cristian M Villegas Lobos <master.villegas@gmail.com> and Lizandra C Fabio <lizandrafabio@gmail.com>
}
