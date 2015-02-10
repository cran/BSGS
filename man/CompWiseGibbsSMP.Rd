\name{CompWiseGibbsSMP}
\alias{CompWiseGibbsSMP}
%\keyword{arith}
\title{Stochastic matching pursuit for variable selection.}
\description{Perform MCMC procedure to generate the posterior samples to estimate posterior quantities
of interest in Bayesian variable selection using stochastic matching pursuit approach.}
\usage{CompWiseGibbsSMP(Y, X, beta.value, r, tau2, rho, sigma2, nu, lambda, 
num.of.inner.iter, num.of.iteration)
}
\arguments{
\item{Y}{vector of observations of length \eqn{n}.}
\item{X}{design matrix of dimension \eqn{n \times p}.}
\item{beta.value}{Initial values of regression coefficients, \eqn{\beta}.}
\item{r}{Initial values of indicator variables for individual regressors.}
\item{tau2}{Variance in the prior distribution for regression coefficients.}
\item{rho}{Prior probability including a variable.}
\item{sigma2}{Initial value of \eqn{\sigma^2}.}
\item{nu}{Given value in the prior distribution of \eqn{\sigma^2}.}
\item{lambda}{Given value in the prior distribution of \eqn{\sigma^2}.}
\item{num.of.inner.iter}{The number of iterations before sampling \eqn{\sigma^2}.}
\item{num.of.iteration}{The number of iterations to be runned for sparse group variable selection.}


}
\value{A list is returned with estimates of regression coefficients, \eqn{\beta}, binary variables for variable selection, \eqn{\gamma},
and variance, \eqn{\sigma^2}.}
\examples{


\dontrun{CompWiseGibbsSMP(Y, X, beta.value, r, tau2, rho, sigma2, nu0, lambda0, 
num.of.inner.iter, num.of.iteration)}

}
