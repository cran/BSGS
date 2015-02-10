\name{CompWiseGibbsSimple}
\alias{CompWiseGibbsSimple}
%\keyword{arith}
\title{Running MCMC procedure to generate the posterior samples from the posterior distribution.}
\description{Perform Bayesian variable selection with component-wise prior to identify the important variables.}
\usage{CompWiseGibbsSimple(Y, X, beta.value, r, tau2, rho, sigma2, nu, lambda,
                           num.of.inner.iter, num.of.iteration)}
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
\value{A list is returned with posterior samples of regression coefficients, \eqn{\beta}, variance \eqn{\sigma^2}, and binary variables, \eqn{\gamma}.}
\examples{


\dontrun{CompWiseGibbsSimple(Y, X, beta.value, r, tau2, rho, sigma2, nu, lambda,
 num.of.inner.iter.default, num.of.iteration)}

}
