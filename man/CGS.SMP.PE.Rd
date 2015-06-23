\name{CGS.SMP.PE}
\alias{CGS.SMP.PE}
%\keyword{arith}
\title{Posterior estimates of parameters.}
\description{Calculate the posterior estimates of parameters based on the samples generated from the posterior distribution by the stochastic matching pursuit (SMP).}
\usage{CGS.SMP.PE(CGS.SMP.Output)}
\arguments{
\item{CGS.SMP.Output}{A list of random samples for parameters}
}
\value{A list is returned with estimates of regression coefficients, \eqn{\beta}, binary variables for variable selection, \eqn{\gamma},
and variance, \eqn{\sigma^2}.}
\examples{
\dontrun{
output = CompWiseGibbsSMP(Y, X, beta.value, r, tau2, rho, sigma2, nu, lambda, 
	num.of.inner.iter, num.of.iteration, MCSE.Sigma2.Given)
CGS.SMP.PE(output)
}

}


