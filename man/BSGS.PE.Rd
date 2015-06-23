\name{BSGS.PE}
\alias{BSGS.PE}
%\keyword{arith}
\title{Posterior estimates of parameters.}
\description{Provide the posterior estimates of parameters.}
\usage{BSGS.PE(BSGS.Output)}
\arguments{
\item{BSGS.Output}{A list of random samples generated from the posterior distribution by MCMC procedures.}
}
\value{A list is returned with estimates of regression coefficients, \eqn{\beta}, the posterior probability of binary variable \eqn{\eta} for group selection equal to 1, binary variable \eqn{\gamma} for variable selection equal to 1,
and variance, \eqn{\sigma^2}.}
\examples{
\dontrun{
output = BSGS.Simple(Y, X, Group.Index, r.value, eta.value, beta.value, tau2.value, 
	rho.value, theta.value, sigma2.value, nu, lambda, 
	Num.of.Iter.Inside.CompWise, Num.Of.Iteration, MCSE.Sigma2.Given)
BSGS.PE(output)
}

}


