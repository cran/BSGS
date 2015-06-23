\name{BSGS.Sample}
\alias{BSGS.Sample}
%\keyword{arith}
\title{Sample version of group-wise Gibbs sampler for sparse group selection.}
\description{Generate the posterior samples by an approximation sampling method to perform Bayesian sparse group selection to identify the important groups of variables and variables within those.}
\usage{BSGS.Sample(Y, X, Group.Index, r.value, eta.value, beta.value, tau2.value, 
rho.value, theta.value, sigma2.value, nu, lambda, Num.of.Iter.Inside.CompWise, 
Num.Of.Iteration, MCSE.Sigma2.Given)
}
\arguments{
\item{Y}{vector of observations of length \eqn{n}.}
\item{X}{design matrix of dimension \eqn{n \times p}.}
\item{Group.Index}{Specify the group label to which the variable belongs.}
\item{r.value}{Initial values of indicator variables for individual variables.}
\item{eta.value}{Initial values of indicator variables for the group variables.}
\item{beta.value}{Initial values of regression coefficients, \eqn{\beta}.}
\item{tau2.value}{Variance in the prior distribution for regression coefficients.}
\item{rho.value}{Prior inclusion probability for a variable.}
\item{theta.value}{Prior inclusion probability for a group.}
\item{sigma2.value}{Initial value of \eqn{\sigma^2}.}
\item{nu}{The hyperparameter in the prior distribution of \eqn{\sigma^2}.}
\item{lambda}{The hyperparameter in the prior distribution of \eqn{\sigma^2}.}
\item{Num.of.Iter.Inside.CompWise}{Specify the number of iterations within component wise Gibbs sampler for variable selection within a group.}

\item{Num.Of.Iteration}{Specify the number of iterations for sparse group variable selection.}
\item{MCSE.Sigma2.Given}{Prespecified value which is used to stop simulating samples. When the MCSE of estimate of \eqn{\sigma^2} less than the given value, the simulation is terminated.}
}
\value{A list is returned with posterior random samples of regression coefficients, \eqn{\beta}, binary variables for group selection, \eqn{\eta}, binary variables for variable selection, \eqn{\gamma}, variance, \eqn{\sigma^2} and the number of iterations performed and the elapsed time in second required for the run.}

\examples{


\dontrun{
output = BSGS.Sample(Y, X, Group.Index, r.value, eta.value, beta.value, tau2.value, 
	rho.value, theta.value, sigma2.value, nu, lambda, 
	Num.of.Iter.Inside.CompWise, Num.Of.Iteration, MCSE.Sigma2.Given)
}

}
