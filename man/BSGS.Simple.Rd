\name{BSGS.Simple}
\alias{BSGS.Simple}
%\keyword{arith}
\title{The group-wise Gibbs sampler for sparse group selection .}
\description{Perform Bayesian sparse group selection to identify the important groups of variables and variables within those.}
\usage{BSGS.Simple(Y, X, Group.Index, r.value, eta.value, beta.value, tau2.value,
rho.value, theta.value, sigma2.value, nu, lambda, Num.of.Iter.Inside.CompWise,
Num.Of.Iteration, Burn.In)

}
\arguments{
\item{Y}{Vector of observations of length \eqn{n}.}
\item{X}{Design matrix of dimension \eqn{n \times p}.}
\item{Group.Index}{Specify the group that each variable within design matrix belongs to.}
\item{r.value}{Initial values of indicator variables for individual variables.}
\item{eta.value}{Initial values of indicator variable for group variables.}
\item{beta.value}{Initial values of regression coefficients, \eqn{\beta}.}
\item{tau2.value}{Variance in the prior distribution for regression coefficients.}
\item{rho.value}{Prior probability including a variable.}
\item{theta.value}{Prior probability including a group.}
\item{sigma2.value}{Variance in the regression model.}
\item{nu}{Given value in the prior distribution of \eqn{\sigma^2}.}
\item{lambda}{Given value in the prior distribution of \eqn{\sigma^2}.}
\item{Num.of.Iter.Inside.CompWise}{Specify the number of iterations to be runned within component wise Gibbs sampler for variable selection within a group.}
\item{Num.Of.Iteration}{Specify the number of iterations to be runned for sparse group variable selection.}
\item{Burn.In}{Specify the number of iterations to be discarded.}
}
\value{A list is returned with estimates of regression coefficients, \eqn{\beta},  binary variables for group selection, \eqn{\eta}, binary variables for variable selection, \eqn{r},
and variance, \eqn{\sigma^2}.}
\examples{


\dontrun{

output = BSGS.Simple(Y, X, Group.Index, r.value, eta.value, beta.value, tau2.value, rho.value,
theta.value, sigma2.value, nu, lambda, Num.of.Iter.Inside.CompWise, Num.Of.Iteration, Burn.In)
output$beta.est
output$eta.est
output$r.est
output$sigma2.est

}

}
