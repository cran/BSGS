\name{TCR.TPR.FPR.BSGS}
\alias{TCR.TPR.FPR.BSGS}
%\keyword{arith}
\title{Evaluate TCR, TPR and FPR for sparse group variable selection problmes.}
\description{Calculate the true classification rate (TCR), the true positive rate (TPR), and the false positive rate (FPR).}
\usage{TCR.TPR.FPR.BSGS(Output, True.r, Critical.Point)}
\arguments{
\item{Output}{A list of random samples for parameters.}
\item{True.r}{The true value of indicator variable.}
\item{Critical.Point}{When the posterior estiamte of \eqn{r=1} greater than this critical point, then it would be assign to 1, and otherwise 0.}
}
\value{A list is returned with TCP, TPR and FPR.}
\examples{
\dontrun{
output = BSGS.Simple.SaveAllSimulatedSamples(Y, X, Group.Index, r.value, eta.value,
         beta.value, tau2.value, rho.value, theta.value, sigma2.value, nu, lambda, 
         Num.of.Iter.Inside.CompWise, Num.Of.Iteration, MCSE.Sigma2.Given)
TCR.TPR.FPR.BSGS(output, r.true, critical.value)
}

}
