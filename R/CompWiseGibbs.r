#=======================================================================================================================================#
# Bayesian Variable Selection Approach with Gibbs Sampler                                                                               #
# Date: 02/28/2013                                                                                                                      #
# Maintainer : Kuo-Jung Lee & Ray-Bing Chen                                                                                             #                                                            #
# Description: Performa Bayesian variable selection with component-wise prior to identify the important variable                        #
#=======================================================================================================================================#

CompWiseGibbs = function(Y, X, beta.value, r, tau2, rho, sigma2, Iter)
{
  for(iter in 1:Iter){
    for(var.index in 1:ncol(X)){
      Res = Y- cbind(X[, -var.index]) %*% cbind(beta.value[-var.index, ])
      ri = t(Res)%*%cbind( X[, var.index] )* tau2[var.index] /( sum(X[, var.index]^2)*tau2[var.index] + sigma2 )
      sigma2i.star = sigma2 * tau2[var.index]/( sum(X[, var.index]^2)*tau2[var.index] + sigma2 )
      log.zi = 0.5* ( log(sigma2i.star)-log(tau2[var.index]) ) + ri^2/(2*sigma2i.star)
      prob.0 = rho[var.index]/( rho[var.index]+ (1-rho[var.index])* exp(log.zi) ) 
      r[var.index] = ifelse( runif(1)< prob.0, 0, 1)
      beta.value[var.index] = ifelse(r[var.index], rnorm(1, ri, sqrt(sigma2i.star)), 0)
    }
  }
    list(beta.sample = beta.value, r.sample = r)
}
CompWiseGibbs = cmpfun(CompWiseGibbs)
 
 
 

CompWiseGibbsSimple = function(Y, X, beta.value, r, tau2, rho, sigma2, nu, lambda, num.of.inner.iter, num.of.iteration)
{
	num.of.obs = length(Y)
	num.of.covariates = ncol(X)
	beta.value = c(beta.value)
	beta.samples = matrix(0, num.of.covariates, num.of.iteration)
	beta.samples[, 1] = beta.value

	sigma2.samples = rep(0, num.of.iteration)

	sigma2.samples[1]  = sigma2

	r.samples = matrix(0, num.of.covariates, num.of.iteration)
	r.samples[, 1] = r

  iter = 1
  while( iter< num.of.iteration){
		for(num.of.inner.iter in 1:num.of.inner.iter){
			for(var.index in 1:ncol(X)){
				Res = Y- X[, -var.index, drop=FALSE] %*% cbind(beta.value[-var.index])
				ri = ( t(Res)%*%X[, var.index, drop=FALSE] )* tau2 /( sum(X[, var.index]^2)*tau2 + sigma2 )
				sigma2i.star = sigma2 * tau2/( sum(X[, var.index]^2)*tau2+ sigma2 )
				log.zi = 0.5* ( log(sigma2i.star)-log(tau2) ) + ri^2/(2*sigma2i.star)
				prob.0 = rho/( rho+ (1-rho)* exp(log.zi) ) #(1-rho[var.index])*exp(log.zi) / (rho[var.index] + (1-rho[var.index])*exp(log.zi) )#1/(1 + 1 + (rho[var.index]/(1-rho[var.index]))* exp(log.zi)) #
				r[var.index] = ifelse( (runif(1)< prob.0), 0, 1)
				beta.value[var.index] = ifelse(r[var.index], rnorm(1, ri, sqrt(sigma2i.star)), 0)
			}
		}
    Residual2 = sum((Y-X%*%beta.value)^2)
    sigma2.samples[iter+1] = 1/rgamma(1, (length(Y)+nu)/2, (Residual2+nu*lambda)/2) 
		
		beta.samples[, iter+1] = beta.value
		r.samples[, iter+1] = r

		iter = iter+1
		
  }
  list(beta.sim = beta.samples, sigma2.sim = sigma2.samples, r.sim = r.samples)
}
CompWiseGibbsSimple = cmpfun(CompWiseGibbsSimple)
