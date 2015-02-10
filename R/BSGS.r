#=======================================================================================================================================#
# Bayesian Sparse Group Selection                                                                                                       #
# Date: 01/15/2014                                                                                                                      #
# Maintainer : Kuo-Jung Lee & Ray-Bing Chen                                                                                             #                                                            #
# Description: Perform Bayesian sparse group selection                                                                                  #
#=======================================================================================================================================#
BSGS = function(Y, X, Group.Index, r.value, eta.value, beta.value, tau2.value, rho.value, theta.value, sigma2.value, nu, lambda, Num.of.Iter.Inside.CompWise, Num.Of.Iteration, Burn.In)
{               
  if(Num.Of.Iteration<Burn.In) stop("Numer of iterations should be greater than the number of burn-in.") 
  num.of.obs = length(Y)
  num.of.covariates = ncol(X)
  group.inf = (count(Group.Index))
  num.of.groups = dim(group.inf)[1]
  group.lables = group.inf[, 1]
  num.of.var.in.groups = group.inf[, 2]
  

  r.sample = r.value# r.true = as.numeric(beta.value!=0)
  beta.sample = beta.value
  sigma2.sample = sigma2.value
  eta.sample = eta.value #rep(0, num.of.group.var) 

  
  
  colnames(X) = names(r.sample) = names(beta.sample)= names(tau2.value) = names(rho.value) = Group.Index
  names(eta.sample) = names(theta.value) =  group.lables
  beta.sample = cbind(beta.sample)

	if(Burn.In){
		r.est = r.value/Num.Of.Iteration #rep(0, num.of.covariates)/num.of.iter
		eta.est = rep(0, num.of.groups)/Num.Of.Iteration
		beta.est = beta.value/Num.Of.Iteration
		sigma2.est = sigma2.value/Num.Of.Iteration
	}
	else{
		r.est = 0
		eta.est = 0
		beta.est = 0
		sigma2.est = 0
	}
  iter = 1
  while(iter < Num.Of.Iteration){
  
  ## Step I
    group.selection = sample(group.lables, 1)

    Res = Y- (X[, Group.Index != group.selection]) %*% cbind(beta.sample[Group.Index != group.selection])

  ## Step II & III
    
    X.cws = cbind(X[, Group.Index == group.selection])
    beta.sample.cws = cbind(beta.sample[Group.Index == group.selection, ])
    r.sample.cws = r.sample[Group.Index == group.selection]
    rho.value.cws = rho.value[Group.Index == group.selection]
    theta.value.cws = theta.value[names(eta.sample) == group.selection]
    tau2.value.cws = tau2.value[Group.Index == group.selection]

 
    if(eta.sample[names(eta.sample) == group.selection] == 0){ # Add a group
      beta.r.cand = CompWiseGibbs(Res, X.cws, beta.sample.cws, r.sample.cws, tau2.value.cws, rho.value.cws, sigma2.sample, Num.of.Iter.Inside.CompWise)
      if(all(beta.r.cand$r.sample==0))
        eta.cand = 0
      else
        eta.cand = GroupMH(Res, X.cws, beta.r.cand$beta.sample, beta.r.cand$r.sample, tau2.value.cws, rho.value.cws, sigma2.sample, theta.value.cws, 0)

      if(eta.cand == 1){
        beta.sample[Group.Index == group.selection, ] = beta.r.cand$beta.sample #beta.true[group.index == group.selection]#
        r.sample[Group.Index == group.selection] = beta.r.cand$r.sample#r.true[group.index == group.selection]#
      }
      else{
        beta.sample[Group.Index == group.selection, ] = 0 #beta.true[group.index == group.selection]#
        r.sample[Group.Index == group.selection] = 0 #r.true[group.index == group.selection]
      }

    }
    else{ # Delete a group
      eta.cand = GroupMH(Res, X.cws, beta.sample.cws, r.sample.cws, tau2.value.cws, rho.value.cws, sigma2.sample, theta.value.cws, 1)
      if(eta.cand == 0){
        beta.sample[Group.Index == group.selection, ] = 0 #
        r.sample[Group.Index == group.selection] = 0 #
      }
      else{
        beta.r.cand = CompWiseGibbs(Res, X.cws, beta.sample.cws, r.sample.cws, tau2.value.cws, rho.value.cws, sigma2.sample, Num.of.Iter.Inside.CompWise)
        beta.sample[Group.Index == group.selection, ] = beta.r.cand$beta.sample # beta.true[group.index == group.selection]#
        r.sample[Group.Index == group.selection] =  beta.r.cand$r.sample#r.true[group.index == group.selection]#
      }

    }

   ## Step IV
   
    SSE = sum( (Y- X%*% beta.sample)^2 )
    sigma2.sample = 1/rgamma(1, ((num.of.obs+nu)/2), (SSE+nu*lambda)/2)#sigma2.true #
    
    eta.sample[names(eta.sample) == group.selection] = eta.cand

    if(iter > Burn.In){
      eta.est = eta.est +  eta.sample/(Num.Of.Iteration-Burn.In)
      beta.est = beta.est + beta.sample/(Num.Of.Iteration-Burn.In)
      r.est = r.est + r.sample/(Num.Of.Iteration-Burn.In)
      sigma2.est = sigma2.est + sigma2.sample/(Num.Of.Iteration-Burn.In)
    }
    iter = iter + 1

  }
  list(beta.est = beta.est, eta.est = eta.est, r.est = r.est, sigma2.est = sigma2.est)
}


