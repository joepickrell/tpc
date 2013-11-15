
library(MCMCpack)
library(mvtnorm)

run_mcmc = function(ca, cb, init_theta, dprior, cprior, cpropose, tpropose, fpropose, N, Nburn, Nsamp, outfile){
	toreturn = vector()
	toreturns = vector()	
	theta = init_theta
	theta_bk = theta
	lik = llk(ca, cb, theta)
	lik_bk = lik
	total = 0
	acceptfa =0
	acceptc = 0
	accepttb = 0
	for (i in 1:N){
		total = total+1
		if (i %% Nsamp ==0 & i > Nburn){
			lik = llk(ca, cb, theta)
			print(paste(i, lik, theta[2]))
			print(theta[1][[1]])
			print(theta[3][[1]])
        		
			s = sampab(theta, sum(ca), sum(cb))
                        sa = s[1][[1]]
                        sb = s[2][[1]]
			sp = sa+sb
			tmpsa = sa[sp>0]
			tmpsb = sb[sp>0]
			c = chisq.test(rbind(tmpsa, tmpsb))
			c = c$statistic
			toreturn = append(toreturn, c)
			toreturns = append(toreturns, sum(sa >0))
			sink(file = outfile, append = T)
			id = paste("iter", i, sep = "")
			cat(id, "llik:", lik, "\n")
			cat(id, "c:", exp(theta[2][[1]]), "\n")
			cat(id, "fa:", theta[1][[1]], "\n")
			cat(id, "tb:", theta[3][[1]], "\n")	
			cat(id, "sa:", sa, "\n")
			cat(id, "sb:", sb, "\n")
			cat(id, "schisq:", c, "\n")		
			sink()	
		}
		theta= theta_bk

		# propose a new fa
		# proposal is a 99%/1% mixture of dirichlet distributions sampled from the current fa and rep(0.5, N) (this avoids hitting 0 and getting stuck)
		# NB: not a symmetric proposal

		theta[1][[1]] = propose_fa(theta[1][[1]], fpropose)
		lik = llk(ca, cb, theta) + log(ddirichlet(theta[1][[1]], dprior)) + ldirmix(theta_bk[1][[1]], theta[1][[1]], fpropose)
		if ( is.nan(lik)){ lik = -Inf}
                oldlik = llk(ca, cb, theta_bk)+ log(ddirichlet(theta_bk[1][[1]], dprior)) + ldirmix(theta[1][[1]], theta_bk[1][[1]], fpropose)

		
		accept = mhupdate(lik, oldlik)
                if (accept){
			acceptfa = acceptfa+1
                        theta_bk = theta
                        lik_bk = lik
                }
                else{
                        theta = theta_bk
                        lik = lik_bk
                }


		#update c
		#new logc is normally dsitributed around the old logc
		#NB: symmetric
		theta[2][[1]] = propose_c(theta[2][[1]], cpropose)
		lik = llk(ca, cb, theta) + cprob(theta[2][[1]], cprior)
		oldlik = llk(ca, cb, theta_bk)+ cprob(theta_bk[2][[1]], cprior)
		accept = mhupdate(lik, oldlik)
		if (accept){
			acceptc = acceptc+1
			theta_bk = theta
			lik_bk = lik
		}
		else{
			theta = theta_bk
			lik = lik_bk
		}

		# update ta
		# tas are normally distributed around old tas
		# symmetric (also no prior, since it enters through fa)
		theta[3][[1]] = propose_tb(theta[3][[1]], tpropose)
		
		# need to make sure this is a legal move (ie. frequences sum to 1)
		ok = tbok(theta[3][[1]])
		if (ok == FALSE){
			theta = theta_bk
			lik = lik_bk
		}
		else{
			lik =  llk(ca, cb, theta)
                	oldlik = llk(ca, cb, theta_bk)
               		 accept = mhupdate(lik, oldlik)
                	if (accept){
				accepttb = accepttb+1
                        	theta_bk = theta
                        	lik_bk = lik
                	}
                	else{
                        	theta = theta_bk
                        	lik = lik_bk
                	}
		}	
	}
	sink(file = outfile, append = T)
	cat("acceptfa:", acceptfa/total, "\n")
	cat("acceptc:", acceptc/total, "\n")
	cat("accepttb:", accepttb/total, "\n")
	sink()	
	return(list(toreturn, toreturns))

}


sampab = function(theta, Na, Nb){
	fa = theta[1][[1]]
	tb = theta[3][[1]]
	fb = tb
	fb[fb <0] = 0
	fb[fb >1] = 1
	ca = rmultinom(1, Na, prob = fa)
	cb = rmultinom(1, Nb, prob = fb)
	return(list(ca, cb))

}
propose_fa = function(fa, fpropose){
	toreturn  = 0.99*rdirichlet(1, fpropose*fa)[1,]+0.01* rdirichlet(1, rep(0.5, length(fa)))[1,]
	return(toreturn)

}

ldirmix = function(y, x, frpopose){
	toreturn = 0.99* ddirichlet(y, x*fpropose) + 0.01* ddirichlet(y, rep(0.5, length(y)))
	return(log(toreturn))

}
tbok = function(tb){
	fb = tb[1:(length(tb)-1)]
	fb[fb < 0] = 0
	fb[fb > 1] = 1
	if (sum(fb) > 1){ return(FALSE)}
	else{ return(TRUE)}

}
mhupdate = function(lnum, ldenom){
	test = exp(lnum - ldenom)
	if (test > 1) { return(TRUE)}
	else{
		ru = runif(1)
		if (ru < test){ return(TRUE)}
		else{ return(FALSE)}
	}

}

cprob = function(logc, cprior){
 	toreturn = dnorm(logc, log(cprior), sd = 1, log = T)
	return(toreturn)
}

simdata = function(theta, Na, Nb){
        fa = theta[1][[1]]
        logc = theta[2][[1]]
        tb = theta[3][[1]]
        fb = tb
        fb[fb<0] = 0
        fb[fb>1] = 1
	sa = rmultinom(1, Na, fa)
	sb = rmultinom(1, Nb, fb)
	return(list(sa, sb))
	

}

llk = function(ca, cb, theta){ # theta is fa, c, tb
	fa = theta[1][[1]]
	logc = theta[2][[1]]
	tb = theta[3][[1]]
	fb = tb
	fb[fb<0] = 0
	fb[fb>1] = 1

        toreturn = dmultinom(ca, prob = fa, log = T)
	# make sure all fas are greater than 0
	tmpw= which(fa>1e-10)
	tmpfa = fa[tmpw]
	tmptb = tb[tmpw]

	sigma = get_sigma(tmpfa, exp(logc))
	#print(tb)
	#print(logc)
	#print(fa)
	#print(sigma)
	toreturn = toreturn+ dmvnorm(tmptb[1:(length(tmptb)-1)], tmpfa[1: (length(tmpfa)-1)], sigma, log = T) 
	#print("nothere")
	toreturn = toreturn  + dmultinom(cb, prob = fb, log = T)
	return(toreturn)
}

get_sigma = function(fa, c){
	sigma = matrix(nrow = length(fa)-1, ncol = length(fa)-1) #naive sigma is singuar because sum(fa) == 1
        for (i in 1: (length(fa)-1)){
                for (j in i:(length(fa)-1)){
                        if (i == j){ sigma[i,j] = c*fa[i]*(1-fa[i])}
                        else{ 
                                sigma[i,j] = -c*fa[i]*fa[j]
                                sigma[j,i] = -c*fa[i]*fa[j]
                        }
                }
        }
        return(sigma)
}
propose_c = function(logc, var){
	newc = rnorm(1, mean = logc, sd = sqrt(var))
	return (newc)
	
	
}

gibbs_fa = function(ca, dprior){
	newfa = rdirichlet(1, ca+dprior)
	return(newfa)
}

propose_tb = function(tb, var){
	fb = tb
	fb[fb <0] = 0
	fb[fb>1] = 1
	sigma = get_sigma(fb, var)
	for (i in 1:(length(fb)-1)){
		sigma[i,] = 0
		sigma[i,i] = var
	}
	newtb = mvrnorm(n = 1, mu = tb[1:(length(tb)-1)], Sigma = sigma)
	tmp = newtb
	tmp[tmp < 0] = 0
	tmp[tmp > 1] = 1
	newtb = append(newtb, 1-sum(tmp))
	return(newtb)
} 


