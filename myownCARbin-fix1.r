
## prep the state vars for MCMC

st.names = levels(as.factor(cnew3$state))
nst = length(st.names)
mem = NULL 
for (i in 1:nst){
  mem[[i]] = which(cnew3$state == st.names[i])	  	  
}

## this is a Gibbs sampler for Leroux model coded on my own 
 #   ...CAR1a.r is specifically for deaths data
 #  ...CARbin is for Binomial model
 
  require(MASS)
  require(beepr)
 
   # list of neighbors 
   
   # Cmat = as(co.mat2, "CsparseMatrix")     use Cmat3a which can be read from a file
   K = dim(Cmat3a)[1]
   Nlist = NULL 
   nneighb = numeric(K)
   for (k in 1:K){
      Nlist[[k]] = which(Cmat3a[k,] == 1)	
  	  nneighb[k] = length(Nlist[[k]])
	}  

	
	logfcp.rho = function(r){
	  add1 = sum(log(r*nneighb + (1-r)))
	  add2 = 0 
	  for (k in 1:K){
	    d1 = r*nneighb[k] + 1 - r
		mu1 = r*sumphiki[k]/d1
	    add2 = add2 + d1*(phi[k] - mu1)^2
	  }	  
	  add2 = 1/(tau^2)*add2
	  return((add1 - add2)/2)
	}

	
	logfcp.phik = function(phik,k, mu1, d1){   # can't do vector form?
	   log1e = log(1 + exp(Yhat[k] + phik + Ste[k]))
	   add1 = Y[k]*phik - nk[k]*log1e 
	   add2 = - d1/2/tausq*(phik - mu1)^2
	   add3 = dnorm(phik,0,2, log=TRUE)
	  return(add1 + add2 + add3)
	}
	
	logfcp.st = function(stp,j){   # can't do vector form?
	   memj = mem[[j]]
	   log1e = log(1 + exp(Yhat[memj] + phi[memj] + stp))   # this one's a vector
	   add1 = sum(Y[memj]*stp) - sum(nk[memj]*log1e) 
	   add2 = dnorm(stp,0,nu, log=TRUE)            # should I restrict total st closer to 0? or just to 0, period?
	     meanst = (sum(st) + stp - st[j])/nst
	   add3 = dnorm(meanst,0,0.01, log=TRUE)  
	  return(add1 + add2 + add3)
	}
	
	
	logfcp.beta = function(bet){   # vector form
	   Yhat1 = Xmat %*% bet
	   log1e = log(1 + exp(Yhat1 + phi + Ste))    # note, here Y are the counts and Yhat is in the log-odds space
	   add1 = sum(Y*Yhat1) 
	   add2 = -sum(nk*log1e)
	  return(add1 + add2)
	}

 nk = cnew3$TotalPop

  lm1 = lm(Y ~ X)    # + Popul.mil)   # "state" effect will be dealt with separately; using lm here just to get the Xmat
 
  Xmat = model.matrix(lm1)      #cbind(rep(1,K), X)
  XtXi = solve(t(Xmat) %*% Xmat)
  XtXiXt = XtXi %*% t(Xmat)
  bet = c(-6,-1)    # rep(0, dim(Xmat)[2])
        
  tau = 0.9; tausq = tau^2 	
  nu = 0.5
  rho = 0.5 
  nu0 = 0.1;   df.nu = 2
  tau0 = 1; df.tau = 2   # priors on variances
  rho.step = 0.05
  phi.step = 0.5   # same step for all k, for now
  st.step  = 0.1;    # = 0.05 for cases, 0.2 for deaths?
    beta.step = 0.2   #   = 0.05 for cases, 0.2 for deaths?
  phi = rep(0,K)
  st = rep(0,nst); Ste = numeric(K)

  M = 1000      # M*Nskip = 1000 takes about 45 sec;    a quality run would be M = 5000, Nskip = 5?
  Nskip = 100
  burnin = M*0.2 
 
  phihist = matrix(0,M,K); predhist = matrix(0,M,K); sthist = matrix(0,M,nst)
  nuhist = numeric(M);  tauhist = numeric(M)
  rhohist = numeric(M); betahist = matrix(0,M,dim(Xmat)[2])
  acc.ctr.rho = 0;  acc.ctr.phi = 0; acc.ctr.beta = 0; acc.ctr.st = 0
  
  ptm = proc.time()
  
  # Y = Ytest    # <<<<<   uncomment this line to run in the testing mode

  for(mc in 1:M){              # main MCMC loop
   for (iskip in 1:Nskip){
      
  # sampling beta 
   dbet = beta.step*mvrnorm(1,0*bet,XtXi)   
   betap = bet + dbet
   
   acc.prob = exp(min(logfcp.beta(betap) - logfcp.beta(bet),0))
		if (runif(1) < acc.prob){
		  bet = betap
		  acc.ctr.beta = acc.ctr.beta + 1	
		}
   # bet = bet.true 
   Yhat = Xmat %*% bet    
	  
    # sampling state effects   # logfcp.st = function(stp,j){ 
	
	 for (j in 1:nst){
	   stp = st[j] + runif(1,-1,1)*st.step
	   acc.prob = exp(min(logfcp.st(stp,j) - logfcp.st(st[j],j),0))
		if (runif(1) < acc.prob){
		  st[j] = stp
		  acc.ctr.st = acc.ctr.st + 1	
		}
		Ste[ mem[[j]] ] = st[j] 
	 }
	# sampling nu
      ssq = sum(st^2) 
      nusq = (ssq + nu0^2*df.nu)/rchisq(1,df.nu + nst)   # inverse chisq
	  nu = sqrt(nusq)	
	 # nu = 0.01
	 
    # sampling phik 
# if (0 == 1){	
	   for (k in 1:K){
		 neighb = Nlist[[k]]
		 d1 = rho*length(neighb) + 1 - rho 
		 mu1 = rho*sum(phi[neighb])/d1

		 phiprop = phi[k] + runif(1,-1,1)*phi.step 	 
		 acc.prob = exp(min(logfcp.phik(phiprop,k, mu1, d1) - logfcp.phik(phi[k],k, mu1, d1),0))
			if (runif(1) < acc.prob){
			  phi[k] = phiprop
			  acc.ctr.phi = acc.ctr.phi + 1	
			}
	   }
 # } 
   #  phi = phitrue
   # sampling tau 
   
   SStau = 0
   sumphiki = numeric(K)
   for (k in 1:K){
     neighb = Nlist[[k]]   
	 denom = rho*length(neighb) + 1 - rho    ## or nneighb[k]
	 sumphiki[k] = sum(phi[neighb])
	 mu1 = rho*sumphiki[k]/denom 
	 SStau = SStau + (phi[k] - mu1)^2*denom       
   }
   tausq = (SStau + tau0^2*df.tau)/rchisq(1, df.tau + K)
   tau = sqrt(tausq)
   # tau = tautrue
   
   # sampling rho 
   
    ## stepr = min(rho.step, rho/2, (1-rho)/2)  # "adaptive" version ??
	stepr = rho.step
    rhop = rho + runif(1,-stepr, stepr)
	
	if ( (rhop < 0.9999) & (rhop > 0)){
		acc.prob = exp(min(logfcp.rho(rhop) - logfcp.rho(rho),0))
		if (runif(1) < acc.prob){
		  rho = rhop
		  acc.ctr.rho = acc.ctr.rho + 1	
		}
    }  
    # rho = rhotrue
   
   # posterior predictive 
   
    logitpk = Xmat %*% bet + phi + Ste
	pk = 1 - 1/(1+exp(logitpk))
	ypred = rbinom(K,nk,pk)   
   }   # end skip
   predhist[mc,] = ypred
   phihist[mc,] = phi
   nuhist[mc] = nu
   tauhist[mc] = tau
   rhohist[mc] = rho
   betahist[mc,] = bet
   sthist[mc,] = st
  }   ## of the main MCMC loop
   
   proc.time() - ptm 
plot(rhohist[50:M])
acc.ctr.rho
 beepr::beep("mario")

 print(paste("Acc.pct.beta=",round(acc.ctr.beta/M/Nskip,4)," Acc.pct.st=",round(acc.ctr.st/nst/M/Nskip,4)," Acc.pct.rho=",round(acc.ctr.rho/M/Nskip,4)))
plot(betahist[50:M,1],betahist[50:M,2], col=rgb((50:M)/M,0, 1 - (50:M)/M))
     
  ##  CI output 
   rnge = burnin:M 
   CI.rho = quantile(rhohist[rnge], probs = c(0.025,0.975))
   CI.tau = quantile(tauhist[rnge], probs = c(0.025,0.975))
   CI.nu = quantile(nuhist[rnge], probs = c(0.025,0.975))
   CI.bet = NULL
   for (i in 1:dim(Xmat)[2]){
      CI.bet = rbind(CI.bet, quantile(betahist[rnge,i], probs = c(0.025,0.975)))   
   }
   row.names(CI.bet) = colnames(Xmat)
   
   round(rbind(CI.rho, CI.tau, CI.nu, CI.bet),4)
    
   mean(betahist[rnge,1])
   mean(betahist[rnge,2])
   mean(rhohist[rnge])
   mean(tauhist[rnge])
   mean(nuhist[rnge])
   
   
   yhat = apply(predhist,2, "median")
   plot(yhat,Y,log="xy")
   
   # plot(sthist[,2]) 
   
     
  ## run to here with the main dataset # =========================================
   # ======== #
  
  
  
  
  
   
   
   # -------------------------------------
   ## creating the test data 
   
	 bet.true = c(-6.6, -1.87) 
	 
	 phi = rep(0,K)
	 rho = 0.5;  tau = 0.6

  M = 500 
  ptm = proc.time() 

  for(mc in 1:M){              # main MCMC loop
  
    # sampling phik 
   for (k in 1:K){
     neighb = Nlist[[k]]
	 denom = rho*length(neighb) + 1 - rho 
	 mu1 = rho*sum(phi[neighb])/denom 
     phi[k] = rnorm(1,mu1, tau/sqrt(denom))           
   }
  }   ## of the main MCMC loop
  phitrue = phi 
  rhotrue = rho;   tautrue = tau
   proc.time() - ptm 
	
	pk = 1 - 1/(1 + exp(Xmat %*% bet.true + phi))
	Ytest = 0*Y
	
   for (k in 1:K){
     Ytest[k] = rbinom(1,nk[k], pk[k])
    }   
	
   
    ## running test data through CARBAyes: seems to work ok!
	
	 formula <-  Ytest ~ vaxx0 
     model <- S.CARleroux(formula=formula, family="binomial", trials= cnew3$TotalPop, W = co.mat3a, burnin=2000, n.sample=10000, thin = 100)
 
 
     
     # Assuming your MCMC output is stored in a variable named 'mcmc_output'
     
     # Extract the number of iterations and parameters
     num_iterations <- nrow(mcmc_output)
     num_parameters <- ncol(mcmc_output)
     
     # Create the first chain
     chain_1 <- mcmc(mcmc_output)
     
     # Simulate a new chain starting from the last state of the original chain
     start_value_chain_2 <- as.vector(tail(chain_1, 1))
     new_chain_2_values <- matrix(rnorm(num_iterations * num_parameters, mean = start_value_chain_2, sd = 0.01),
                                  ncol = num_parameters)
     new_chain_2 <- mcmc(new_chain_2_values)
     
     # Combine the two chains into a list
     mcmc_combined <- mcmc.list(chain_1, new_chain_2)
     
     # Check the resulting MCMC object
     print(mcmc_combined)
     
     # Check the Gelman-Rubin diagnostic
     gelman_diag_result <- gelman.diag(mcmc_combined)
     print(gelman_diag_result)
     