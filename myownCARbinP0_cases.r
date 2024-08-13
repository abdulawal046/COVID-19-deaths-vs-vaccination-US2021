
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
	   add3 = dnorm(meanst,0,0.002, log=TRUE)  
	  return(add1 + add2 + add3)
	}
	
	
	logfcp.beta = function(bet){   # vector form
	   Yhat1 = Xmat %*% bet
	   log1e = log(1 + exp(Yhat1 + phi + Ste))    # note, here Y are the counts and Yhat is in the log-odds space
	   add1 = sum(Y*Yhat1) 
	   add2 = -sum(nk*log1e)
	  return(add1 + add2)
	}
	
 
 #Y =  round(totDeth * cnew3$TotalPop/1e5)           # round(cnew3$case.pAug * cnew3$TotalPop) 
            # round(totDeth * cnew3$TotalPop/1e5)     #

 # Y= round(cnew3$dethSep * cnew3$TotalPop/1e5)   
 vaxx = cnew3$vacJul/100    # rescaled to between 0 and 1
 meanvaxx = mean(vaxx)
 X = vaxx - meanvaxx
 Popul.mil = cnew3$TotalPop/1e+6
 nk = cnew3$TotalPop

  lm1 = lm(Y ~ X)    # + Popul.mil)   # "state" effect will be dealt with separately; using lm here just to get the Xmat
 
  Xmat = model.matrix(lm1)      #cbind(rep(1,K), X)
  XtXi = solve(t(Xmat) %*% Xmat)
  XtXiXt = XtXi %*% t(Xmat)
  bet = c(b0s,-1)    # rep(0, dim(Xmat)[2])
        
  tau = 0.9; tausq = tau^2 	
  nu = 0.5
  rho = 0.5 
  nu0 = 0.1;   df.nu = 2
  tau0 = 1; df.tau = 2   # priors on variances
  rho.step = 0.1
  phi.step = 0.5   # same step for all k, for now
  st.step  = 0.02    # = 0.02 for cases, 0.2 for deaths?
    beta.step = 0.025   #   = 0.02 for cases, 0.2 for deaths?
  phi = rep(0,K)
  st = rep(0,nst); Ste = numeric(K)

  M = 1000      # M*Nskip = 1000 takes about 47 sec;    a quality run would be M = 5000, Nskip = 20?
  Nskip = 300  # 100
  burnin = M*0.2;   M10 = M %/% 10
 
  pred0hist = matrix(0,M,K); pred10hist = matrix(0,M,K); pred20hist = matrix(0,M,K);
  pred30hist = matrix(0,M,K); pred40hist = matrix(0,M,K); pred50hist = matrix(0,M,K); 
  pred60hist = matrix(0,M,K); pred70hist = matrix(0,M,K); pred80hist = matrix(0,M,K);
  pred90hist = matrix(0,M,K); pred100hist = matrix(0,M,K);
  
  phihist = matrix(0,M,K); predhist = matrix(0,M,K); sthist = matrix(0,M,nst)
  nuhist = numeric(M);  tauhist = numeric(M)
  rhohist = numeric(M); betahist = matrix(0,M,dim(Xmat)[2])
  acc.ctr.rho = 0;  acc.ctr.phi = 0; acc.ctr.beta = 0; acc.ctr.st = 0
  
  ptm = proc.time()
  
  # Y = Ytest    # to run in the testing mode

  for(mc in 1:M){              # main MCMC loop
   if ((mc %% M10) == 0)  print (paste(mc/M10, "finished"))
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

    logitpk0 = bet[1] + (0 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 0 vaxx
	pk0 = 1 - 1/(1+exp(logitpk0))
    ypred0 = rbinom(K,nk,pk0)  # new: predicted deaths for 0 vaccination
   
    logitpk10 = bet[1] + (0.1 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 10 vaxx
    pk10 = 1 - 1/(1+exp(logitpk10))
    ypred10 = rbinom(K,nk,pk10)  # new: predicted deaths for 10 vaccination
    
    logitpk20 = bet[1] + (0.2 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 20 vaxx
    pk20 = 1 - 1/(1+exp(logitpk20))
    ypred20 = rbinom(K,nk,pk20)  # new: predicted deaths for 20 vaccination
    
    logitpk30 = bet[1] + (0.3 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 30 vaxx
    pk30 = 1 - 1/(1+exp(logitpk30))
    ypred30 = rbinom(K,nk,pk30)  # new: predicted deaths for 30 vaccination
    
    logitpk40 = bet[1] + (0.4 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 40 vaxx
    pk40 = 1 - 1/(1+exp(logitpk40))
    ypred40 = rbinom(K,nk,pk40)  # new: predicted deaths for 40 vaccination
    
    logitpk50 = bet[1] + (0.5 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 50 vaxx
    pk50 = 1 - 1/(1+exp(logitpk50))
    ypred50 = rbinom(K,nk,pk50)  # new: predicted deaths for 50 vaccination
    
    logitpk60 = bet[1] + (0.6 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 60 vaxx
    pk60 = 1 - 1/(1+exp(logitpk60))
    ypred60 = rbinom(K,nk,pk60)  # new: predicted deaths for 60 vaccination
    
    logitpk70 = bet[1] + (0.7 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 70 vaxx
    pk70 = 1 - 1/(1+exp(logitpk70))
    ypred70 = rbinom(K,nk,pk70)  # new: predicted deaths for 70 vaccination
    
    logitpk80 = bet[1] + (0.8 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 80 vaxx
    pk80 = 1 - 1/(1+exp(logitpk80))
    ypred80 = rbinom(K,nk,pk80)  # new: predicted deaths for 80 vaccination
    
    logitpk90 = bet[1] + (0.9 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 90 vaxx
    pk90 = 1 - 1/(1+exp(logitpk90))
    ypred90 = rbinom(K,nk,pk90)  # new: predicted deaths for 90 vaccination
    
    logitpk100 = bet[1] + (1 - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 100 vaxx
    pk100 = 1 - 1/(1+exp(logitpk100))
    ypred100 = rbinom(K,nk,pk100)  # new: predicted deaths for 10 vaccination
    
    }   # end skip
   predhist[mc,] = ypred
   
   pred0hist[mc,] = ypred0
   pred10hist[mc,] = ypred10
   pred20hist[mc,] = ypred20
   pred30hist[mc,] = ypred30
   pred40hist[mc,] = ypred40
   pred50hist[mc,] = ypred50
   pred60hist[mc,] = ypred60
   pred70hist[mc,] = ypred70
   pred80hist[mc,] = ypred80
   pred90hist[mc,] = ypred90
   pred100hist[mc,] = ypred100
   
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
      
   
   stmean = apply(sthist,1,mean)
    plot(stmean)
   
   yhat = apply(predhist,2, "median")
   plot(yhat,Y,log="xy")
   
   y0hat = apply(pred0hist,2, "median")
   sum(y0hat)
   plot(y0hat,Y,log="xy")
   
   TotD0hist = apply(pred0hist,1, "sum")
   plot(TotD0hist)
   (CI.TotD0 = quantile(TotD0hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD0hist[rnge])
   
   TotD10hist = apply(pred10hist,1, "sum")
   plot(TotD10hist)
   (CI.TotD10 = quantile(TotD10hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD10hist[rnge])
   
   TotD20hist = apply(pred20hist,1, "sum")
   plot(TotD20hist)
   (CI.TotD20 = quantile(TotD20hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD20hist[rnge])
   
   TotD30hist = apply(pred30hist,1, "sum")
   plot(TotD30hist)
   (CI.TotD30 = quantile(TotD30hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD30hist[rnge])
   
   TotD40hist = apply(pred40hist,1, "sum")
   plot(TotD40hist)
   (CI.TotD40 = quantile(TotD40hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD40hist[rnge])
   
   TotD50hist = apply(pred50hist,1, "sum")
   plot(TotD50hist)
   rnge = burnin:M 
   (CI.TotD50 = quantile(TotD50hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD50hist[rnge])#  one run, 5000x10 
  
   TotD60hist = apply(pred60hist,1, "sum")
   plot(TotD60hist)
   (CI.TotD60 = quantile(TotD60hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD60hist[rnge])
   
   TotD70hist = apply(pred70hist,1, "sum")
   plot(TotD70hist)
   (CI.TotD70 = quantile(TotD70hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD70hist[rnge])
   
   TotD80hist = apply(pred80hist,1, "sum")
   plot(TotD80hist)
   (CI.TotD80 = quantile(TotD80hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD80hist[rnge])
   
   TotD90hist = apply(pred90hist,1, "sum")
   plot(TotD90hist)
   (CI.TotD90 = quantile(TotD90hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD90hist[rnge])
   
   TotD100hist = apply(pred100hist,1, "sum")
   plot(TotD100hist)
   (CI.TotD100 = quantile(TotD100hist[rnge], probs = c(0.025,0.975)) )
   mean(TotD100hist[rnge])
  
   ######## Making graphs for vaccination and death projection
   
   TotDhist = 1.1476* c (mean(TotD0hist[rnge]), mean(TotD10hist[rnge]), mean(TotD20hist[rnge]),
                 mean(TotD30hist[rnge]), mean(TotD40hist[rnge]), mean(TotD50hist[rnge]),
                 mean(TotD60hist[rnge]), mean(TotD70hist[rnge]), mean(TotD80hist[rnge]),
                 mean(TotD90hist[rnge]), mean(TotD100hist[rnge]))
   
   CI.TotD = 1.1476* t(cbind(CI.TotD0, CI.TotD10, CI.TotD20, CI.TotD30, CI.TotD40,
                     CI.TotD50, CI.TotD60, CI.TotD70, CI.TotD80, CI.TotD90, CI.TotD100))
   
   
   xaxis = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
   graph_data <- data.frame (xaxis, TotDhist, CI.TotD)
   
   ggplot(graph_data,
          aes(x = xaxis,
              y = TotDhist)) +
     geom_smooth(aes(ymin = X2.5., 
                     ymax = X97.5.),
                 stat = "identity",
                 color = "red") +
     geom_point(color = "blue")+
     ylim(0,1000000) +
     scale_x_continuous(breaks = xaxis) +
     xlab("Vaccination coverage (%)") +
     ylab("Projected total death")
   
   
   
   #  2.5%  97.5%   652004 784202     mean(TotD0hist[rnge])  718920.6   time ~43 min
   #  one run, 5000x20
  #    2.5%  97.5% 662015 780258        mean 719867.2
  #
   mean(TotD0hist[rnge]) - sum(Y)
   # 435285.2   "total lives saved"
   
   
    plot(Y,y0hat, asp=1, xlim = c(-10,15000), ylim = c(-10,15000))
   lines(c(1,10000),c(1,10000), col="red")


 
 
 temporary:   (still no convergence yet)
 
               2.5%   97.5%
CI.rho       0.7785  0.9844
CI.tau       0.7553  0.8463
CI.nu        0.2540  0.6311
(Intercept) -2.8422 -2.5070
X           -0.3664 -0.0681


cases Aug-Oct only (about 20 hr total runs)

               2.5%   97.5%
CI.rho       0.3720  0.5655
CI.tau       0.5675  0.6624
CI.nu        0.4354  0.6633
(Intercept) -3.3548 -3.2821
X           -0.9933 -0.6011









