
 ## this is a Gibbs sampler for Leroux model coded by OM
 #     
 
  require(MASS)
  require(Matrix)
  require(beepr)
 
 setwd("C:/localH/research/Awal/codeFeb15")

  cnew3 = read.csv( "cnew3.csv")
  Cmat3a = readMM(file='Cmat3a.txt')
  co.mat3a = as.matrix(Cmat3a) + 0
  
   totDeth = cnew3$dethAug + cnew3$dethSep + cnew3$dethOct + cnew3$dethNov + cnew3$dethDec + cnew3$dethJan + cnew3$dethFeb
   # Deaths are given per 100,000
   # totCase = cnew3$case.pAug + cnew3$case.pSep + cnew3$case.pOct + cnew3$case.pNov + cnew3$case.pDec + cnew3$case.pJan
      
	Ndeath = round(totDeth * cnew3$TotalPop/1e5)  # back to integers 
  	Y = Ndeath
	 
    b0s = log(sum(Y)/sum(cnew3$TotalPop))   #  this is an approximate value to start b0 with  
    st = cnew3$state
		
   ## prep the state vars for MCMC
	  
	  st.names = levels(as.factor(cnew3$state))
	  nst = length(st.names)
	  mem = NULL 
	  for (i in 1:nst){
	    mem[[i]] = which(cnew3$state == st.names[i])	  	  
	  }
 
   # list of neighbors 
 
   K = dim(Cmat3a)[1]
   Nlist = NULL 
   nneighb = numeric(K)
   for (k in 1:K){
      Nlist[[k]] = which(Cmat3a[k,] == 1)	
  	  nneighb[k] = length(Nlist[[k]])
	}  

	logfcp.rho = function(r){                     # FCPs 
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
	
	logfcp.phik = function(phik,k, mu1, d1){   # can't do vector form
	   log1e = log(1 + exp(Yhat[k] + phik + Ste[k]))
	   add1 = Y[k]*phik - nk[k]*log1e 
	   add2 = - d1/2/tausq*(phik - mu1)^2
	   add3 = dnorm(phik,0,2, log=TRUE)
	  return(add1 + add2 + add3)
	}
	
	logfcp.st = function(stp,j){   # vector form
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
	  
 vaxx = cnew3$vacJul/100    # rescaled to between 0 and 1
 meanvaxx = mean(vaxx)
 X = vaxx - meanvaxx
 nk = cnew3$TotalPop

  lm1 = lm(Y ~ X)      # "state" effect will be dealt with separately; using lm here just to get the Xmat
  Xmat = model.matrix(lm1)      
  XtXi = solve(t(Xmat) %*% Xmat)
  XtXiXt = XtXi %*% t(Xmat)
  bet = c(b0s,-1)    
        
  tau = 0.9; tausq = tau^2 	   # initializing model variables
  nu = 0.5
  rho = 0.5 
  nu0 = 0.1;   df.nu = 2
  tau0 = 1; df.tau = 2   # priors on variances
  rho.step = 0.1
  phi.step = 0.5   # same step for all k, for now
  st.step  = 0.15    
    beta.step = 0.2   
  phi = rep(0,K)
  st = rep(0,nst); Ste = numeric(K)
  ypred = matrix(0,K,11)    # to store simulated counterfactual counts
  vaccpct = (0:10)/10        # vacc.proportion to use in ypred
  
  M = 2000      # M*Nskip = 1000 takes about 47 sec;    a quality run would be M = 2000, Nskip = 100
  Nskip = 100     
  burnin = M*0.2;   M10 = M %/% 10
  
  phihist = matrix(0,M,K); predhist = matrix(0,M,K); sthist = matrix(0,M,nst)
  ypredhist = array(0,c(M,K,11))
  nuhist = numeric(M);  tauhist = numeric(M)                            # storage vars for MCMC
  rhohist = numeric(M); betahist = matrix(0,M,dim(Xmat)[2])
  acc.ctr.rho = 0;  acc.ctr.phi = 0; acc.ctr.beta = 0; acc.ctr.st = 0     
  
  ptm = proc.time()
  
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
   Yhat = Xmat %*% bet    
	  
    # sampling state effects   
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
	 
    # sampling phik 
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
   
   # sampling rho   
	stepr = rho.step
    rhop = rho + runif(1,-stepr, stepr)
	
	if ( (rhop < 0.9999) & (rhop > 0)){
		acc.prob = exp(min(logfcp.rho(rhop) - logfcp.rho(rho),0))
		if (runif(1) < acc.prob){
		  rho = rhop
		  acc.ctr.rho = acc.ctr.rho + 1	
		}
    }  
 
    }   # end skip
	
  for (i in 1:11){  # run counterfactual scenarios 
   	logitpk = bet[1] + (vaccpct[i] - meanvaxx)*bet[2] + phi + Ste   # projected log odds for 0 vaxx
	pk = 1 - 1/(1+exp(logitpk))
    ypred[,i] = rbinom(K,nk,pk)  # new: predicted deaths for 0 vaccination
   }  		
		
   ypredhist[mc,,] = ypred		    
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
   MC.bet = NULL
   for (i in 1:dim(Xmat)[2]){
      mean.bet = mean(betahist[rnge,i])
      MC.bet = rbind(MC.bet, 
	        c(mean.bet,quantile(betahist[rnge,i], probs = c(0.025,0.975))) 
		 )   
   }
   mean.rho = mean(rhohist[rnge])
   mean.tau = mean(tauhist[rnge])
   mean.nu = mean(nuhist[rnge])
   
   OUT = rbind(c(mean.rho,CI.rho), c(mean.tau, CI.tau),  c(mean.nu, CI.nu), MC.bet)
   row.names(OUT) = c("rho","tau","nu", colnames(Xmat))
   
   print(round(OUT, 4))   
   
   TotDhist = matrix(0,M,11)
   CI.TotD = matrix(0,2,11)
   for (i in 1:11){ 
        TotDhist[,i] = apply(ypredhist[,,i],1, "sum")
        CI.TotD[,i] = quantile(TotDhist[rnge,i], probs = c(0.025,0.975)) 
      print(c(vaccpct[i],mean(TotDhist[rnge,i]), CI.TotD[,i]))
	 }
   
  
   ######## Making graphs for vaccination and death projection
    
   PopFactor = 1.1476	
   PostMean = PopFactor*apply(TotDhist[rnge,],2,mean)      
   
   
   xaxis = seq(0,100,by=10)
   graph_data <- data.frame (xaxis, PostMean, X2.5. = PopFactor*CI.TotD[1,], X97.5. = PopFactor*CI.TotD[2,] )
   
   require(ggplot2)
   ggplot(graph_data,
          aes(x = xaxis,
              y = PostMean)) +
     geom_smooth(aes(ymin = X2.5., 
                     ymax = X97.5.),
                 stat = "identity",
                 color = "red") +
     geom_point(color = "blue")+
     ylim(0,1000000) +
     scale_x_continuous(breaks = xaxis) +
     xlab("Vaccination coverage (%)") +
     ylab("Projected total death")
   
  
    PopFactor*(mean(TotDhist[rnge,1]) - sum(Y))   # total lives saved
   
 
   # xaxis PostMean    X2.5.    X97.5.
# 1      0 995854.7 891489.0 1104260.9
# 2     10 809696.4 737195.3  901680.8
# 3     20 658411.4 601578.8  718577.8
# 4     30 533786.2 494922.0  581735.7
# 5     40 434932.8 404298.3  476054.3
# 6     50 354091.6 327107.3  385677.4
# 7     60 287838.3 265429.6  316609.1
# 8     70 233621.1 214481.8  259242.8
# 9     80 190133.7 171505.4  212118.9
# 10    90 154407.4 138184.8  172849.2
# 11   100 125798.6 110310.8  143277.9
 


