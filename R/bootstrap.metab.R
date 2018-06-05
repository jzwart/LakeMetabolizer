# =============================================
# = Bootstrapping function for MLE metabolism =
# =============================================
bootstrap.metab <- function(guesses, pars, fn, do.obs, do.sat, k.gas, z.mix, irr, wtr, error.type='OE', n.boot=1000, ar1.resids, freq, ...){
  #Calculate DOHat given max likelihood parameter estimates
  predix <- predictDO(par=pars,do.obs=do.obs, do.sat=do.sat, k.gas=k.gas,
                      z.mix=z.mix, irr=irr, wtr=wtr, error.type = error.type) # returns MLE DO
  DOHat<-predix$DOHat
  res <- predix$res

  n.obs = length(do.obs)

  #If we are maintaining the ar1 component of the residuals,
  # we must estimate ar1 coeff and the ar1 residual standard deviation
  if(ar1.resids){
    ar1.lm    = lm(res[1:n.obs-1] ~ res[2:n.obs]-1)
    ar1.coeff = ar1.lm$coefficients
    ar1.sd    = sd(ar1.lm$residuals)
  }

  #Pre-allocate the result data frame
  result <- data.frame(boot.iter = 1:n.boot,gppCoeff = rep(NA,n.boot),rCoeff = rep(NA,n.boot),doInit = rep(NA,n.boot),nll = rep(NA,n.boot),
                       GPP=rep(NA,n.boot),R=rep(NA,n.boot),NEP=rep(NA,n.boot))

  for(i in 1:n.boot){
    #Randomize the residuals using one of two methods
    if(ar1.resids){ #residual randomization keeping the ar1 data structure
      simRes = rep(NA, n.obs)
      simRes[1] = sample(res[!is.na(res)],1)
      for(j in 2:n.obs){
        simRes[j] = ar1.coeff*simRes[j-1] + rnorm(n=1, sd=ar1.sd)
      }
    }else{ #Raw residual randomization
      #Randomize residuals without replacement
      simRes = sample(res[!is.na(res)], length(res), replace=FALSE)
    }

    do.sim = DOHat + simRes

    #Run metab model again with new simulated DO signal
    if(error.type=='OE'){
      boot.temp <- optim(guesses, fn=mleNllOE, do.obs=do.sim, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr,logged=logged)
      pars0 <- boot.temp$par

      if(logged){
        pars <- c("gppCoeff"=exp(pars0[1]), "rCoeff"=-exp(pars0[2]), "Q"=exp(pars0[3]), "nll"=boot.temp$value, "doInit"=exp(pars0[4]))
      }else{
        pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]), "nll"=boot.temp$value, "doInit"=pars0[4])
      }

      GPP <- mean(pars[1]*irr, na.rm=TRUE) * freq
      R <- mean(pars[2]*log(wtr), na.rm=TRUE) * freq
      NEP<-GPP+R
      # store result of current iteration in output data frame
      result[i,2:3]<-pars[1:2]
      result[i,4]<-pars[5]
      result[i,5]<-pars[4]
      result[i,6]<-GPP
      result[i,7]<-R
      result[i,8]<-NEP

    }else if(error.type=='PE'){
      boot.temp <- optim(guesses, fn=mleNllPE, do.obs=do.sim, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr,logged=logged)
      pars0 <- boot.temp$par

      if(logged){
        pars <- c("gppCoeff"=exp(pars0[1]), "rCoeff"=-exp(pars0[2]), "Q"=exp(pars0[3]), "nll"=boot.temp$value)
      }else{
        pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]), "nll"=boot.temp$value)
      }

      GPP <- mean(pars[1]*irr, na.rm=TRUE) * freq
      R <- mean(pars[2]*log(wtr), na.rm=TRUE) * freq
      NEP<-GPP+R
      # store result of current iteration in output data frame
      result[i,2:3]<-pars[1:2]
      result[i,4]<-do.sim[1]
      result[i,5]<-pars[4]
      result[i,6]<-GPP
      result[i,7]<-R
      result[i,8]<-NEP
    }
  }
  return(result)
}
