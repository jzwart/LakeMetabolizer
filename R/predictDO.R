# ================================================
# = prediction function given modeled parameters =
# = returns DO predictions =
# ================================================
predictDO<-function(par, do.obs, do.sat, k.gas, z.mix, irr, wtr, error.type='OE', ...){
  match.arg(error.type,choices = c('OE','PE'))
  if(error.type=='PE'){
    c1 <- par[1] #PAR coeff
    c2 <- par[2] #log(Temp) coeff
    # Set first true value equal to first observation
    alpha <- rep(0, length(do.obs))
    alpha[1] <- do.obs[1]#Let's give this model some starting values
  }
  if(error.type=='OE'){
    c1 <- par[1] #PAR coeff
    c2 <- par[2] #log(Temp) coeff
      # Set first true value equal to first observation
    alpha <- rep(0, length(do.obs))
    alpha[1] <- par[3]#Let's give this model some starting values
  }

  # See KalmanDO_smooth.R comments for explanation of beta
  kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
  beta <- exp(-kz) # This beta is for using the differential equation form


  # R version of C loop
  for(i in 2:length(do.obs)){
    if(kz[i-1]!=0){
      a1 <- c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
      alpha[i] <- a1/kz[i-1] + -exp(-kz[i-1])*a1/kz[i-1] + beta[i-1]*alpha[i-1] # NOTE: beta==exp(-kz); kz=K/Zmix
    }else{
      alpha[i]<-alpha[i-1] + c1*irr[i-1]+c2*log(wtr[i-1])
    }
  }

  #Compare observed and predicted DO; calculate residuals
  #Exclude from calculation any cases where DOObs=NA
  if (any(is.na(do.obs)))
  {
    NAObs <- which(is.na(do.obs))
    res <- do.obs[-NAObs] - alpha[-NAObs]
  } else

  {
    res <- do.obs - alpha
  }

  nRes <- length(res)
  SSE <- sum(res^2)
  sigma2 <- SSE/nRes

  #Set up output structure
  dataOut <- list(DOHat=alpha,res=res)

  #Return dataOut
  return(dataOut)

}
