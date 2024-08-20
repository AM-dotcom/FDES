

#================================ Saddlepoint functions ================================


KHat <- function(Tau,ThetaV,param){log(mean(exp(
  crossprod(Tau,PPsiW(x = ThetaV,parms = param))
)))}

MeanPsi <- function(t1,parms){
  Px <- parms$data
  T1 <- matrix(t1,nrow=length(t1),ncol=dim(Px)[2])
  return(Px - T1)}

MeanK <- function(TauM,t2,parms){ log(mean(exp(crossprod(TauM,MeanPsi(t2,parms))))) }

MeanDK <- function(TauM,parms){
  P1 <- MeanPsi(parms$t2,parms)
  E1 <- exp(crossprod(TauM,P1))
  apply( P1 * rbind(E1,E1), MARGIN = 1, FUN = mean)
}

outerProd <- function(v){crossprod(t(v),v)}

MeanSigma <- function(TauM,t2,parms){
  d1 <- length(t2)
  nPsi <- dim(parms$data)[2]
  P1 <- MeanPsi(t2,parms=parms)
  E1 <- exp(-MeanK(TauM,t2,parms=parms))
  E2 <- exp(crossprod(TauM,P1))
  EMS <- matrix(E2,nPsi,d1^2)
  
  SigHatM <- EMS*t(apply(P1,MARGIN = 2,FUN = outerProd))
  SigHat1 <- apply(SigHatM,MARGIN = 2, FUN = mean)
  SigHat <- E1*matrix(SigHat1,nrow = d1)
  return(SigHat)
}

