
#============================== Whittle's functions ============================


PSDARMA1d1 <- function(x,parms){
  CosV <- 2*cos(parms$freq)
  Scale <- parms$Sig^2
  return(Scale/((x[1]^2 - CosV*x[1] + 1)*(2-CosV)^{x[2]}))
}

PPSDARMA1d1 <- function(x,parms){
  CosV <- 2*cos(parms$freq)
  return(1/((x[1]^2 - CosV*x[1] + 1)*(2-CosV)^{x[2]}))
}

D_PPSDARMA1d1 <- function(x,parms){ #x=Parameters vector
  CosV <- 2*cos(parms$freq)
  Denom1 <- (x[1]^2 - CosV*x[1] + 1)
  D1 <- (CosV-2*x[1]) / (Denom1^2 * (2-CosV)^{x[2]})
  D2 <- -log(2-CosV)/(Denom1*(2-CosV)^{x[2]})
  return(rbind(D1,D2))
}

LD_PPSDARMA1d1 <- function(x,parms){ #x=Parameters vector
  CosV <- 2*cos(parms$freq) 
  D1 <- (CosV-2*x[1])/(x[1]^2 - CosV*x[1] + 1)
  D2 <- -log(2-CosV)
  return(rbind(D1,D2))
}

LDD1_PPSDARMA1d1 <- function(x,parms){ #x=Parameters vector
  CosV <- 2*cos(parms$freq)
  Denom1 <- (x[1]^2 - CosV*x[1] + 1)
  LDD1 <- ((2*x[1]-CosV)^2 / Denom1^2) - (2/Denom1)
  return(LDD1)
}

LDD2_PPSDARMA1d1 <- function(x,parms){ #x=Parameters vector
  LDD2 <- rep(0,length(parms$per))
  return(LDD2)
}

PPsiW <- function(x,parms){ #x=Parameters vector
  LD <- LD_PPSDARMA1d1(x,parms)
  CLD <- t(scale(t(LD),center=TRUE,scale = FALSE))
  SP <- parms$per / PPSDARMA1d1(x,parms)
  n1 <- length(SP)
  DW1 <- SP*CLD[1,]
  DW2 <- SP*CLD[2,]
  return(rbind(DW1,DW2))
}

DPPsiW <- function(x,parms){ #x=Parameters vector
  CosV <- 2*cos(parms$freq)
  Per <- parms$per
  PPSD <- PPSDARMA1d1(x,parms)
  PPSD2 <- PPSD^2
  DPPSD <- D_PPSDARMA1d1(x,parms)
  
  LD <- LD_PPSDARMA1d1(x,parms)
  CLD <- t(scale(t(LD),center=TRUE,scale = FALSE))
  
  LDD1 <- LDD1_PPSDARMA1d1(x,parms)
  CLDD1 <- LDD1-mean(LDD1)
  
  M11 <- (CLDD1*PPSD - DPPSD[1,]*CLD[1,]) / PPSD2
  M12 <-  - DPPSD[2,]*CLD[1,] / PPSD2
  M21 <- - DPPSD[1,]*CLD[2,] / PPSD2
  M22 <- - DPPSD[2,]*CLD[2,] / PPSD2
  
  return(cbind(Per*M11,Per*M12,
               Per*M21,Per*M22))
}

WEF <- function(x,parms){apply(PPsiW(x,parms),MARGIN = 1,FUN = mean)}
