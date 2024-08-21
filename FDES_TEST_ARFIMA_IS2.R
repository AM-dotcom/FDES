
# rm(list = ls())

library(foreach)
library(doParallel)
library(doRNG)
library(rootSolve)
library(tictoc)
library(numDeriv)
library(longmemo)
library(mcmc)
library(cubature)
library(calculus)
library(MASS)
library(mvtnorm)
library(fda)
source('WhittleF3_1d0.R')
source('SaddlePointF4_1d0.R')

#=============================== Simulate density ==============================#

N <- 1e3      # sample size
Phi <- 0      # autoregressive parameter
d <- 0.25     # long-memory parameter
Sig <- 1      # standard deviation of the innovations
freqs <- 2*pi*(1:floor((N-1)/2))/N
s1 <- 1e4     # set MC sample size
s <- 210
d1 <- 2
acc <- 1e4
gridCeil <- 30
MCn <- 2e3    # set the sample size for importance sampling
Innov <- function(n1){simARMA0(n = n1, H = d + 0.5)} 
# define long-memory Gaussian innovations

TvalF <- function(param1){optim(par = c(0,0),fn = KHat,
                                ThetaV = c(Phi,d),
                                param = param1)}
# define the test statistic

tic()
set.seed(2444)
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)
List1 <- foreach(j = 1:s1, .packages = c('rootSolve','longmemo','MASS','calculus'),
                 .errorhandling = 'pass') %dorng% {
                   
   X <- arima.sim(n = N, model = list(ar = Phi), rand.gen = Innov, n.start = 1e3)
   
   Per1 <- Mod(fft(X))^2
   Per2 <- Per1[2:(floor((N-1)/2)+1)]/N
   NPsi <- length(Per2)
   
   test <- TvalF(list(per=Per2,freq=freqs,Sig=Sig))
   Tval <- -2*N*test$value
   list(Tval=Tval)
   
 }
stopCluster(cl)
toc()

TVAL <- unlist(List1)
nC1 <- length(TVAL)

TVAL <- unlist(List1)


#=============================== replicate s FDES ==============================#

tic()
set.seed(2444)
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)
List <- foreach(j = 1:s,
                .packages = c('rootSolve','longmemo','MASS','calculus'),
                .errorhandling="pass") %dorng% {
   
   X <- arima.sim(n = N, model = list(ar = Phi), rand.gen = Innov, n.start = 1e3)
   
   Per1 <- Mod(fft(X))^2
   Per2 <- Per1[2:(floor((N-1)/2)+1)]/N
   NPsi <- length(Per2)
   
   WhitHat <- multiroot(f = WEF, parms=list(per=Per2,freq=freqs),
                        start = c(Phi,d))$root
   
   #=========== Empirical saddlepoint for the mean of PsiX0:
   
   PsiX <- list(data = PPsiW(x = WhitHat, parms = list(per=Per2,freq=freqs)))
   PsiX0 <- list(data = PPsiW(x = c(Phi,d), parms = list(per=Per2,freq=freqs)))

   d1 <- dim(PsiX$data)[1]
   
   SN <- MeanSigma(rep(0,d1),rep(0,d1),PsiX)/NPsi
   invSN <- solve(SN)
   SNscale <- det(SN)^{1/2}
   
   MCp <- mvrnorm(n = MCn,mu = rep(0,d1),Sigma = SN)
  
   ToBeInt <- function(xBar){
     
     sp1 <- multiroot(f = MeanDK, start = rep(0,d1),
                      parms= append(PsiX,list(t2=xBar)))$root
     sp2 <- multiroot(f = MeanDK, start = rep(0,d1),
                      parms= append(PsiX0,list(t2=xBar)))$root
     
     MK <- MeanK(sp1,xBar,PsiX)
     
     invSN2 <- solve(MeanSigma(sp2,xBar,PsiX0))
     
     MK2 <- -tcrossprod(crossprod(xBar,invSN2),xBar)/2
     
     ChiSq1 <- tcrossprod(crossprod(xBar,invSN),xBar)
     
     E1 <- exp(NPsi*MK + ChiSq1/2)*SNscale
     
     S1 <- det(MeanSigma(sp1,xBar,PsiX))^(-1/2)
     
     return(c(-2*NPsi*MK2,E1*S1*NPsi^(3/2)))
   }
   
   SimV <- 1:MCn
   MKs <- SimV
   for(i in 1:MCn){
     getSP <- ToBeInt(MCp[i,])
     MKs[i] <- getSP[1]
     SimV[i] <- getSP[2]
     print(i)
   }
   
   Nconst <- mean(SimV)
   CDF1 <- 1:acc
   grid <- seq(0,gridCeil,length.out=acc)
   for (i in 1:acc){
     CDF1[i] <- mean(SimV*(MKs <= grid[i]))/Nconst
   }
   
   list(WhitHat = WhitHat, CDF = CDF1)
   
 }
stopCluster(cl)
toc()

CDFs <- matrix(NA,nrow=s,ncol=acc)
WhitHats <- matrix(NA,nrow=s,ncol=2)
for (i in 1:s){
  if(is.vector(List[[i]]$CDF)){
    CDFs[i,] <- List[[i]]$CDF}else{
    CDFs[i,] <- rep(NA,acc)
  }
  if(is.vector(List[[i]]$WhitHat)){
    WhitHats[i,] <- List[[i]]$WhitHat}else{
    WhitHats[i,] <- rep(NA,2)
    }
}
  
Mask <- (abs(WhitHats[,1])<0.999)&(abs(WhitHats[,2])<0.49999)&(!is.na(CDFs[,250]))
# Get rid of the estimates that did not converge


#=============================== Plot the results =============================#


grid <- seq(0,gridCeil,length.out=acc)
nC <- sum(Mask,na.rm=TRUE)

plot(CDFs[1,],x = grid,type='l')
for(i in 2:s){
  points(CDFs[i,],x = grid,type='l')
}

CleanCDF <- CDFs[Mask,]
plot(CleanCDF[1,],x = grid,type='l')
for(i in 2:nC){
  points(CleanCDF[i,],x = grid,type='l')
}

MyCDF <- CleanCDF[fbplot(fit = t(CleanCDF),x = grid,plot = FALSE)$medcurve,]

nC1 <- length(TVAL)
TCDF <- 1:acc
for (i in 1:acc){
  TCDF[i] <- sum(TVAL <= grid[i],na.rm=TRUE)/nC1
}

plot(MyCDF,x=grid,type='l')
points(TCDF,x=grid,type='l',col='blue3')
curve(expr = pchisq(q = x,df = 2),n = 200,add = TRUE,col='red3')

MyQF1 <- function(p1){grid[which.min(abs(MyCDF-p1))]}
MyQF <- Vectorize(MyQF1)
MyRG <- function(n1){
  MyQF(runif(n1))
}

MyRGV <- MyRG(1e4)

TestV <- TVAL
Margin1 <- 6
xMax <- round(quantile(TestV,0.99)+Margin1,0)

qqChisq <- qqplot(x = TestV, y = rchisq(1e4,df = 2),plot.it = FALSE,
                  xlim = c(0,xMax), ylim = c(0,xMax))
qqTrue <- qqplot(x = TestV, y = MyRGV,col = 'blue4',
                 main = 'ARFIMA(1,d,0) with and n = 1000',
                 #main = NULL,
                 xlab = 'True quantiles',
                 ylab = 'Approximate quantiles',axes=FALSE,cex=0.5,
                 xlim = c(0,xMax), ylim = c(0,xMax))

a1 <- axis(side = 2, at = c(0,5,10,15,20,28),labels = TRUE)
axis(side = 1, at=c(0,quantile(TestV,c(0.9,0.95,0.99)),xMax),
     labels= c(0,'q(0.9)','q(0.95)','q(0.99)',xMax))
points(qqChisq$x,qqChisq$y,col = 'red3',pch = 2,cex=0.5)
points(qqTrue$x,qqTrue$y,col = 'blue4',cex=0.5)
abline(v=quantile(TestV,c(0.9,0.95,0.99)),lty=3)
abline(0,1,lty=2)
legend(x = 0, y = xMax-1, legend = c('FDES','Chi-square','Identity'),
       col=c('blue4','red3','black'), pch = c(1,2,NA),
       lty = c(NA,NA,2),bty = 'n')

