source('~/Documents/Rstuff/DynamicFactorAnalysis_TMB_UsingUnstructuredCorr_Phased.R')

RsimUNC<-as.matrix(read.csv("~/Documents/Rstuff/BigCovarianceMatrix.csv",header=F))

Zscore<-function(x){
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
}

#These are the loadings for the simulation
Zsim<-matrix(runif(nrow(RsimUNC)*2,0,1),ncol=2)#c(0.3,0.5,0.8)
Zsim[,1]<-sort(Zsim[,1])
Zsim[,2]<-rev(sort(Zsim[,2]))

#Simulated randomwalk common trend
uSim<-matrix(rep(rep(NA,100),2),nrow=2)
uSim[,1]<--3
for(i in 2:100){
  uSim[,i]<-uSim[,i-1]+rnorm(2,0,1)
}
uSim<-t(apply(uSim,1,FUN=Zscore))

plot(uSim)

#simulate a ts for air temp. I used an autocorrelated process but not random walk
airSim<-rep(NA,100)
airSim[1]<-0
for(i in 2:100){
  airSim[i]<-0.7*airSim[i-1]+rnorm(1,0,1)
}
airSim<-Zscore(airSim)

plot(airSim,type='l')

#AirTemperature loadings for simulation
Dsim<-runif(nrow(RsimUNC),0,1)#c(0.1,0.6,0.9)

library(MASS)

#R matrices for simulation. Testing all three common forms.
#RsimDE<-diag(0.1,3)
#RsimDUE<-diag(c(0.1,0.2,0.4),3)
#RsimUNC<-matrix(c(0.004328536,0.025618901,0.023676532,0.025618901,0.1516287,0.1405050,0.023676532,0.1405050,0.38084122),nrow=3,ncol=3)

#RsimUNC<-matrix(c(0.038364445,-0.001262234,-0.079377036,-0.001262234,0.067458194,0.013427599,-0.07937704,0.01342760,0.16810893),nrow=3,ncol=3)
#Simulate data series
simObs<- matrix(Zsim,ncol=2) %*% uSim + matrix(Dsim,ncol=1) %*% airSim + t(mvrnorm(100,rep(0,nrow(RsimUNC)),Sigma=RsimUNC))

simObs[sample.int(300,size=50)]<-NA

plot(simObs[1,],type='l',ylim=c(-4,4))
points(simObs[2,],type='l',col='red')
points(simObs[3,],type='l',col='blue')

library(MARSS)
#Fit the model using MARSS
marssFit<-MARSS(simObs,model=list(m=1,R='unconstrained'),covariates=airSim,form='dfa',control=list(maxit=500,MCInit=T))
marssPred<-coef(marssFit,type='matrix')$Z %*% marssFit$states + coef(marssFit,type='matrix')$D %*% airSim


par(mfrow=c(3,2),mar=c(2,2,1,1))
ScaleFacM<-max(Mod(marssFit$states[1,]))/max(Mod(uSim))
plot(uSim)
points(marssFit$states[1,]/ScaleFacM,type='l')

plot(marssFit$states[1,]/ScaleFacM~uSim)
abline(0,1)

coef(marssFit)$Z*ScaleFacM
Zsim

plot(coef(marssFit)$Z*ScaleFacM~Zsim,ylim=c(0,1),xlim=c(0,1))
abline(0,1)

coef(marssFit)$D
Dsim

plot(coef(marssFit)$D~Dsim,ylim=c(0,1),xlim=c(0,1))
abline(0,1)

plot(simObs[2,])
points(marssPred[2,],type='l')

plot(coef(marssFit)$R~RsimUNC[lower.tri(RsimUNC,diag=T)],ylim=c(min(RsimUNC),max(RsimUNC)),xlim=c(min(RsimUNC),max(RsimUNC)))
abline(0,1)


# Fit the model with TMB code
myFit<-runDFA(simObs,NumStates=2,ErrStruc='UNC',EstCovar=T,Covars=matrix(airSim,nrow=1))

par(mfrow=c(2,2),mar=c(4,4,1,1))
ScaleFac1<-max(Mod(myFit$Estimates$u[1,]))/max(Mod(uSim[2,]))
ScaleFac2<-max(Mod(myFit$Estimates$u[2,]))/max(Mod(uSim[1,]))
plot(uSim[1,],main='Shared Trend')
points(myFit$Estimates$u[2,]/ScaleFac1,type='l')
points(uSim[2,],col='red')
points(myFit$Estimates$u[1,]/ScaleFac2,type='l',col='red')
legend('top',legend=c('Simulated','ModelFit'),pch=c(1,NA),lty=c(NA,1))



#plot(myFit$Estimates$u[1,]/ScaleFac~uSim,xlab='si')
#abline(0,1)

myFit$Estimates$Z*ScaleFac
Zsim

plot(myFit$Estimates$Z[,1]*ScaleFac2~Zsim[,2],ylab='Estimated',xlab='Simulated',main='Loadings(Z)')
abline(0,1)
points(myFit$Estimates$Z[,2]*ScaleFac1~Zsim[,1])


myFit$Estimates$D
Dsim

plot(myFit$Estimates$D~Dsim,ylab='Estimated',xlab='Simulated',main='Covariates(D)')
abline(0,1)

plot(myFit$Estimates$R~RsimUNC,ylab='Estimated',xlab='Simulated',main='CovarianceMatrix(R)')
abline(0,1)



myFit$Estimates$R
coef(marssFit,type='matrix')$R
RsimUNC

 

