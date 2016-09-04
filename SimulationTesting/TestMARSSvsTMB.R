NaCols<-rep(NA,100)
for(i in 1:ncol(Z_StreamTemps13)){
  if(NA %in% Z_StreamTemps13[,i]){
    NaCols[i]<-i
  }
}
NaCols<-na.omit(NaCols)
noNA_Stream<-Z_StreamTemps13[,-NaCols]
noNA_Air<-Z_AirTemp13[-NaCols]
streamSet<-c(3,8,23)#c(3,8,10,15,23)
mod13<-runDFA(obs=noNA_Stream[streamSet,],NumStates=1,ErrStruc='UNC',EstCovar=TRUE,Covars=matrix(noNA_Air,nrow=1))
mod13_marss<-MARSS(noNA_Stream[streamSet,],model=list(m=1,R='unconstrained'),covariates=matrix(noNA_Air,nrow=1),form='dfa')

plot(mod13$Estimates$u[1,],type='l')
points(mod13_marss$states[1,],type='l',col='red')

plot(mod13$Estimates$Z~coef(mod13_marss)$Z)
abline(lm(mod13$Estimates$Z~coef(mod13_marss)$Z))

plot(mod13$Estimates$D~coef(mod13_marss)$D)
abline(lm(mod13$Estimates$D~coef(mod13_marss)$D))

mod13$Estimates$R;coef(mod13_marss,type='matrix')$R

MARSSfit<-coef(mod13_marss,type='matrix')$Z %*% mod13_marss$states[1,] + coef(mod13_marss,type='matrix')$D %*% matrix(noNA_Air,nrow=1)

par(mfrow=c(5,1),mar=c(2,2,1,1))
for(i in 1:5){
  Allp<-c(noNA_Stream[streamSet[i],],MARSSfit[i,],mod13$Fits[i,])
  plot(noNA_Stream[streamSet[i],],ylim=c(min(c(Allp)),max(Allp)))
  points(MARSSfit[i,],type='l',col='red',lwd=2)
  points(mod13$Fits[i,],type='l',col='blue',lwd=2)
  title(main=paste('MARSS:',round(summary(lm(noNA_Stream[streamSet[i],]~MARSSfit[i,]))$r.squared,2),' ::: TMB:',round(summary(lm(noNA_Stream[streamSet[i],]~mod13$Fits[i,]))$r.squared,2)))
}


#####

mod13_DUE<-runDFA(obs=noNA_Stream[streamSet,],NumStates=1,ErrStruc='UNC',EstCovar=TRUE,Covars=matrix(noNA_Air,nrow=1))
mod13_marss_DUE<-MARSS(noNA_Stream[streamSet,],model=list(m=1,R='unconstrained'),covariates=matrix(noNA_Air,nrow=1),form='dfa',control=list(maxit=2000))

plot(mod13_DUE$Estimates$u[1,]/max(Mod(mod13_DUE$Estimates$u[1,])),type='l')
points(mod13_marss_DUE$states[1,]/max(Mod(mod13_marss_DUE$states[1,])),type='l',col='red')

plot(mod13_DUE$Estimates$Z~coef(mod13_marss_DUE)$Z)
abline(lm(mod13_DUE$Estimates$Z~coef(mod13_marss_DUE)$Z))

plot(mod13_DUE$Estimates$D~coef(mod13_marss_DUE)$D)
abline(lm(mod13_DUE$Estimates$D~coef(mod13_marss_DUE)$D))

mod13_DUE$Estimates$R;coef(mod13_marss_DUE,type='matrix')$R

mod13_DUE$AIC
mod13_marss_DUE$AIC

MARSSfit<-coef(mod13_marss_DUE,type='matrix')$Z %*% mod13_marss_DUE$states[1,] + coef(mod13_marss_DUE,type='matrix')$D %*% matrix(noNA_Air,nrow=1)

par(mfrow=c(5,1),mar=c(2,2,1,1))
for(i in 1:5){
  Allp<-c(noNA_Stream[streamSet[i],],MARSSfit[i,],mod13_DUE$Fits[i,])
  plot(noNA_Stream[streamSet[i],],ylim=c(min(c(Allp)),max(Allp)))
  points(MARSSfit[i,],type='l',col='red',lwd=2)
  points(mod13_DUE$Fits[i,],type='l',col='blue',lwd=2)
  title(main=paste('MARSS:',round(summary(lm(noNA_Stream[streamSet[i],]~MARSSfit[i,]))$r.squared,2),' ::: TMB:',round(summary(lm(noNA_Stream[streamSet[i],]~mod13_DUE$Fits[i,]))$r.squared,2)))
}

resids<-noNA_Stream[streamSet,]-MARSSfit


cor(resids[1,],resids[2,])
cor(resids[1,],resids[3,])
cor(resids[1,],resids[4,])
cor(resids[1,],resids[5,])


