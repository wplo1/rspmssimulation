library(Surrogate)
library(doParallel); library(foreach)
finalperf <- function(resdata,truers,truepms){
  resdf.df <- data.frame()
  for(i in 1:nrow(resdata)){resdf.df<-rbind(resdf.df, resdata[i,]$resulX)}
  resdf.df[,2]<-as.numeric(resdf.df[,2])
  resdf.df[,3]<-as.numeric(resdf.df[,3])
  resdf.df[,4]<-as.numeric(resdf.df[,4])
  resdf.df[,5]<-as.numeric(resdf.df[,5])
  resdf.df[,6]<-as.numeric(resdf.df[,6])
  resdf.df[,7]<-as.numeric(resdf.df[,7])
  resdf.df[,8]<-as.numeric(resdf.df[,8])
  resdf.df[,9]<-as.numeric(resdf.df[,9])
  resdf.df[,10]<-as.numeric(resdf.df[,10])
  resdf.df[,11]<-as.numeric(resdf.df[,11])
  resdf.df[,12]<-as.numeric(resdf.df[,12])
  resdf.df[,13]<-as.numeric(resdf.df[,13])
  resdf.df<-transform(resdf.df,covrs=ifelse(lowerRs<=truers&upperRs>=truers,1,0),covpms=ifelse(lowerPMS<=truepms&upperPMS>=truepms,1,0))
  perfdf<-data.frame(method=c("parast","kuroki"),
                     bias=c(mean(resdf.df$Rsemp)-truers,mean(resdf.df$PMSemp)-truepms),
                     empSE=c(sqrt(var(resdf.df$Rsemp)),sqrt(var(resdf.df$PMSemp))),
                     avgModSE=c(sqrt(mean(resdf.df$varRsemp)),sqrt(mean(resdf.df$varPMSemp))),
                     coverage=c(sum(resdf.df$covrs)/length(resdf.df$covrs),sum(resdf.df$covpms)/length(resdf.df$covpms)),
                     extrapolate=c(median(resdf.df$TCprop),median(resdf.df$CTprop)))
  perfdf<-data.frame(perfdf,MSE=c(mean((resdf.df$Rsemp-truers)^2),mean((resdf.df$PMSemp-truepms)^2)))
  perfdf2<-perfdf[,c(1,2,3,7,4,5,6)]
  montedf<-data.frame(method=c("parast_mc","kuroki_mc"),
                      bias=c(perfdf$empSE[1]/sqrt(nsim),perfdf$empSE[2]/sqrt(nsim)),
                      empSE=perfdf$empSE/sqrt(2*nsim-2),
                      MSE=c(sqrt(mean(((resdf.df$Rsemp-truers)^2-perfdf$MSE[1])^2)/(nsim-1)),sqrt(mean(((resdf.df$PMSemp-truepms)^2-perfdf$MSE[2])^2)/(nsim-1))),
                      avgModSE=c(sqrt(var(resdf.df$varRsemp)/(4*nsim*perfdf$avgModSE[1]^2)),sqrt(var(resdf.df$varPMSemp)/(4*nsim*perfdf$avgModSE[2]^2))),
                      coverage=sqrt(perfdf$coverage*(1-perfdf$coverage)/nsim),
                      extrapolate=c(IQR(resdf.df$TCprop),IQR(resdf.df$CTprop)))
  return(list(rbind(perfdf2,montedf)))
}

#### XNormalOverlap200LL1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,ttemp,ccemp,tcemp,ctemp,CPMemp,CPIemp,NPSemp))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-10;bYT<-0.05;aYC<-5; bYC<-0.97;
  ts1<-rlnorm(200,mean=muST,sd=sigmaST);ts0<-rlnorm(200,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-10;bYT<-0.05;aYC<-5; bYC<-0.97;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap200LL1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap200LL3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-10;bYT<-0.05;aYC<-5; bYC<-0.97;
  ts1<-rlnorm(300,mean=muST,sd=sigmaST);ts0<-rlnorm(100,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-10;bYT<-0.05;aYC<-5; bYC<-0.97;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap200LL3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap200ML1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-0.5;bYT<-0.6;aYC<-2; bYC<-1;
  ts1<-rlnorm(200,mean=muST,sd=sigmaST);ts0<-rlnorm(200,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-0.5;bYT<-0.6;aYC<-2; bYC<-1;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap200ML1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap200ML3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-0.5;bYT<-0.6;aYC<-2; bYC<-1;
  ts1<-rlnorm(300,mean=muST,sd=sigmaST);ts0<-rlnorm(100,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-0.5;bYT<-0.6;aYC<-2; bYC<-1;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap200ML3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time





#### XNormalOverlap200MM1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-10;bYT<-0.2;aYC<-1; bYC<-0.6;
  ts1<-rlnorm(200,mean=muST,sd=sigmaST);ts0<-rlnorm(200,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-10;bYT<-0.2;aYC<-1; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap200MM1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap200MM3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-10;bYT<-0.2;aYC<-1; bYC<-0.6;
  ts1<-rlnorm(300,mean=muST,sd=sigmaST);ts0<-rlnorm(100,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-10;bYT<-0.2;aYC<-1; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap200MM3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap200MH1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,ttemp,ccemp,tcemp,ctemp,CPMemp,CPIemp,NPSemp))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-2;bYT<-0.4;aYC<-0.2; bYC<-0.6;
  ts1<-rlnorm(200,mean=muST,sd=sigmaST);ts0<-rlnorm(200,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-2;bYT<-0.4;aYC<-0.2; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap200MH1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap200MH3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-2;bYT<-0.4;aYC<-0.2; bYC<-0.6;
  ts1<-rlnorm(300,mean=muST,sd=sigmaST);ts0<-rlnorm(100,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-2;bYT<-0.4;aYC<-0.2; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap200MH3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap200HM1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-1;bYT<-0.5;aYC<-10; bYC<-0.6;
  ts1<-rlnorm(200,mean=muST,sd=sigmaST);ts0<-rlnorm(200,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-1;bYT<-0.5;aYC<-10; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap200HM1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap200HM3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-1;bYT<-0.5;aYC<-10; bYC<-0.6;
  ts1<-rlnorm(300,mean=muST,sd=sigmaST);ts0<-rlnorm(100,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-1;bYT<-0.5;aYC<-10; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap200HM3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap200HH1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-1;bYT<-0.7;aYC<-10; bYC<-0.6;
  ts1<-rlnorm(200,mean=muST,sd=sigmaST);ts0<-rlnorm(200,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-1;bYT<-0.7;aYC<-10; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap200HH1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap200HH3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-1;bYT<-0.7;aYC<-10; bYC<-0.6;
  ts1<-rlnorm(300,mean=muST,sd=sigmaST);ts0<-rlnorm(100,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-1;bYT<-0.7;aYC<-10; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap200HH3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time


#### XNormalOverlap1000LL1 ####
cs <- makeCluster(1)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-10;bYT<-0.05;aYC<-5; bYC<-0.97;
  ts1<-rlnorm(1000,mean=muST,sd=sigmaST);ts0<-rlnorm(1000,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-10;bYT<-0.05;aYC<-5; bYC<-0.97;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap1000LL1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap1000LL3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-10;bYT<-0.05;aYC<-5; bYC<-0.97;
  ts1<-rlnorm(1500,mean=muST,sd=sigmaST);ts0<-rlnorm(500,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-10;bYT<-0.05;aYC<-5; bYC<-0.97;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap1000LL3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time


#### XNormalOverlap1000ML1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-0.5;bYT<-0.6;aYC<-2; bYC<-1;
  ts1<-rlnorm(1000,mean=muST,sd=sigmaST);ts0<-rlnorm(1000,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-0.5;bYT<-0.6;aYC<-2; bYC<-1;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap1000ML1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap1000ML3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-0.5;bYT<-0.6;aYC<-2; bYC<-1;
  ts1<-rlnorm(1500,mean=muST,sd=sigmaST);ts0<-rlnorm(500,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-0.5;bYT<-0.6;aYC<-2; bYC<-1;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap1000ML3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time





#### XNormalOverlap1000MM1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-10;bYT<-0.2;aYC<-1; bYC<-0.6;
  ts1<-rlnorm(1000,mean=muST,sd=sigmaST);ts0<-rlnorm(1000,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-10;bYT<-0.2;aYC<-1; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap1000MM1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap1000MM3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-10;bYT<-0.2;aYC<-1; bYC<-0.6;
  ts1<-rlnorm(1500,mean=muST,sd=sigmaST);ts0<-rlnorm(500,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-10;bYT<-0.2;aYC<-1; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap1000MM3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap1000MH1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-2;bYT<-0.4;aYC<-0.2; bYC<-0.6;
  ts1<-rlnorm(1000,mean=muST,sd=sigmaST);ts0<-rlnorm(1000,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-2;bYT<-0.4;aYC<-0.2; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap1000MH1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap1000MH3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-2;bYT<-0.4;aYC<-0.2; bYC<-0.6;
  ts1<-rlnorm(1500,mean=muST,sd=sigmaST);ts0<-rlnorm(500,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-2;bYT<-0.4;aYC<-0.2; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap1000MH3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap1000HM1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-1;bYT<-0.5;aYC<-10; bYC<-0.6;
  ts1<-rlnorm(1000,mean=muST,sd=sigmaST);ts0<-rlnorm(1000,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-1;bYT<-0.5;aYC<-10; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap1000HM1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap1000HM3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-1;bYT<-0.5;aYC<-10; bYC<-0.6;
  ts1<-rlnorm(1500,mean=muST,sd=sigmaST);ts0<-rlnorm(500,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-1;bYT<-0.5;aYC<-10; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap1000HM3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap1000HH1 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-1;bYT<-0.7;aYC<-10; bYC<-0.6;
  ts1<-rlnorm(1000,mean=muST,sd=sigmaST);ts0<-rlnorm(1000,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-1;bYT<-0.7;aYC<-10; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
XNormalOverlap1000HH1<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time

#### XNormalOverlap1000HH3 ####
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()
nsim=625
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  set.seed(i)
  WPLsimfast <- function(sone,szero,yone,yzero){
    # Miscellaneous
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34) 
    # Which
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    # Estimate
    ttemp <- mean(yone)
    ccemp <- mean(yzero)
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    extraTC<-which(is.nan(tcpre))
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    tcemp <- mean(as.numeric(tcpre))
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    extraCT<-which(is.nan(ctpre))
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    ctemp <- mean(as.numeric(ctpre))
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    # Variance
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    extracolTC<-which(is.na(colSums(TCperturb3)))
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;
  aYT<-1;bYT<-0.7;aYC<-10; bYC<-0.6;
  ts1<-rlnorm(1500,mean=muST,sd=sigmaST);ts0<-rlnorm(500,mean=muSC,sd=sigmaSC)
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-exp(bYC))^3-log(s))))};
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-exp(bYT))^3-log(s))))};
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}
muST<-0.8; sigmaST<-0.2; muSC<-0.5; sigmaSC<-0.2;aYT<-1;bYT<-0.7;aYC<-10; bYC<-0.6;
YTYT<-function(x){1/(1+exp(-aYT*(x-exp(bYT))^3-log(x)))};YCYC<-function(x){1/(1+exp(-aYC*(x-exp(bYC))^3-log(x)))}
STST<-function(x){dlnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dlnorm(x,mean=muSC,sd=sigmaSC)}
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},0,Inf)
TTTT<-integrate(function(x){YTYT(x)*STST(x)},0,Inf)
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},0,Inf)
CTCT<-integrate(function(x){YCYC(x)*STST(x)},0,Inf)
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(3))
XNormalOverlap1000HH3<-finalperfdf[[1]]
end.time<-Sys.time()
end.time-start.time