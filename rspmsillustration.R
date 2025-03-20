library(Surrogate)
library(doParallel); library(foreach)

# Function of resdata (vector of estimates, variances and extrapolation proportions described in resdf function below in line 64),
#             truers (true value of RS) and truepms (true value of PMS)
# Function generates table of performance measures (bias, empirical SE, MSE, average model SE, coverage) for each of RS and PMS with associated Monte Carlo SEs, 
# as well as median extrapolation proportions of TC and CT (not their perturbed versions) with interquartile ranges (IQRs) as final column
finalperf <- function(resdata,truers,truepms){
  resdf.df <- data.frame()
  # The 2nd-13th columns (total of 12 columns) of resdf.df below correspond to the vector of 12 elements of resdf (defined below in line 64)
  # This will be a dataframe, as the resdf function will be applied to each of the 625 simulations
  for(i in 1:nrow(resdata)){resdf.df<-rbind(resdf.df, resdata[i,]$resulX)}
  # Converting to a numerical environment
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
  
  # Computing coverage - what proportion of empirical 2.5%-97.5% quantiles of RS/PMS values contain their true values respectively
  resdf.df<-transform(resdf.df,covrs=ifelse(lowerRs<=truers&upperRs>=truers,1,0),covpms=ifelse(lowerPMS<=truepms&upperPMS>=truepms,1,0))
  
  # Point estimates of performance measures (bias, empirical SE, average model SE, coverage) and median extrapolation proportions of TC and CT (not their perturbed versions)
  perfdf<-data.frame(method=c("parast","kuroki"),
                     bias=c(mean(resdf.df$Rsemp)-truers,mean(resdf.df$PMSemp)-truepms),
                     empSE=c(sqrt(var(resdf.df$Rsemp)),sqrt(var(resdf.df$PMSemp))),
                     avgModSE=c(sqrt(mean(resdf.df$varRsemp)),sqrt(mean(resdf.df$varPMSemp))),
                     coverage=c(sum(resdf.df$covrs)/length(resdf.df$covrs),sum(resdf.df$covpms)/length(resdf.df$covpms)),
                     extrapolate=c(median(resdf.df$TCprop),median(resdf.df$CTprop)))
  # Point estimates of MSEs (mean squared error) appended to dataframe above
  perfdf<-data.frame(perfdf,MSE=c(mean((resdf.df$Rsemp-truers)^2),mean((resdf.df$PMSemp-truepms)^2)))
  perfdf2<-perfdf[,c(1,2,3,7,4,5,6)]
  
  # Monte Carlo SEs of performance measures, except final column with IQR indicating spread of extrapolation proportions of TC and CT (not their perturbed versions)
  montedf<-data.frame(method=c("parast_mc","kuroki_mc"),
                      bias=c(perfdf$empSE[1]/sqrt(nsim),perfdf$empSE[2]/sqrt(nsim)),
                      empSE=perfdf$empSE/sqrt(2*nsim-2),
                      MSE=c(sqrt(mean(((resdf.df$Rsemp-truers)^2-perfdf$MSE[1])^2)/(nsim-1)),sqrt(mean(((resdf.df$PMSemp-truepms)^2-perfdf$MSE[2])^2)/(nsim-1))),
                      avgModSE=c(sqrt(var(resdf.df$varRsemp)/(4*nsim*perfdf$avgModSE[1]^2)),sqrt(var(resdf.df$varPMSemp)/(4*nsim*perfdf$avgModSE[2]^2))),
                      coverage=sqrt(perfdf$coverage*(1-perfdf$coverage)/nsim),
                      extrapolate=c(IQR(resdf.df$TCprop),IQR(resdf.df$CTprop)))
  
  # Putting both point estimates and Monte Carlo SEs into one dataframe
  return(list(rbind(perfdf2,montedf)))
}

#### NormalSubset200LL1 ####

# Code runs simulations for the data-generating mechanism where both RM and RS (as defined in Proposition 2.3) are low, 
# in the case of a Normal surrogate endpoint and linear logit link for true endpoint, such that 
# the range of surrogate values in the control arm lies within the range of surrogate values in the treatment arm (subsetted supports), 
# under equal allocation (r=1) and a sample size of 200 per arm (total sample size 400).


# Six cores
cs <- makeCluster(6)
registerDoParallel(cs)

start.time<-Sys.time()

# Number of simulations
nsim=625

# resdf generates dataframe of 625 rows of RS and PMS estimates, variance, quantiles and extrapolation proportions
resdf <- foreach(i=1:nsim, .combine = 'rbind') %dopar% {
  library(Rmpfr);library(stats);library(dplyr);library(kit)
  
  # Reproducibility for each of the 625 simulations
  set.seed(i)
  
  # Function of vectors sone (treatment arm of surrogate), szero (control arm of surrogate), yone (treatment arm of true endpoint), yzero (control arm of true endpoint)
  # Function generates vectors of estimates, variances and extrapolation proportions
  WPLsimfast <- function(sone,szero,yone,yzero){
    
    # Allocation proportions (pc = control arm proportion vs pt = treatment arm proportion)
    szerolength<-length(szero); sonelength<-length(sone); m<-szerolength+sonelength;
    pc<-szerolength/m; pt<-1-pc;
    
    # Matrix of dimension m x 500 to generate perturbed values of TC and CT
    m1<-matrix(rexp(500*m,rate=1),ncol=500)
    
    # Bandwidth for kernel density estimation of TC
    bwd<-1.06 * sonelength^(-0.45) * 
      min(sqrt(var(sone)),(quantile(sone,0.75)-quantile(sone,0.25))/1.34)
    
    # Bandwidth for kernel density estimation of CT
    bwdct<-1.06 * szerolength^(-0.45) * 
      min(sqrt(var(szero)),(quantile(szero,0.75)-quantile(szero,0.25))/1.34)
    
    # Identify index numbers of yone/yzero vectors that have value 1
    yoneind <- which(yone==1)
    yzeroind <- which(yzero==1)
    
    ## Estimates
    
    # Empirical TT
    ttemp <- mean(yone)
    
    # Empirical CC
    ccemp <- mean(yzero)
    
    # tcpre = KDE computation for TC
    tcpre <- sapply(szero,function(s){sum(dnorm((sone[yoneind]-s)/bwd))/sum(dnorm((sone-s)/bwd))})
    # Identify index numbers of NA values
    extraTC<-which(is.nan(tcpre))
    # Update tcpre to account for these NA values by using the mpfr function to precisely estimate
    tcpre[extraTC]<-sapply(szero[extraTC],function(s){as.numeric(sum(dnorm(mpfr((sone[yoneind]-s)/bwd,53)))/sum(dnorm(mpfr((sone-s)/bwd,53))))})
    
    # Empirical TC (tcemp) 
    tcemp <- mean(as.numeric(tcpre))
    
    # ctpre = KDE computation for CT
    ctpre<-sapply(sone,function(s){sum(dnorm((s-szero[yzeroind])/bwdct))/sum(dnorm((s-szero)/bwdct))})
    # Identify index numbers of NA values
    extraCT<-which(is.nan(ctpre))
    # Update ctpre to account for these NA values by using the mpfr function to precisely estimate
    ctpre[extraCT]<-sapply(sone[extraCT],function(s){as.numeric(sum(dnorm(mpfr((s-szero[yzeroind])/bwdct,53)))/sum(dnorm(mpfr((s-szero)/bwdct,53))))})
    
    # Empirical CT (ctemp)
    ctemp <- mean(as.numeric(ctpre))
    
    # Empirical RS
    Rsemp <- (ttemp-tcemp)/(ttemp-ccemp)
    
    # Empirical CPM
    CPMemp <- pt*(ttemp-tcemp)+pc*(ctemp-ccemp)
    
    # Empirical CPI
    CPIemp <- if(pc==pt){0}else{(pc-pt)*(ttemp-tcemp-ctemp+ccemp)}
    
    # Empirical NPS
    NPSemp <- pt*(ttemp-ctemp)+pc*(tcemp-ccemp)
    
    # Empirical PMS
    PMSemp <- CPMemp^2/(CPMemp^2+CPIemp^2+NPSemp^2)
    
    ## Variances
    
    # Perturbed TT (vector of 500 values)
    TTperturb<-colSums(m1[1:sonelength,]*yone)/colSums(m1[1:sonelength,])
    
    # Perturbed CC (vector of 500 values)
    CCperturb<-colSums(m1[(1+sonelength):m,]*yzero)/colSums(m1[(1+sonelength):m,])
    
    # tcmat1 = component of numerator of KDE computation for perturbed TC (prior to using m1 matrix to perturb)
    tcmat1<-matrix(rep(sone[yoneind]/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sum(yone)),nrow=szerolength,byrow=FALSE)
    # tcmat2 = component of denominator of KDE computation for perturbed TC (prior to using m1 matrix to perturb)
    tcmat2<-matrix(rep(sone/bwd,szerolength),nrow=szerolength,byrow=TRUE)-matrix(rep(szero/bwd,sonelength),nrow=szerolength,byrow=FALSE)
    # TCperturb3 = KDE computation for perturbed TC
    TCperturb3<-(t(m1[yoneind,])%*%t(dnorm(tcmat1)))/(t(m1[1:sonelength,])%*%t(dnorm(tcmat2)))
    # Identify index numbers of NA values
    extracolTC<-which(is.na(colSums(TCperturb3)))
    
    # Perturbed TC (vector of 500 values) - if(no NA values then perturbed TC = ...) else(perturbed TC = ...)
    if(sum(extracolTC)==0){
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }else{
      for(j in 1:length(extracolTC)){
        TCperturb3[,extracolTC[j]]<-t(t(rep(yone[topn(abs(sone-szero[extracolTC[j]]),1L,decreasing=F)[1L]],500)))
      }
      TCperturb<-rowSums(t(m1[(1+sonelength):m,])*TCperturb3)/colSums(m1[(1+sonelength):m,])
    }
    
    # ctmat1 = component of numerator of KDE computation for perturbed CT (prior to using m1 matrix to perturb)
    ctmat1<-matrix(rep(szero[yzeroind]/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,sum(yzero)),nrow=sonelength,byrow=FALSE)
    # ctmat2 = component of denominator of KDE computation for perturbed CT (prior to using m1 matrix to perturb)
    ctmat2<-matrix(rep(szero/bwdct,sonelength),nrow=sonelength,byrow=TRUE)-matrix(rep(sone/bwdct,szerolength),nrow=sonelength,byrow=FALSE)
    # CTperturb3 = KDE computation for perturbed CT
    if(length(yzeroind)==1){
      CTperturb3<-(t(t(m1[yzeroind,]))%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }else{
      CTperturb3<-(t(m1[yzeroind,])%*%t(dnorm(ctmat1)))/(t(m1[1:szerolength,])%*%t(dnorm(ctmat2)))
    }
    extracolCT<-which(is.na(colSums(CTperturb3)))
    
    # Perturbed CT (vector of 500 values) - if(no NA values then...) else(...)
    if(sum(extracolCT)==0){
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }else{
      for(K in 1:length(extracolCT)){
        CTperturb3[,extracolCT[K]]<-t(t(rep(yzero[topn(abs(szero-sone[extracolCT[K]]),1L,decreasing=F)[1L]],500)))
      }
      CTperturb<-rowSums(t(m1[(1+szerolength):m,])*CTperturb3)/colSums(m1[(1+szerolength):m,])
    }
    
    # Perturbed RS (vector of 500 values)
    Rperturb<-as.numeric((TTperturb-TCperturb)/(TTperturb-CCperturb))
    
    # Perturbed CPM (vector of 500 values)
    CPMperturb<-as.numeric(pt*(TTperturb-TCperturb)+pc*(CTperturb-CCperturb))
    
    # Perturbed CPI (vector of 500 values)
    CPIperturb<-if(pc==pt){0}else{as.numeric((pc-pt)*(TTperturb-TCperturb-CTperturb+CCperturb))}
    
    # Perturbed NPS (vector of 500 values)
    NPSperturb<-as.numeric(pt*(TTperturb-CTperturb)+pc*(TCperturb-CCperturb))
    
    # Perturbed PMS (vector of 500 values)
    PMSperturb<-CPMperturb^2/(CPMperturb^2+CPIperturb^2+NPSperturb^2)
    
    # Generate the vector with the following 12 elements:
    # [1] empirical RS,[2] variance of 500 perturbed RS values, 
    # [3] lower 2.5% quantile of perturbed RS, [4] upper 2.5% quantile of perturbed RS,
    # [5] empirical PMS, [6] variance of 500 perturbed PMS values,
    # [7] lower 2.5% quantile of perturbed PMS, [8] upper 2.5% quantile of perturbed PMS,
    # [9] proportion of estimated TC values requiring mpfr function due to NA,
    # [10] proportion of estimated CT values requiring mpfr function due to NA,
    # [11] number of perturbed TC values requiring mpfr function due to NA,
    # [12] number of perturbed CT values requiring mpfr function due to NA
    return(c(Rsemp,var(Rperturb),
             quantile(Rperturb,0.025,names=F),quantile(Rperturb,0.975,names=F),
             PMSemp,var(PMSperturb),
             quantile(PMSperturb,0.025,names=F),quantile(PMSperturb,0.975,names=F),
             length(extraTC)/szerolength,length(extraCT)/sonelength,length(extracolTC),length(extracolCT)))
  }
  resul <- list()
  resul$seed<-i
  
  # Fixed parameters for NormalSubset200LL1
  muST<-5; sigmaST<-2; muSC<-4; sigmaSC<-1;
  aYT<-1;bYT<-2.5;aYC<-1; bYC<-4.9
  
  # 200 values of true arm of surrogate endpoint with mean muST=5 and standard deviation sigmaST=2
  # 200 values of true arm of surrogate endpoint with mean muST=4 and standard deviation sigmaST=1
  # Note allocation ratio (r) = 1
  ts1<-rnorm(200,mean=muST,sd=sigmaST);ts0<-rnorm(200,mean=muSC,sd=sigmaSC);
  
  # Logit link function for control arm of true endpoint
  functnorm0<-function(s){rbinom(1,1,1/(1+exp(-aYC*(s-bYC))))};
  # Logit link function for treatment arm of true endpoint
  functnorm1<-function(s){rbinom(1,1,1/(1+exp(-aYT*(s-bYT))))};
  # Apply logit link based on surrogate endpoint dataset
  ty1<-sapply(ts1,functnorm1);ty0<-sapply(ts0,functnorm0);
  
  # Apply WPLsimfast function to this particular dataset
  singlesim <- WPLsimfast(sone=ts1,szero=ts0,yone=ty1,yzero=ty0);
  resul$resulX<-data.frame(rep=i,Rsemp=singlesim[1],varRsemp=singlesim[2],lowerRs=singlesim[3],
                           upperRs=singlesim[4],PMSemp=singlesim[5],varPMSemp=singlesim[6],
                           lowerPMS=singlesim[7],upperPMS=singlesim[8],TCprop=singlesim[9],
                           CTprop=singlesim[10],TCSEprop=singlesim[11],CTSEprop=singlesim[12],row.names=NULL)
  colnames(resul$resulX) <- c("rep","Rsemp","varRsemp","lowerRs","upperRs","PMSemp","varPMSemp","lowerPMS","upperPMS",
                              "TCprop","CTprop","TCSEprop","CTSEprop")
  resul
}

# True distributions of surrogate and true endpoints under control and treatment arms based on fixed parameters for NormalSubset200LL1
muST<-5; sigmaST<-2; muSC<-4; sigmaSC<-1;aYT<-1;bYT<-2.5;aYC<-1; bYC<-4.9
YTYT<-function(x){1/(1+exp(-aYT*(x-bYT)))};YCYC<-function(x){1/(1+exp(-aYC*(x-bYC)))}
STST<-function(x){dnorm(x,mean=muST,sd=sigmaST)};SCSC<-function(x){dnorm(x,mean=muSC,sd=sigmaSC)}

# True value of TC
TCTC<-integrate(function(x){YTYT(x)*SCSC(x)},-Inf,Inf)

# True value of TT
TTTT<-integrate(function(x){YTYT(x)*STST(x)},-Inf,Inf)

# True value of CC
CCCC<-integrate(function(x){YCYC(x)*SCSC(x)},-Inf,Inf)

# True value of CT
CTCT<-integrate(function(x){YCYC(x)*STST(x)},-Inf,Inf)

# True value of RS
RSRS<-(TTTT$value-TCTC$value)/(TTTT$value-CCCC$value)

# True value of RS'
RSPRSP<-(CTCT$value-CCCC$value)/(TTTT$value-CCCC$value)

# True value of PMS as a function of r (allocation ratio)
PMSPMS<-function(r){(r*RSRS+RSPRSP)^2/((r*RSRS+RSPRSP)^2+(r-1)^2*(RSRS-RSPRSP)^2+(1-RSRS+r*(1-RSPRSP))^2)}

# Generate table of performance measures in this case - equal allocation (r=1) to compute PMS
finalperfdf<-finalperf(resdf,RSRS,PMSPMS(1))
NormalSubset200LL1<-finalperfdf[[1]]
end.time<-Sys.time()

# Time how long simulation ran
end.time-start.time

# Final output of table of performance measures (bias, empirical SE, MSE, average model SE, coverage) for each of RS and PMS with associated Monte Carlo SEs, and median extrapolation proportions of TC and CT (not their perturbed versions) with interquartile ranges (IQRs) as final column
NormalSubset200LL1