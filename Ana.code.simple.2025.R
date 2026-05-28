# environment ====

# install.packages(pkgs=c("zoo","tidyverse","lubridate"))
library(zoo)
library(tidyverse)
library(lubridate)

# data in


#==== FUNCTIONS  all start with "z" ====
source("AnaFunctions.R")

#==== Generate dataframe ====
source("GenerateData.R")

head(df)
head(sweApr1)
head(metricz)
#=== Import dataframe ====
# ANA.DF2 <- read.csv("./data/anatone.dataframe.final.csv")



#==== WAS Fit 5 param Air regime ====

airP5.init <- c(log(2), 11.64,11.04,164.3,20,0)  # initial air regime parameters, ADDED log(sigma)

I <- air2$year >= 1988
Airfit5 <- optim(par=airP5.init,fn=zSineAir5fit,df=air2[I,],method="BFGS",control=list(trace=3), hessian=TRUE)
# This is the air regime: Modeled air temperature
AirModelpar1 <- Airfit5$par
names(AirModelpar1) <- c("log.sigma", "M","N","P","J","K")
AirModelpar <- Airfit5$par[-1]
names(AirModelpar) <- c("M","N","P","J","K")

# CI for sigma included: 
Airfit5.cov = solve(Airfit5$hessian)  #returns the inverse
Airmodelpar1.SE = sqrt(diag(Airfit5.cov))  #SE's of estimated parameters
Airfit5.summary = data.frame(AirModelpar1, zgetFitSummary(Airfit5))
# For returning CI of sigma after fitting sigma on the log scale
# The easy way is to do the following:
Airfit5.sigma = exp(Airfit5.summary$ParmEst[1])
Airfit5.sigma.CI = exp(c(Airfit5.summary$CI.low[1], Airfit5.summary$CI.up[1]))


Mair <- Mregimeair <- cbind.data.frame("doy"=1:365,AirRegime=zSineAir5(AirModelpar,1:365))

#=== Assemble df ====
usedf <- left_join(df,metricz,by="year")
usedf <- left_join(usedf,Mair,by="doy")
usedf <- usedf[usedf$doy <= 365,]

# usedf is a long timeseries. Not all years have SWE for fitting 
print(head(usedf))
print(length(unique(usedf$year)))
print(length(unique(usedf$year[!is.na(usedf$swe)])))

ANA.DF <- usedf[!is.na(usedf$swe) & usedf$doy != 366,]
print(t(head(ANA.DF)))

#==== Fit 6 param Water regime ====

AnaModel6 <- optim(par=c(14,12,140,9,200,300), fn=zHiatusFit, 
                  data=ANA.DF$ana.temp.fill,
                  method="BFGS",
                  control=list(maxit=100000,trace=2), hessian=TRUE) 

AnaModel6par <- AnaModel6$par
names(AnaModel6par) <- c("AA","BB","CC","DD","EE","FF")
Regime <- zSinFuncHiatusParVec(1:365,AnaModel6par)

# diagnostics
AnaModel6.cov = solve(AnaModel6$hessian)  #returns the inverse
AnaModel6par1.SE = sqrt(diag(AnaModel6.cov))  #SE's of estimated parameters
AnaModel6.summary = data.frame(AnaModel6par, zgetFitSummary(AnaModel6))
# For returning CI of sigma after fitting sigma on the log scale
# The easy way is to do the following:
AnaModel6.sigma = exp(AnaModel6.summary$ParmEst[1])
AnaModel6.sigma.CI = exp(c(AnaModel6.summary$CI.low[1], AnaModel6.summary$CI.up[1]))

ANA.DF <- left_join(ANA.DF,cbind.data.frame("doy"=1:365,"Regime"=Regime),by=c("doy"))

#==== Fit annual Water regime  ====
Eachyear6water <- NULL
anawateryears <- unique(df$year)
for(y in anawateryears){
  I <- df$year==y
  x <- df[I,]
  fit = optim(par=c(12,10,160,4,170,311), fn=zHiatusFit, data=x$ana.temp.fill,method="L-BFGS-B",
              lower=c(8,6,140,0,125,280), upper=c(16,15,175,10,250,350),
              control=list(maxit=100000,trace=0)) 
  x2x <- zSinFuncHiatusParVec(1:365,fit$par)
  Eachyear6water <- rbind.data.frame(Eachyear6water,cbind.data.frame(year=y,t(fit$par),maxtemp=max(x2x),maxtempday=x$doy[x2x == max(x2x)]) )
}
names(Eachyear6water)[1:7] <- c("year","A","B","C","D","E","F")


globalfixedpars <- c(AnaModel6par['CC'],
                     AnaModel6par['FF'],
                     AirModelpar['P'],
                     AirModelpar['J'],
                     AirModelpar['K'],
                     AnaModel6par['EE'])


#==== Fit annual air regime ====
Eachyear5air <- NULL
airyears <- unique(df$year)
airyears <- airyears[airyears <= 2022]
for(y in airyears){
  I <- df$year==y
  x <- df[I & df$doy %in% 1:365,]
  # Remember to chope off the log(sigma) in the front
  fit <- optim(par=airP5.init[-1],fn=zSineAir5fit.v0,df=x,method="BFGS",control=list(trace=3))
  
  x2x <- zSineAir5(fit$par,days=1:365)
  Eachyear5air <- rbind.data.frame(Eachyear5air,cbind.data.frame(year=y,t(fit$par),maxtemp=max(x2x),maxtempday=x$doy[x2x == max(x2x)]) )
}
rm(x);rm(x2x)

names(Eachyear5air)[1:6] <- c("year","M","N","P","J","K")


#==== Fit Air regime with covariates ====
#

AnaAir <- optim(par=c(log(2),1,1,1,1),fn=zAnaFormAirfit.wsigma, mydf=ANA.DF,
                fixedpars=globalfixedpars,
                method="BFGS",control=list(maxit=1000,trace=3),hessian=TRUE)
  AnaAirpar <- AnaAir$par; 
  names(AnaAirpar) <- c("log.sigma","m0","m1","n0","n1")
  
  # Diagnostics
  ## JF Example usage. CI for sigma included: 
  AnaAir.cov = solve(AnaAir$hessian)  #returns the inverse
  AnaAirpar1.SE = sqrt(diag(AnaAir.cov))  #SE's of estimated parameters
  AnaAir.summary = data.frame(AnaAirpar, zgetFitSummary(AnaAir))
  # For returning CI of sigma after fitting sigma on the log scale
  # The easy way is to do the following:
  AnaAir.sigma = exp(AnaAir.summary$ParmEst[1])
  AnaAir.sigma.CI = exp(c(AnaAir.summary$CI.low[1], AnaAir.summary$CI.up[1]))
  
  
# Add the air parameters to the fixed set.
globalfixedpars <- c(globalfixedpars,AnaAirpar[-1])




#==== Fit Daily Water model with 6 parameters AKA type=2 ==== 

alignside <- "right" # needed in code 
waterpar.init    <-    c(1, 1,  4, .5, 170,  -.3,  .5,  .3)

# AnaFit2 does not try to fit E so the globalfixedpars value is used.

AnaWater3 <- optim(par=c(log(2),waterpar.init[c(1:4,7,8)]), fn=zAnaFit2.wsigma, daf=ANA.DF, fixedpars=globalfixedpars,useairlag=3,
                   method="BFGS",control=list(maxit=1000,trace=2),hessian=TRUE)

AllPar3 <- AnaWater3$par[-1]
names(AllPar3) <- c("a0","a1","b0","b1","d","p0")
AllPar3.1 <- AnaWater3$par
names(AllPar3.1) <- c("log.sigma","a0","a1","b0","b1","d","p0")

AnaWater <- AnaWater3 ;useairlag <- 3 ; AllPar <- AllPar3
print(AllPar)

# Diagnostics
## JF Example usage. CI for sigma included: 
AnaWater.cov = solve(AnaWater$hessian)  #returns the inverse
AnaWaterpar1.SE = sqrt(diag(AnaWater.cov))  #SE's of estimated parameters
AnaWater.summary = data.frame(AllPar3.1, zgetFitSummary(AnaWater))
# For returning CI of sigma after fitting sigma on the log scale
# The easy way is to do the following:
AnaWater.sigma = exp(AnaWater.summary$ParmEst[1])
AnaWater.sigma.CI = exp(c(AnaWater.summary$CI.low[1], AnaWater.summary$CI.up[1]))



# For printing to manuscript table of coeffiients etc.:
print(Airfit5.summary)
print(AnaModel6.summary)
print(AnaWater.summary)
print(AnaAir.summary)


#==== Leave-one-out evaluation ====
# useairlag = 3 is preferred 

if(1){  ##  slow. Do it only when the time for the 34 calibrations is no matter
  # set up for eiher the 8 parameter or 6 parameter version
  LOYOpars <- NULL
  years <- unique(ANA.DF$year)
  for(y in years){
    AnaModeloneout <- optim(par=waterpar.init[c(1:4,7,8)],fn=zAnaFit2,daf=ANA.DF[ANA.DF$year != y,],fixedpars=globalfixedpars,useairlag=3,
                            method="BFGS",control=list(maxit=10000,trace=2))
    LOYOpars <- rbind.data.frame(LOYOpars,c(y,AnaModeloneout$par))
  } ; cat("\n");
  # names(LOYOpars ) <- c("year","a0","a1","b0","b1","e0","e1","d","p0")
  names(LOYOpars ) <- c("year","a0","a1","b0","b1","d","p0")
} # end leaveoneout run

#==== leave one out analysis ====
if(1){ # leaveoneout analysis
  leaveoutresults <- NULL
  mod.obs.diffs.leaveone <- NULL
  mod.one.out.vals <- NULL
  useairlag <- 3

  for(y in years){
    legwords <- leglines <- legcols <- NULL
    AnaModeloneout <- LOYOpars[LOYOpars$year == y,][-1]
    
    K <- ANA.DF$year==y 
    # 8 param w/ zAnaFit
    # MD <- zAnaFit(AnaModeloneout,ANA.DF[K,],fixedpars=globalfixedpars,retype="pred",useairlag=useairlag)
    # 6 param w/ zAnaFit2
    MD <- zAnaFit2(AnaModeloneout,ANA.DF[K,],fixedpars=globalfixedpars,retype="pred",useairlag=useairlag)
    # Regime <- zAnaForm2(AnaModeloneout,ANA.DF[K,],fixedpars=globalfixedpars)
    airOb <- ANA.DF$Tair[K]
    Ob <- ANA.DF$ana.temp.fill[K]
    x <- ANA.DF$doy[K]
    flow <- ANA.DF$ana.flow[K]
    airRegime <- zAnaFormAir(globalfixedpars[c("m0","m1","n0","n1")],ANA.DF[K,],fixedpars=globalfixedpars)
    Q <- ANA.DF$ana.flow[K]
    M1 <- zSinFuncHiatusParVec(1:365,Eachyear6water[Eachyear6water$year==y,2:7])
    air5fit <- zSineAir5(as.numeric(Eachyear5air[Eachyear5air$year==y,2:6]),days=1:365) 
    airR <- airOb - air5fit  # the one day residual
    air5R <- rollmean(airOb - air5fit,k=useairlag,align=alignside,fill=NA) # the smoothed residual using the useairlag
    Q <- ANA.DF$ana.flow[K]
    
    mad <- mean(abs(Ob-MD),na.rm=T)
    mre <- mean(Ob-MD,na.rm=TRUE)
    rmse <- sqrt(mean((MD-Ob)^2))
    diffs <- Ob - MD 
    cumerror <- sum(Ob - MD)
    cumsumerror <- cumsum(Ob-MD)
    worstcumsum <- max(abs(cumsumerror))
    R2 <- summary(lm(Ob ~ MD))$r.squared
    # meandiff7 <- mean(zrunavg(Ob-MD,7))
    meandiff7 <- mean(rollmean(Ob-MD,7,align=alignside,fill=NA),na.rm=T)
    
    # diff7days <- quantile((zrunavg(MD - Ob,7)),probs=c(0.5,0.75,0.80,0.825,0.85,0.875,0.90,1))
    diff7days <- quantile((rollmean(MD - Ob,k=7,align=alignside,fill=NA)),probs=c(0.5,0.75,0.80,0.825,0.85,0.875,0.90,1),na.rm=TRUE)
    summerMAD <- mean(abs(diffs[274:335]))
    winterMAD <- mean(abs(diffs[193:151]))
    leaveoutresults <- rbind.data.frame(leaveoutresults,
                                        cbind.data.frame(leaveoutyear=y,mad=mad,mre=mre,rmse=rmse,cumerror=cumerror,varExp=R2,
                                                         worstcumsum=worstcumsum,error7max=diff7days[["100%"]],
                                                         error7.90=diff7days[["90%"]],error7.80=diff7days[["80%"]],
                                                         error7.50=diff7days[["50%"]],summerMAD,winterMAD,meandiff7))
    mod.obs.diffs.leaveone <- c(mod.obs.diffs.leaveone,diffs)
    mod.one.out.vals <- c(mod.one.out.vals,MD)
    
  } # end for(year) loop
} # end leaveoneout analysis
print(leaveoutresults)
for(i in 2:14) print(sprintf("%11s %5s %3.4f %5s %3.4f",names(leaveoutresults[i])," min:",min(leaveoutresults[,i])," max:",max(leaveoutresults[,i])))
OB <- ANA.DF$ana.temp.fill
probss <- c(0,0.1,0.125,0.15,0.2,0.25,0.5,0.75,0.80,0.825,0.85,0.875,0.90,1)
diff7days <- quantile((rollmean(OB-mod.one.out.vals,7,align=alignside,fill=NA)),probs=probss,na.rm=TRUE)
mediandiff7 <- median(rollmean(OB-mod.one.out.vals,7,align=alignside,fill=NA),na.rm=TRUE)
meandiff7 <- mean(rollmean(OB-mod.one.out.vals,7,align=alignside,fill=NA),na.rm=TRUE)
diff10days <- quantile(abs(rollmean(OB-mod.one.out.vals,10,align=alignside,fill=NA)),probs=probss,na.rm=TRUE)

cat(paste(
  "   Mean7day = ", meandiff7,"\n",
  "Median7day = ", mediandiff7,"\n",
  "   Max7day = ", diff7days[["100%"]],"\n",
  "       MAD = ",  mean(abs(OB-mod.one.out.vals),na.rm=TRUE),"\n",
  "      RMSE = ", sqrt(mean((OB-mod.one.out.vals)^2,na.rm=TRUE)),"\n",
  "       MRE = ", mean(OB-mod.one.out.vals,na.rm=TRUE),"\n",
  "        R2 = ", summary(lm(OB ~ mod.one.out.vals))$r.squared,"\n"
))

print(summary(leaveoutresults))

#==== Put the annual regimes into the DF ====
ANA.DF2 <- cbind.data.frame(ANA.DF,annualRegime=NA,annualFit=NA,dailyFittype2 = NA,dailyFit=NA,dailyFittype2.a1=NA )

for(y in unique(ANA.DF2$year)){
  # type 1 is old method wit regression for 'E' is depricated
    K <- ANA.DF2$year==y 
    Regime <- zAnaForm2(AllPar3,ANA.DF2[K,],fixedpars=globalfixedpars)
    ANA.DF2$annualRegime[K] <- Regime
    ANA.DF2$annualFit[K] <- MD
    ANA.DF2$dailyFittype2[K] <- zAnaFit2(AllPar3,ANA.DF[K,],fixedpars=globalfixedpars,retype="pred",useairlag=useairlag)
}

#==== Generate results and write out files. ====

metaresults <- z5(pdff="Results.type2.pdf",type=2)

z5(years=2017,type=2)


# Summary ?
dim(ANA.DF2)   #  [1] 12775    28
print(names(ANA.DF2))
zf(2,2)
par(mar=c(4,4,1,1))
zplotfitDF(ANA.DF2,"annualRegime","ana.temp")
zplotfitDF(ANA.DF2,"annualFit","ana.temp")
zplotfitDF(ANA.DF2,"dailyFittype2","ana.temp")


#==== References to key R objects
print(AnaModel6par)
print(AirModelpar)
print(AllPar3)


zf(2,3)
zplotfitDF(df=Eachyear6water,"year","A",pch=16,showANAlines=TRUE)
zplotfitDF(df=Eachyear6water,"year","B",pch=16,showANAlines=TRUE)
zplotfitDF(df=Eachyear6water,"year","E",pch=16,showANAlines=TRUE)
zplotfitDF(df=Eachyear5air,"year","M",pch=16,showlines=TRUE)
zplotfitDF(df=metrics,"year","winterairDJ",pch=16,showlines=TRUE,ylimm=c(-4.5,5.2))
zplotfitDF(df=metricz,"year","swe",pch=16,showlines=TRUE)

#
#==== PLOT: COMBINED air and water ====
zf(1,1)
airmod <- zSineAir5(AirModelpar,1:365)
dayairmeans <- usedf %>% group_by(doy) %>% summarise(avg=mean(Tair,na.rm=T),median=median(Tair,na.rm=T))
regimeairdayerrormeans <- airmod - dayairmeans$avg
regimeairdayerrormedians <- airmod - dayairmeans$median
regimedayRMSE <- sqrt(mean(regimeairdayerrormeans^2))

for(type in c("png","","pdf")){    
# for(type in c("")){      
  if(type =="pdf")  pdf(file="allregime.pdf",width=11,height=13)
  if(type =="png")  png(file="allregime.png",width=800,height=1000)
  {
  par(mar=c(3,5,1,5),cex.lab=1.2,cex.axis=1.1)
  plot(1:365,1:365,ylim=c(-25,32),xlab="",ylab="",type="n",axes=F);box()
  par(mgp=c(3,0.5,0))
  zxWYdates() ; par(mgp=c(3,1,0))
  par(cex.lab=1.5,cex.axis=1.3)
  x <- c(daylength$time[274:365],daylength$time[1:273])
  lines(1:365,(x-min(x))*1 + min(x)-35,type="l",col="grey40",lwd=8)
  lines(1:365,(x-min(x))*1 + min(x)-35,type="l",col="white",lwd=1)
  # All ref lines:
  zvabline(v=c(82,264,301,314),f=c(0,0.12,0.85,0.49),coll="grey40",ltyy=3)
  text(c(82,264,301,314)-7, -24.5,c("","  Insolation Peak","Air Peak","Water Peak"),srt=90,las=0)
  
  legend("bottomleft",inset=c(0.15,0.06),bty="n",lty=1,col= "grey40",lwd=8, legend=c("Solar\nInsolation"),cex=1.5)
  legend("bottomleft",inset=c(0.15,0.06),bty="n",lty=1,col= "white",lwd=1, legend=c("Solar\nInsolation"),cex=1.5)
   y <- c(8,9,10,11,12,13,14,15,16)
  yat <- (y-min(x))*1 + min(x)-35
  axis(2,at = yat,labels=y,col.axis="grey40")
  mtext("Daylength (hours)",side=2,line=2.5,adj=0,las=3,cex=1.8,col="grey40")

  palegreen <- rgb(10,20,10,5,NULL,25)
  palegreen2 <- rgb(10,20,10,10,NULL,25)
  paleblue <- rgb(10,10,25,5,NULL,25)
  paleblue2 <- rgb(10,10,25,10,NULL,25)
  
  points(ANA.DF$doy,ANA.DF$Tair,col=palegreen,pch=16,cex=0.7);
  medz2 <- rep(NA,365)
  for(i in 1:365)medz2[i] <- median(ANA.DF$Tair[ANA.DF$doy==i],na.rm=T)
  mtext(expression(paste("Air temperature ",degree,"C")),2,line=2,col="darkgreen",cex=2,las=3)
  par(las=1)
  axis(2,at=c(-10,-5,0,5,10,15,20,25,30),col.axis="darkgreen")
  axis(2,at=seq(-10,30,by=1),labels=NA)
 

  maxday <- c(1:length(airmod))[airmod==max(airmod)]
  lines(1:365,dayairmeans$avg,col="white",lwd=7)
  lines(1:365,dayairmeans$avg,col="brown",lwd=5)
  # lines(1:365,dayairmeans$median,col="white",lwd=7)
  # lines(1:365,dayairmeans$median,col="orange",lwd=5)
  lines(1:365,airmod,lwd=9,col="white")
  lines(1:365,airmod,lwd=7,col="darkgreen")
  AF5p <- AirModelpar
  lines(1:365,zSineAir5(c(AF5p[1],  AF5p[2], AF5p[3], 0 , 0 )),col="darkgreen",lwd=2,lty=1)
 
  legend("topleft",bty="n",inset=c(0.0,0),legend=c("Air Thermal Regime and data\nLewiston, ID (1962-2022)"),adj=0,cex=1.4,text.col="darkgreen")
  legend("topleft",bty="n",inset=c(0.04,0.08),
         legend=c("Multi-year Air regime","Air regime w/o asymmetry","Day of year mean","Data Cloud 1962-2022"),
         lwd=c(7,2,5,10),col=c("darkgreen","darkgreen","brown",palegreen2),lty=c(1,1,1,1),cex=1)
  
  cexx <- 1.1
  print(AF5p)
  legend("topleft",inset=c(0.12,0.2),bty="n",legend=c(
    paste("M =",round(AF5p[1],2)),
    paste("N =",round(AF5p[2],2)),
    paste("P =",round(AF5p[3],1)),
    paste("J =",round(AF5p[4],1)),
    paste("K =",round(AF5p[5],1))
  ),cex=cexx,text.col="darkgreen")
  
  
  segments(-2,AF5p[1],70,AF5p[1],lwd=2)
  rect(70,12,110,14,col="white",border=NA)
  text(90,AF5p[1],"Mean: M",cex=cexx,col="darkgreen")
  arrows(1,(AF5p[1]-AF5p[2]),1,(AF5p[1]+AF5p[2]),code=3,angle=20,length=0.1,lwd=2)
  text(1,0,expression(Range: M %+-% N),adj=0.1,cex=cexx,col="darkgreen")
  arrows(365-AF5p[3],27.5,365,27.5,code=3,angle=20,length=0.1,lwd=2)
  rect(233,26.7,273,28.2,col="white",border=NA)
  text(253,27.5," Phase: P",cex=cexx,col="darkgreen")
  
  
  segments(365-AF5p[5],7,365-AF5p[5],zSineAir5(c(AF5p),365-AF5p[5]),lwd=2)
  segments(365,7,365,29,lwd=2)
  arrows(365-AF5p[5],10,365,10,code=3,angle=20,length=0.1,lwd=2)
  rect(272,9,320,12,col="white",border=NA)
  text(295,10.25,adj=0,paste("Asymmetry (K)\noffset to P"),cex=cexx,col="darkgreen")
  
  arrows(365-AF5p[5]-91.5,zSineAir5(c(AF5p),365-AF5p[5]-91.5),365-AF5p[5]-91.5-AF5p[4],zSineAir5(c(AF5p),365-AF5p[5]-91.5),code=3,angle=20,length=0.1,lwd=4,col="white")
  arrows(365-AF5p[5]-91.5,zSineAir5(c(AF5p),365-AF5p[5]-91.5),365-AF5p[5]-91.5-AF5p[4],zSineAir5(c(AF5p),365-AF5p[5]-91.5),code=3,angle=20,length=0.1,lwd=2)
  # rect(200,21.2,262,22.8,col="white",border=NA)
  text(190,21.5,paste("Asymmetry: J"),cex=cexx,col="darkgreen")
  print(paste("RMSE:",regimedayRMSE))
  }

  {
mod <- zSinFuncHiatusParVec(paramvals = AnaModel6par)
# usedf is 52 years. ANA.DF only encompasses the SNOTEL/SWE available years
daymeans <- usedf %>% group_by(doy) %>% summarise(avg=mean(ana.temp.fill),median=median(ana.temp.fill))
regimedayerrormeans <- mod - daymeans$avg
regimedayerrormedians <- mod - daymeans$median
regimedayRMSE <- sqrt(mean(regimedayerrormeans^2))
regimeQs <- NULL

{

axis(4,at=c(-20,-15,-10,-5, 0),labels=c(0,5,10,15,20),col.axis="blue")
axis(4,at=seq(-20,0,by=1),labels=NA)
mtext(expression(paste("Water temperature ",degree,"C")),4,line=3,adj=0.2,cex=2,col="blue",las=3)
points(ANA.DF$doy,ANA.DF$ana.temp-20,col=paleblue,pch=16,cex=0.7)

lines(1:365,daymeans$avg-20,col="brown",lwd=3)
  # lines(1:365,mod,col="white",lwd=5)
  lines(1:365,mod-20,col="darkblue",lwd=5)
  lines(1:365,zSine(AnaModel6par)-20,col="darkblue",lwd=2,lty=1)

  legend("topright",bty="n",inset=c(0.0,0.0),legend=c("Water Thermal Regime and data\nAnatone, WA gauge (1962-2022)"),adj=0,cex=1.4,text.col="darkblue")
  rect(270,-17,330,-12,col="white",border=NA)
  legend("bottomright",border="white", bty="n", inset=c(0.10,0.16),
         legend=c("Multi-year Water Thermal Regime","Regime w/o snowmelt effect","Day of year mean ","Data Cloud 1962-2022"),
         lwd=c(7,2,5,10),col=c("darkblue","darkblue","brown",paleblue2),
         cex=1.1)
  cexx <- 1.1
  legend("bottomleft",inset=c(0.28,0.25),bty="n",legend=c(
    paste("A =",round(AnaModel6par[1],2)),
    paste("B =",round(AnaModel6par[2],2)),
    paste("C =",round(AnaModel6par[3],0)),
    paste("D =",round(AnaModel6par[4],2)),
    paste("E =",round(AnaModel6par[5],0)),
    paste("F =",round(AnaModel6par[6],0))
  ),cex=cexx,text.col="darkblue")
  print(paste("RMSE",sqrt(mean(regimedayerrormeans^2))))
  # arrows(1,AnaModel6par[1],300,AnaModel6par[1],code=3,angle=20,length=0.1,lwd=2)
  segments(-2,AnaModel6par[1]-20,70,AnaModel6par[1]-20,lwd=2)
  rect(38,12.3-20,70,13.7-20,col="white",border=NA)
  text(55,AnaModel6par[1]-20,"Mean: A",cex=cexx,col="darkblue")
  arrows(368,(AnaModel6par[1]-AnaModel6par[2])-20,368,(AnaModel6par[1]+AnaModel6par[2])-20,code=3,angle=20,length=0.1,lwd=2)
  text(370,0.5-20,expression(Range: A %+-% B),adj=1,cex=cexx,col="darkblue")
  arrows(365-AnaModel6par[3],25-20,365,25-20,code=3,angle=20,length=0.1,lwd=2)
  rect(280,24.5-20,320,25.5-20,col="white",border=NA)
  text(300,25-20,"Phase: C",cex=cexx,col="darkblue")
  junk <- (AnaModel6par[6]+AnaModel6par[5])/2
  junk2 <- (zSine(AnaModel6par,junk)+mod[junk])/2
  arrows(junk,mod[junk]-20,junk,zSine(AnaModel6par,junk)-20,code=3,angle=20,length=0.14,lwd=5,col="white")
  arrows(junk,mod[junk]-20,junk,zSine(AnaModel6par,junk)-20,code=3,angle=20,length=0.1,lwd=2)
  rect(205,junk2-0.5-20,240,junk2+0.5-20,col="white",border=NA)
  text(210,junk2-20,"Snow Effect: D",cex=cexx,col="darkblue")
  segments(AnaModel6par[5],2.5-20,AnaModel6par[5],mod[AnaModel6par[5]]-20,col="white",lwd=4)
  segments(AnaModel6par[5],2.5-20,AnaModel6par[5],mod[AnaModel6par[5]]-20,lwd=2)
  segments(AnaModel6par[6],16-20,AnaModel6par[6],mod[AnaModel6par[6]]-20,col="white",lwd=4)
  segments(AnaModel6par[6],16-20,AnaModel6par[6],mod[AnaModel6par[6]]-20,lwd=2)
  text(AnaModel6par[5]+10,1-20,"Begin Melt: E",cex=cexx,adj=0.35,col="darkblue")
  text(AnaModel6par[6],15.5-20,"End Melt: F",cex=cexx,col="darkblue",adj=0.5)
}
  

if(type =="pdf")  dev.off()
if(type =="png")  dev.off()
  } 
  }
# end plot fig type  pdf

# Summarize metrics  ====

for(i in 2:ncol(metricz) )print(paste(names(metricz)[i],signif(var(metricz[,i],na.rm=TRUE),3),signif(sqrt(var(metricz[,i],na.rm=TRUE)))))
