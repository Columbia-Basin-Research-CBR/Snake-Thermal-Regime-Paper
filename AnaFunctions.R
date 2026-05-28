# AnaFunctions.R

#==== FUNCTIONS  all start with "z" ====

zgetFitSummary = function(fitobj, alpha=0.05) {
  fxhess <- fitobj$hess
  fvcb <- solve(fxhess)  #vcov of param ests
  if ( (any(diag(fvcb)<0) | any(eigen(fvcb)$values<0)) ){
    cat("WARNING: Hessian failure!\n")
  }
  fpar <- fitobj$par
  fse <- round(sqrt(diag(fvcb)), 6)
  fz <- round(fpar/fse, 4)
  fpval <- round(2*pnorm(-abs(fz)), 6)
  fmult = qnorm(1-alpha/2)
  fcil = round(fpar - fmult*fse, 6)
  fciu = round(fpar + fmult*fse, 6)
  fstab <- 	cbind(ParmEst=fpar, SE=fse, Zstat=fz, Pval=fpval, CI.low=fcil, CI.up=fciu)
  return(fstab)
}

zvabline <- function(v,f,coll=1,lwdd=2,label=FALSE,yoff=0,ltyy=1){  # truncates the abline(v=...) values at frctions of tth plot
  for( i in 1:length(v)) {
    segments(v[i],par()$usr[3],v[i],f[i]*par()$usr[4] + par()$usr[3]*(1-f[i]),col=coll,lwd=lwdd,lty=ltyy)
    if(label){text(v[i],(0.985-yoff)*par()$usr[4] + par()$usr[3]*(1- (0.985-yoff)),v[i],cex=0.75,col=coll) }
  }}

zSine <- function(pv,days=1:365){
  pv[1] + pv[2]*sin(2*pi/365*(days + pv[3]))  
} 
zFiller <- function(x){ # linear interpolation of values to fill NA in a vector
  y <- x
  if(any(!is.na(x))){
    j <- length(x)
    ok <- 0
    if(is.na(x[1])){
      ok <- which(!is.na(x))[1]
      y[1:ok] <- rep(x[ok])
    }
    i <- ok+1
    while(i <= j){
      k <- 0
      if(is.na(x[i])){ # find next good value
        k <- 0
        while(is.na(x[i+k] && (i+k) < j)){k <- k+1}
        # print(paste(x[ok],x[i+k],i,k,j,ok))
        if((i+k) == j & is.na(x[i+k]) ){ # lat one is NA
          y[i:(i+k)] <- x[ok]
        } else  {
          # print(paste(i,j,k,ok))
          y[ok:(i+k)] <- seq(x[ok],x[i+k],length.out=length(ok:(i+k)))
        }
        ok <- i+k
        i <- i+k+1
      } else {
        y[i] <- x[i]
        ok <- i
        i <- i+1
      }
    }
    return(y)
  } else{return(NA)}
}
zDOY.WY <- function(x){ # if x has these columns and is sanitized so day 366 is ONLY in calendar years like (1992,1996,2000,2004,2008,2012,2016,2020,2024,2028)
  names(x)[1:2] <- c("year","DOY")
  x <- cbind.data.frame(x,"WY"=0,"WYdoy"=0)
  for( i in 1:dim(x)[1]){
    x$WYdoy[i] <- x$DOY[i] + 92  
    x$WY[i] <- x$year[i]
    test <- 273 ; what <- 365 # non-leap values
    if(match(x$year[i],seq(1940,2080,by=4),nomatch=0) > 0){
      test <- 274; what <- 366;
    }
    if(x$DOY[i] > test){
      x$WYdoy[i] <- x$WYdoy[i] - what ; 
      x$WY[i] <-x$WY[i] + 1 
    }
  }
  return(x)
}
zxWYdates <- function(with=T, axiss=1,labs=TRUE,line=3){
  at <- c(1,32,60,91,121,152,182,213,244,274,305,335)+92
  at[at > 365] <- at[at > 365] -365
  lab <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
  if(labs)axis(axiss,at=at,labels=lab)
  if(with & axiss == 1) abline(v=at,lwd=1,col="grey80",lty=3)
  if(with & axiss == 2) abline(h=at,lwd=1,col="grey80",lty=3)
}


zSineAir5fit.v0 <- function(pv,df,retype="LL"){
  obs <- df$Tair
  pred <-  zSineAir5(pv,days=df$doy)
  diffs = obs-pred
  likes = dnorm(diffs, 0, 1)
  if(retype=="LL") { return (-sum(log(likes))) }
  if(retype=="resid") { return(diffs) } else { return (pred) }
}

zSineAir5fit <- function(pv,df,retype="LL"){
  # NOTE: pv[1] is log(sigma), where sigma is residual standard deviation
  obs <- df$Tair 
  pred <-  zSineAir5(pv[-1],days=df$doy)
  diffs = obs-pred
  sigma = exp(pv[1])
  likes = dnorm(diffs, 0, sigma)
  if(retype=="LL") { return (-sum(log(likes))) } 
  if(retype=="resid") { return(diffs) } else { return (pred) }
}

zSineAir5 <- function(pv,days=1:365){
  pv[1] + pv[2]*sin(2*pi/365*(days + pv[3]+(0+pv[4]*sin(2*pi/365*(days+pv[5])))))  
} 

zHiatusFit = function(paramvals, data, retype="LL"){
  obs = data
  pred = zSinFuncHiatusParVec(day=1:length(data),paramvals[1:6])
  test <- 91.25 - paramvals[3]
  if(paramvals[2] < 0) test <- test+180
  diffs = (obs-pred)
  mu <- 0
  sigma <- 1
  likes <- dnorm(diffs, mu, sigma)
  if(retype=="LL") {
    out <- -sum(log(likes))
    return (out)
  } else {
    return (pred)
  }
}
zSinFuncHiatusParVec <- function(day=1:365, paramvals){
  # print(paramvals)
  aa <- as.numeric(paramvals[1])
  bb <- as.numeric(paramvals[2])
  cc <- as.numeric(paramvals[3])
  mag <- as.numeric(paramvals[4])
  begin <- as.numeric(paramvals[5])
  end <- as.numeric(paramvals[6])
  test <- day %% 366 >= begin & (day %% 366) <= end
  out <- aa + bb * sin((2 * pi)/365 * (day + cc))
  out2 <- out
  out2[test] <- out[test] - (mag/2 - mag/2 * sin((2 * pi)/(end - begin) * (day %% 365 - begin + (end - begin)/4)))[test]
  #  if(begin >= 0 & mag >= 0 & end > begin) {
  out[test] <- out2[test]
  # }
  return(out)
}
zf <- function(a,b){par(mfrow=c(a,b))}

zplotfitDF <- function(df=EachMetrics,x,y,xlab=NULL,ylab=NULL,addlegend=TRUE,col=1,pch=1,
                       ylimm=NULL,xlimm=NULL,illustrate=NULL,showmean=FALSE,force0int = FALSE,
                       hilite=1.0,linelwd=1,predinterval=FALSE,showlines=FALSE,showANAlines=FALSE){
  # hilite will color the legend if the R2 value is greater than hilite
  # par(mar=c(5,5,1,1))
  if(is.null(xlab))xlab <- substitute(x)
  if(is.null(ylab))ylab <- substitute(y) # deparse(substitute(y))
  usex <- df[,names(df)==x]
  usey <- df[,names(df)==y]
  if(!is.null(illustrate))useillustrate <- df[,names(df)==illustrate]
  if(!is.null(ylimm) & is.null(xlimm)){  plot(usex,usey,ylab=ylab,xlab=xlab,col=col,pch=pch,ylim=ylimm) } 
  if(is.null(ylimm) & !is.null(xlimm)){ plot(usex,usey,ylab=ylab,xlab=xlab,col=col,pch=pch,xlim=xlimm)  } 
  if(!is.null(ylimm) & !is.null(xlimm)){ plot(usex,usey,ylab=ylab,xlab=xlab,col=col,pch=pch,xlim=xlimm,ylim=ylimm) }
  if(is.null(ylimm) & is.null(xlimm)){  plot(usex,usey,ylab=ylab,xlab=xlab,col=col,pch=pch) }
  if(showlines)lines(usex,usey,lwd=0.5,col="grey50")
  # special hack for Anatone DF
  if(showANAlines){
  print("Using SPECIAL Hack for Anatone data frame will skip missing years in line drawning. Look in text to change this")
  lines(usex[1:7],usey[1:7],lwd=0.5,col="grey50")
  lines(usex[8:11],usey[8:11],lwd=0.5,col="grey50")
  lines(usex[12:18],usey[12:18],lwd=0.5,col="grey50")
  lines(usex[19:22],usey[19:22],lwd=0.5,col="grey50")
  lines(usex[23:52],usey[23:52],lwd=0.5,col="grey50")
  }
  if(force0int){
    fit <- lm(usey ~ 0 + usex)
    r2 <- summary(fit)$r.squared
    p <- summary(fit)$coefficients[1,4]
    slp <- summary(fit)$coefficients[1,1]
    int <- 0
    
  }
  else {fit <- lm(usey ~ usex,na.action=na.omit)
  r2 <- summary(fit)$r.squared
  p <- summary(fit)$coefficients[2,4]
  slp <- summary(fit)$coefficients[2,1]
  int <- summary(fit)$coefficients[1,1]
  newx <- (min(usex) - 1) : (max(usex) + 1)
  }
  
  # if(y != "F"){
  {
    abline(fit,col="darkred",lwd=linelwd)
    if(showmean) points(mean(usex,na.rm=T),mean(usey,na.rm=T),col="darkgreen",pch=16,cex=2)
    if(addlegend  )legend("topleft",bty="n",legend=c(paste("R2=", round(r2,4)),
                                                     paste("p=", signif(p,6)),
                                                     paste("slope=",signif(slp,3)),
                                                     ifelse(is.null(illustrate),"",paste0("\"",illustrate,"\""," terciles shown")),
                                                     ifelse(is.null(illustrate),"",paste0("     w/ \"X\" at medians"))
    ),text.col=ifelse(round(r2,2) >= hilite,"blue",col))
  }
  if(!is.null(illustrate)){ # necessariy has to be of same dimensions as x and y. uses x and y to locate points
    test <- quantile(useillustrate,probs=c(0.3333333,0.66666667),na.rm=TRUE)
    cols <- rep("grey30",length(usex));
    cols[useillustrate < test[1]] <- "blue";
    cols[useillustrate >= test[2]] <- "orange";
    points(usex,usey,col=cols,pch=16)
    points(median(usex[useillustrate < test[1]],na.rm=T),median(usey[useillustrate < test[1]],na.rm=T),col="blue",cex=1.5,pch="X",font=2)
    points(median(usex[useillustrate < test[2] & useillustrate >= test[1]],na.rm=T),
           median(usey[useillustrate < test[2] & useillustrate >= test[1]],na.rm=T),col="grey30",cex=1.5,pch="X",font=2)
    points(median(usex[useillustrate > test[2]],na.rm=T),median(usey[useillustrate > test[2]],na.rm=T),col="orange",cex=1.5,pch="X",font=2)
  }
  return(c(r2=r2,p=p,meanx=mean(usex,na.rm=T),meany=mean(usey,na.rm=T),int=int))
}

zAnaFit <- function(pav,daf,fixedpars,retype="LL",useairlag=1,myobs=NULL,myair=NULL){  
  # 8 paraemter version
  # will also do the weather as well as climate 
  # fits the air regime too
  # chop eliminates the beginning of year data for which we don't have the lags
  # note that use of myobs still requires a year with air temps, etc. and a substitutes pred for obs if forward looking as in a forecast
  pav <- as.numeric(pav)
  obs <- daf$ana.temp
  # receives FULL parameter vector
  aa <- bb <- 0
  aa <- zAnaForm(pav[1:7],daf,fixedpars)
  # tweak daf$Tair with the myair if it is not null
  if(!is.null(myair)){
    daf$Tair <- myair
    daf$Tair[is.na(myair)] <- zAnaFormAir(pv=fixedpars[c("m0","m1","n0","n1")],mydf=daf,fixedpars)[is.na(myair)]
  }
  bb <- zAnaResid(mypv=pav[8],daf,fixedpars,useairlag)
  if(!is.null(myobs)){obs <- myobs
    obs[is.na(myobs)] <- (aa+bb)[is.na(myobs)]
  }
 
  pred <- aa + bb
  
  diffs = obs-pred
  likes = dnorm(diffs, 0, 1)
  if(retype=="LL") { return (-sum(log(likes))) } 
  if(retype=="resid") { return(diffs) } else { return (pred) }
}

zAnaFit2.wsigma <- function(par,daf,fixedpars,retype="LL",useairlag=1,myobs=NULL,myair=NULL){  
  # notes with zAnaFit2
  par <- as.numeric(par)
  obs <- daf$ana.temp.fill
  aa <- bb <- 0
  # truncated parameter without the e0 and e1
  # in this wsigma, first par is log.sigma
  # ALL par[] refs have been adjusted
  aa <- zAnaForm2(pv=par[2:6],daf,fixedpars)
  # tweak daf$Tair with the myair if it is not null
  if(!is.null(myair)){
    daf$Tair <- myair
    daf$Tair[is.na(myair)] <- zAnaFormAir(pv=fixedpars[c("m0","m1","n0","n1")],mydf=daf,fixedpars)[is.na(myair)]
  }
  bb <- zAnaResid(mypv=par[7],daf,fixedpars,useairlag)
  if(!is.null(myobs)){obs <- myobs
  obs[is.na(myobs)] <- (aa+bb)[is.na(myobs)]
  }
  
  pred <- aa + bb
  diffs = obs-pred
  likes = dnorm(diffs, 0, exp(par[1]))
  if(retype=="LL") { return (-sum(log(likes))) } 
  if(retype=="resid") { return(diffs) } else { return (pred) }
}

zAnaFit2 <- function(pav,daf,fixedpars,retype="LL",useairlag=1,myobs=NULL,myair=NULL){  

  # Uses Anaform2 with fixed E
  pav <- as.numeric(pav)
  obs <- daf$ana.temp.fill
  # receives FULL parameter vector
  aa <- bb <- 0
  # truncated parameter without the e0 and e1
  aa <- zAnaForm2(pav[1:5],daf,fixedpars)
  # tweak daf$Tair with the myair if it is not null
  if(!is.null(myair)){
    daf$Tair <- myair
    daf$Tair[is.na(myair)] <- zAnaFormAir(pv=fixedpars[c("m0","m1","n0","n1")],mydf=daf,fixedpars)[is.na(myair)]
  }
  bb <- zAnaResid(mypv=pav[6],daf,fixedpars,useairlag)
  if(!is.null(myobs)){obs <- myobs
  obs[is.na(myobs)] <- (aa+bb)[is.na(myobs)]
  }

  pred <- aa + bb 
  diffs = obs-pred
  likes = dnorm(diffs, 0, 1)
  if(retype=="LL") { return (-sum(log(likes))) } 
  if(retype=="resid") { return(diffs) } else { return (pred) }
}

zAnaResid <- function(mypv,mydf,fixedpars,useairlag){
  # Only has the needed pars sent to it, not whole vector!
  myfullairfit <- zAnaFormAir(pv=fixedpars[c("m0","m1","n0","n1")],mydf,fixedpars)
  anom <- mypv[1]*rollmean(x=mydf$Tair - myfullairfit,k=useairlag,align=alignside,fill=NA)
  anom[is.na(anom)] <- mypv[1]*(mydf$Tair - myfullairfit)[is.na(anom)]
  return(anom)
}

zAnaFormAir <- function(pv,mydf,fixedpars){
  # Only has the 4 pars from the residual function sent to it
  pv <- as.numeric(pv)
  GM <- pv[1]  + pv[2]*mydf$gMeanAir
  RA <- ( pv[3] + pv[4]*(mydf$summerairJA - mydf$winterairDJ))
  days <- mydf$doy
  out <-   GM  + RA *sin(2*pi/365*(days + fixedpars['P']+(0+fixedpars['J']*sin(2*pi/365*(days+fixedpars['K'])))))  
  return(out)      
}
zAnaFormAirfit <- function(pv,mydf,fixedpars,retype="LL"){
  pv <- as.numeric(pv)
  obs <- mydf$Tair 
  pred <-  zAnaFormAir(pv,mydf,fixedpars)
  
  diffs = obs-pred
  likes = dnorm(diffs, 0, 1)
  if(retype=="LL") { return (-sum(log(likes))) } 
  if(retype=="resid") { return(diffs) } else { return (pred) }
}

zAnaFormAirfit.wsigma <- function(pv,mydf,fixedpars,retype="LL"){
  pv <- as.numeric(pv)
  obs <- mydf$Tair 
  pred <-  zAnaFormAir(pv[-1],mydf,fixedpars)
  
  diffs = obs-pred
  likes = dnorm(diffs, 0, exp(pv[1]))
  if(retype=="LL") { return (-sum(log(likes))) } 
  if(retype=="resid") { return(diffs) } else { return (pred) }
}


zAnaForm <- function(pv,df,fixedpars){
  aa <- bb <- cc <- 0
  aa <- (pv[1]  + pv[2]*df$gMeanAir)
  bb <- (pv[3] + pv[4]*(df$summerairJA - df$winterairDJ)) * sin((2*pi)/365 * (df$doy + fixedpars['CC']))
  onset <- (pv[5]  + pv[6]*df$springairA) ; 
  offset <- fixedpars['FF']
  cc0 <- pv[7]*log(1+ df$swe)
  cc1 <- (1 - sin(2*pi/(offset - onset)*(df$doy - onset + (offset-onset)/4)))
  cc <- cc0*cc1
  cc[ df$doy < onset | df$doy > offset] <- 0
  out <- aa + bb - cc
  return(out)
}
zAnaForm2 <- function(pv,df,fixedpars){
  aa <- bb <- cc <- 0
  aa <- (pv[1]  + pv[2]*df$gMeanAir)
  bb <- (pv[3] + pv[4]*(df$summerairJA - df$winterairDJ)) * sin((2*pi)/365 * (df$doy + fixedpars['CC']))
  onset <- fixedpars['EE']; 
  offset <- fixedpars['FF']
  cc0 <- pv[5]*log(1+ df$swe)
  cc1 <- (1 - sin(2*pi/(offset - onset)*(df$doy - onset + (offset-onset)/4)))
  cc <- cc0*cc1
  cc[ df$doy < onset | df$doy > offset] <- 0
  out <- aa + bb - cc
  return(out)
}

zWaterAuto <- function(pred,obs,ar1=0){
  out <- c(0,lag(obs-pred,1)[-1])*ar1 
  # remove the NA off the front and replace with 0
  # This shifts yesterdays error to be indexed by today!
  # and returns the AR correction to be added
  return(out)
}

#==== Results Wrapper function ====

z5 <- function(pdff=NULL,showw=TRUE,
               years=unique(ANA.DF2$year),
               myobs=NULL,
               myair=NULL,type=2,ar=FALSE){
  # type == 2 is going to use zAnafit2 requireing 2 fewer parameters t a fixed "E" or "EE"
  # uses ANA.DF2 directly
  # add AR option
  results <- NULL
  mod.obs.diffs <- NULL
  if(!is.null(pdff)){pdf(file=pdff,width=9.2,height=7) ; zf(1,1)}
  for(y in years){
    {
      legwords <- leglines <- legcols <- NULL
      K <- ANA.DF2$year==y 
      if(type==2){ 
        MD <- ANA.DF2$dailyFittype2[K]
        if(ar)MD <- ANA.DF2$dailyFittype2.a1[K]
      }else { 
        MD <- ANA.DF2$dailyFit[K]
      }
        Regime <- ANA.DF2$Regime[K]
      airOb <- ANA.DF2$Tair[K]
      Ob <- ANA.DF2$ana.temp[K]
      x <- ANA.DF2$doy[K]
      flow <- ANA.DF2$ana.flow[K]
      airRegime <- ANA.DF2$AirRegime[K]
      Q <- ANA.DF2$ana.flow[K]
      mad <- mean(abs(Ob-MD),na.rm=T)
      mre <- mean(Ob-MD,na.rm=TRUE)
      rmse <- sqrt(mean((MD-Ob)^2))
      diffs <- Ob - MD 
      cumerror <- sum(Ob - MD)
      cumsumerror <- cumsum(Ob-MD)
      worstcumsum <- max(abs(cumsumerror))
      R2 <- summary(lm(Ob ~ MD))$r.squared
      diff7days <-       quantile(abs(rollmean(MD,7,align=alignside,fill=NA)-rollmean(Ob,7,align=alignside,fill=NA)),na.rm=TRUE,probs=c(0.5,0.75,0.80,0.825,0.85,0.875,0.90,1))
      diff7days.type2 <- quantile(abs(rollmean(MD - Ob,7,align=alignside,fill=NA)),na.rm=TRUE,probs=c(0.5,0.75,0.80,0.825,0.85,0.875,0.90,1))
      summerMAD <- mean(abs(diffs[274:335]))
      winterMAD <- mean(abs(diffs[193:151]))
      results <- rbind.data.frame(results,cbind.data.frame(year=y,mad=mad,mre=mre,rmse=rmse,cumerror=cumerror,varExp=R2,
                                                           worstcumsum=worstcumsum,error7max=diff7days[["100%"]],
                                                           error7.90=diff7days[["90%"]],error7.80=diff7days[["80%"]],
                                                           error7.50=diff7days[["50%"]],summerMAD,winterMAD))
      mod.obs.diffs <- c(mod.obs.diffs,diffs)
    } 
    if(showw){
      par(mar=c(4,4,1,3),cex=1.25,las=1)
      axiswordsize <- 1.5
      axisnumsize <- 1.2
      plot(x,MD,type="n",xlab="",ylab="",ylim=c(-15,25),axes=F);box();zxWYdates();
      axis(2,at=c(5,10,15,20,25),las=1,cex.axis=axisnumsize); mtext(bquote("Water"~ degree*"C"),2,adj=0.7,line=2.5,las=3,cex=axiswordsize)
      lines(x,Ob,col="palegreen",lwd=3)
      legwords <- c(legwords,"Water temperature data");legcols <- c(legcols,"palegreen");leglines <- c(leglines,3)
      lines(x,Regime,col="blue",lwd=3)
      legwords <- c(legwords,paste("Thermal Regime"));legcols <- c(legcols,"blue");leglines <- c(leglines,2)
      lines(x,MD,col="darkgreen",lwd=3)
      legwords <- c(legwords,paste("Temperature Model"));legcols <- c(legcols,"darkgreen");leglines <- c(leglines,3)
      
      
      lineoffset <- 10
      lines(x,diffs - lineoffset,col="brown",lwd=1);
      legwords <- c(legwords,"Data - Model difference");leglines <- c(leglines,1);legcols <- c(legcols,"brown")
      axis(2,at=c(-2,0,2)-lineoffset,labels=c(-2,0,2),cex.axis=axisnumsize,col.axis="brown",las=1) ;
      mtext(bquote("Data - Model (error)"~~degree*"C"),2,line=2.5,las=0,adj=0.1,col="brown",cex=axiswordsize); segments(-1,-lineoffset,par()$usr[2])
      
      lines(x,airOb/1.5-5,col="grey80")
      legwords <- c(legwords,"Air temperature data");leglines <- c(leglines,1);legcols <- c(legcols,"grey80")
      airlab <- c(-5,0,5,10,15,20,25,30,35)
      axis(4,at=airlab/1.5-5,labels=airlab,cex.axis=axisnumsize,col.axis="grey50",las=1,line=0) ;
      lines(x,airRegime/1.5-5,lwd=3,col="grey70")
      mtext(bquote("Air"~ degree*"C"),4,line=1,las=0,adj=0.03,col="grey50",cex=axiswordsize)
      legwords <- c(legwords,"Air Regime");leglines <- c(leglines,3);legcols <- c(legcols,"grey70")
      legend("topleft",bty="n",legend=y,cex=1.2)
      legend("topleft",inset=c(0.2,0),bty="n",legend=legwords,col=legcols,lwd=leglines,cex=0.9)
      
    }   
  }
  if(!is.null(pdff)){print(paste0("PDF: ","v2.pdf")) ; dev.off()}
  return(results)
}

