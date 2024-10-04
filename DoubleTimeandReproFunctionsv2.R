#
# File has 5 functions:
#     1) classic daily doubling time calculator for count data, 
#     2) a version for weekly data for count data
#     3) a version to fit to deterministic models for testing or rate data
#     4) the Max likelihood estimator for beta-binomial model 
#     5) Beta-binomial simulator
#     6) gamTVRpred a function that calculates the time varying reproduction number
#     7) simulator function for SIS infection
#     8) simulator function for SEIR infection
#     9) wrapper function for post processing simulations with gamTVRpred

##################
## function for running GAM and creating plots
library(ggplot2)
#library(ggpubr)

DoubleTime <- function(dat, timev, aggregate = 14, npts=200, plt=FALSE, subgroup=FALSE, figtitle=""){
  kval <- floor(length(dat)/20)
  kval <- ifelse(kval<10, 10, kval)
  res<- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts))
  Tv <- timev[1:(length(timev))] 
  DW <- weekdays(as.Date(Tv, origin = "1899-12-30"))
  datfull <- dat
  dat <- dat[1:(length(dat))]
  
  MGAM <- gam(dat~s(Tv, bs='gp', k=kval)+DW, family=nb)
  
# cheap derivative, modified from a weblink to mgcv help files shared by Simon Wood March 2020.
  dow <- rep('Friday', length(timev)) #weekdays(as.Date(xv1, origin = "1899-12-30"))
  newd <- data.frame(Tv=timev, DW=dow) # data.frame(Tv=xv, DW=DW)
  p <- predict(MGAM, newd, type = "link", se.fit = TRUE)
  meanspline <- (exp(mean(c(0,coef(MGAM)[2:7]))+p$fit))
  
  xv<-seq(min(Tv),max(Tv), length=npts)
  dow <- weekdays(as.Date(xv, origin = "1899-12-30"))
  newd <- data.frame(Tv=xv, DW=dow) # data.frame(Tv=xv, DW=DW)
  X0 <- predict(MGAM, newd,type="lpmatrix")
  
  eps <- 1e-7 ## finite difference interval
  xv<-seq(min(Tv),max(Tv), length=npts)+eps
  dow <- weekdays(as.Date(xv, origin = "1899-12-30"))
  newd <- data.frame(Tv=xv, DW=dow) # data.frame(Tv=xv, DW=DW)
  X1 <- predict(MGAM, newd,type="lpmatrix")
  newd <- data.frame(Tv=xv+eps, DW=dow) # data.frame(Tv=xv, DW=DW)
  X2 <- predict(MGAM, newd,type="lpmatrix")

  Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
  Xp2 <- (X2-2*X1+X0)/eps^2 ## maps coefficients to (fd approx.) derivatives
  Xi <- Xp*0 
  Xi[,1:(kval-1)+7] <- Xp[,1:(kval-1)+7] ## Xi%*%coef(MGAM) = smooth deriv i
  Xi2 <- Xp2*0 
  Xi2[,1:(kval-1)+7] <- Xp2[,1:(kval-1)+7] ## Xi%*%coef(MGAM) = smooth deriv i
  df <- Xi%*%coef(MGAM)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%MGAM$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  df2 <- Xi2%*%coef(MGAM)              ## ith smooth derivative 
  df2.sd <- rowSums(Xi2%*%MGAM$Vp*Xi2)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  res$sdt <- df
  res$sdtup <- df+2*df.sd
  res$sdtlow <- df-2*df.sd
  res$time <- as.Date(xv, origin = "1899-12-30")
# CI is approximate.
  res$sdt2 <- df2
  res$sdt2up <- df2+2*df2.sd
  res$sdt2low <- df2-2*df2.sd
  #  MGLM <- glm(dat~(Tv), family=quasipoisson)
  #  Doubling<-c(MGLM$coefficients[2],confint(MGLM)[2,1],confint(MGLM)[2,2])
  
  if(plt==TRUE){

    xv1<-c(timev, max(timev)+1:14)#timev #seq(min(Tv),max(Tv), length=npts)+eps
    dow <- weekdays(as.Date(xv1, origin = "1899-12-30"))
    newd <- data.frame(Tv=xv1, DW=dow) # data.frame(Tv=xv, DW=DW)
    p <- predict(MGAM, newd, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGAM$family$linkinv(upr)
    lwr <- MGAM$family$linkinv(lwr)
    dattmp <- data.frame(xval = as.Date(xv1, origin = "1899-12-30"), central = exp(p$fit), 
                         upper= upr, lower=lwr, obs = c(datfull, rep(NA, 14)), meanspl = c(meanspline, rep(NA, 14)))

    nsim <-1000
    set.seed(42)
    simCIlwr<-apply(array(rnbinom(n=nsim*aggregate, mu=(lwr[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE)), c((aggregate),nsim)),2,sum)
    simCIupr<-apply(array(rnbinom(n=nsim*aggregate, mu=(upr[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE)), c((aggregate),nsim)),2,sum)
    PredInt <- c(round(sum(exp(p$fit)[length(p$fit)-(aggregate-1):0])), 
                     floor(quantile(simCIlwr, c(0.025), na.rm=T)),ceiling(quantile(simCIupr,0.975, na.rm=T)))       
    resprint <- c(PredInt)
        
    p1 <- ggplot(data=dattmp, aes(x=xval, y=central)) + 
          geom_point(aes(x=xval, y= obs), color='red') +
          geom_ribbon(aes(ymin=lower, ymax=upper), fill='lightgreen')+#+ 
          geom_line(color='green') + 
          geom_line(aes(x=xval, y=meanspl), color='blue') +
          labs(x='Date', y='Cases', title=figtitle, 
               caption=paste('Data (red points) from UKHSA \n',   
                             'Projected new events in next ', aggregate, ' days: ',round(resprint[1],0),' (',round(resprint[2],0),',',round(resprint[3],0),')', sep=''))
        
        if(df[npts]-2*df.sd[npts]>0){
          resstring <- paste('Increasing recent trend', sep='')
        }else if(df[npts]+2*df.sd[npts]<0){
          resstring <- paste('Decreasing recent trend', sep='')
        }else{
          if(df[npts]>0){
            resstring <- paste('Plateauing recent trend, may be increasing with probability ',round(1-pnorm(0, df[npts],df.sd[npts]), 3), sep='')
          }else{
            resstring <- paste('Plateauing recent trend, may be decreasing with probability ',round(pnorm(0, df[npts],df.sd[npts]),3), sep='')
          }
        }
        
    p2 <- ggplot(data=res, aes(x=time, y=sdt))+
          geom_ribbon(aes(ymin=sdtlow, ymax=sdtup), fill='lightblue')+#+ 
          geom_line(color='darkblue') +
          labs(x='Date', title='Derivative of spline arising from GAM', caption=resstring)+
          scale_y_continuous(
            name = expression("Instantaneous growth rate"), 
            sec.axis = sec_axis(~., name = "Doubling Time", 
                                breaks = c(-0.1,-0.05,-0.033,-0.0247, 0,0.0247, 0.033, 0.0495, 0.099, 0.173, 0.3466), 
                                labels=c(-7,-14,-21,-28, 'Infinite',28,21, 14, 7,4,2)), 
            limits = c(min(res$sdtlow), max(res$sdtup)))   # theme(
    print(p1+p2+plot_layout(ncol=1))
       # print(ggarrange( p1, p2, ncol=1))  
        
    if(subgroup==FALSE){
      res <- list(splinederiv=res, PI=PredInt, splinedaymean=meanspline )
    }else{
        res <- list(splinederiv=res, PI=PredInt, simlow =simCIlwr, simupper =simCIupr, splinedaymean=meanspline )
    }
          
  }
  
  res
}




## function for running GAM and creating plots
DoubleTimeWeek <- function(dat, timev=seq(1,7*length(dat),7), aggregate = 2, npts=200, plt=FALSE, figtitle=""){
  kval <- floor(length(dat)/3) #*7/20)
  res<- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts))
  Tv <- timev[1:(length(timev))] 
  datfull <- dat
  dat <- dat[1:(length(dat))]
  
#  MGAM <- gam(dat~s(Tv, bs='gp', k=kval), family=poisson)
  MGAM <- gam(dat~s(Tv, bs='gp', k=kval), family=nb)
  
  xv<-seq(min(Tv),max(Tv), length=npts)
  newd <- data.frame(Tv=xv) # data.frame(Tv=xv, DW=DW)
  X0 <- predict(MGAM, newd,type="lpmatrix")
  
  eps <- 1e-7 ## finite difference interval
  xv<-seq(min(Tv),max(Tv), length=npts)+eps
  newd <- data.frame(Tv=xv) # data.frame(Tv=xv, DW=DW)
  X1 <- predict(MGAM, newd,type="lpmatrix")
  Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
  
  Xi <- Xp*0 
  Xi[,1:(kval-1)+1] <- Xp[,1:(kval-1)+1] ## Xi%*%coef(MGAM) = smooth deriv i
  df <- Xi%*%coef(MGAM)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%MGAM$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  res$sdt <- df
  res$sdtup <- df+2*df.sd
  res$sdtlow <- df-2*df.sd
  res$time <- as.Date(xv, origin = "1899-12-30")
  
  #  MGLM <- glm(dat~(Tv), family=quasipoisson)
  #  Doubling<-c(MGLM$coefficients[2],confint(MGLM)[2,1],confint(MGLM)[2,2])
  
  if(plt==TRUE){
    par(mfrow=c(2,1))
    
    plot(as.Date(c(timev), origin = "1899-12-30"), c(datfull), xlim=as.Date(c(min(timev), max(timev)+2), origin = "1899-12-30"), pch=16, col=2, main='', ylab='Number', xlab='Time', xaxt='n')
    axis.Date(1, at=as.Date(c('2020-03-01','2020-04-01','2020-05-01','2020-06-01','2020-07-01','2020-08-01','2020-09-01','2020-10-01','2020-11-01','2020-12-01','2021-01-01','2021-02-01','2021-03-01','2021-04-01','2021-05-01','2021-06-01','2021-07-01','2021-08-01','2021-09-01','2021-10-01','2021-11-01','2021-12-01','2022-01-01','2022-02-01','2022-03-01','2022-04-01','2022-05-01','2022-06-01','2022-07-01','2022-08-01','2022-09-01','2022-10-01','2022-11-01','2022-12-01')), format='%m-%y')
    
    xv1<-c(timev, max(timev)+1:2)#timev #seq(min(Tv),max(Tv), length=npts)+eps
 #   dow <- weekdays(as.Date(xv1, origin = "1899-12-30"))
    newd <- data.frame(Tv=xv1) # data.frame(Tv=xv, DW=DW)
    p <- predict(MGAM, newd, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGAM$family$linkinv(upr)
    lwr <- MGAM$family$linkinv(lwr)
    lines(as.Date(xv1, origin = "1899-12-30"), exp(p$fit))
    lines(as.Date(xv1, origin = "1899-12-30"), upr, col=1, lty=2)
    lines(as.Date(xv1, origin = "1899-12-30"), lwr, col=1, lty=2)
    abline(v=as.Date("2021-02-25"))
    abline(v=as.Date("2021-05-11"))
    abline(v=as.Date("2021-12-18"))
    polygon(x = c(as.Date(xv1, origin = "1899-12-30"), rev(as.Date(xv1, origin = "1899-12-30"))),
            y = c(upr, 
                  rev(lwr)),
            col =  adjustcolor("blue", alpha.f = 0.10), border = NA)  
    #   text(as.Date(min(xv1)+40, origin = "1899-12-30"), 3/4*max(datfull), pos=4,cex=0.75, paste('Projected new deaths \n in next 2 weeks:\n',round(sum(exp(p$fit)[length(p$fit)-13:0])),' (',
    #                                                                qnbinom(c(0.025*14), mu=sum(lwr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)),',',
    #                                                               qnbinom(c(1-0.025*14), mu=sum(upr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)),')', sep=''))
    #    res <- c(round(sum(exp(p$fit)[length(p$fit)-13:0])), 
    #            qnbinom(c(0.025*14), mu=sum(lwr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)),
    #           qnbinom(c(1-0.025*14), mu=sum(upr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)))
    PredInt <- cbind((rnbinom(1000, mu=sum(exp(p$fit)[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE))),
                 (rnbinom(1000, mu=sum(lwr[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE))),
                 (rnbinom(1000, mu=sum(upr[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE))))
    #print((is.finite(res[,3])))
    resprint <- c(quantile(PredInt[,1], probs = 0.5),quantile(PredInt[,2], probs = 0.025),
                  ifelse(sum(is.finite(PredInt[,3]))==0,NA,quantile(PredInt[,3], probs = 0.975)))
    res <- list(splinederiv=res, PI=resprint)
    #    text(as.Date(min(xv1)+70, origin = "1899-12-30"), 3/4*max(datfull), pos=4,cex=0.75, paste('Projected new deaths \n in next 2 weeks:\n',round(resprint[1],0),' (',round(resprint[2],0),',',round(resprint[3],0),')', sep=''))
    #                                                                qnbinom(c(0.025*14), mu=sum(lwr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)),',',
    #                                                               qnbinom(c(1-0.025*14), mu=sum(upr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)),')', sep=''))
    
    plot(as.Date(xv, origin = "1899-12-30"),df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)), ylab='Instantaneous growth rate', xlab='Time', main='', xaxt='n')
    axis.Date(1, at=as.Date(c('2020-03-01','2020-04-01','2020-05-01','2020-06-01','2020-07-01','2020-08-01','2020-09-01','2020-10-01','2020-11-01','2020-12-01','2021-01-01','2021-02-01','2021-03-01','2021-04-01','2021-05-01','2021-06-01','2021-07-01','2021-08-01','2021-09-01','2021-10-01','2021-11-01','2021-12-01','2022-01-01','2022-02-01','2022-03-01','2022-04-01','2022-05-01','2022-06-01','2022-07-01','2022-08-01','2022-09-01','2022-10-01','2022-11-01','2022-12-01')), format='%m-%y')
    polygon(x = c(c(as.Date("2020-03-23"), as.Date("2020-06-01")),rev(c(as.Date("2020-03-23"), as.Date("2020-06-01")))), y = c(0.5,0.5,-0.5,-0.5), col = adjustcolor(rgb(0.8,0.8,0.8), alpha.f = 0.50), border = NA)
    polygon(x = c(c(as.Date("2020-10-31"), as.Date("2020-12-02")),rev(c(as.Date("2020-10-31"), as.Date("2020-12-02")))), y = c(0.5,0.5,-0.5,-0.5), col = adjustcolor(rgb(0.8,0.8,0.8), alpha.f = 0.50), border = NA)
    polygon(x = c(c(as.Date("2021-01-04"), as.Date("2021-03-08")),rev(c(as.Date("2021-01-04"), as.Date("2021-03-08")))), y = c(0.5,0.5,-0.5,-0.5), col = adjustcolor(rgb(0.8,0.8,0.8), alpha.f = 0.50), border = NA)
#    polygon(x = c(c(as.Date("2020-12-20"), as.Date("2021-02-15")),rev(c(as.Date("2020-12-20"), as.Date("2021-02-15")))), y = c(0.5,0.5,-0.5,-0.5), col = adjustcolor(rgb(1,0.9,0.9), alpha.f = 0.50), border = NA)
 #   polygon(x = c(c(as.Date("2021-03-16"), as.Date("2021-05-10")),rev(c(as.Date("2021-03-16"), as.Date("2021-05-10")))), y = c(0.5,0.5,-0.5,-0.5), col = adjustcolor(rgb(1,0.8,0.8), alpha.f = 0.50), border = NA)
abline(v=as.Date("2021-02-25"))
abline(v=as.Date("2021-05-11"))
abline(v=as.Date("2021-12-18"))
lines(as.Date(xv, origin = "1899-12-30"),df+2*df.sd,lty=2);
    lines(as.Date(xv, origin = "1899-12-30"),df-2*df.sd,lty=2)
    polygon(x = c(as.Date(xv, origin = "1899-12-30"), rev(as.Date(xv, origin = "1899-12-30"))),
            y = c(df+2*df.sd, 
                  rev(df-2*df.sd)),
            col =  adjustcolor("blue", alpha.f = 0.10), border = NA)  
    abline(h=0, col=4)
    axis(4,at=c(-0.1,-0.05,-0.033,-0.0247, 0,0.0247, 0.033, 0.0495, 0.099, 0.173, 0.3466), labels=c(c(-7,-14,-21,-28)/7, 'Infinite',c(28,21, 14, 7,4,2)/7))
    mtext(figtitle, outer=TRUE,  cex=1.5, line=-4)
    mtext(paste('Projected new events in next ', aggregate, ' weeks: ',round(resprint[1],0),' (',round(resprint[2],0),',',round(resprint[3],0),')', sep=''), side=3, line=2)
    if(df[npts]-2*df.sd[npts]>0){
      mtext(paste('Increasing recent trend', sep=''), side=3, line=3)
    }else if(df[npts]+2*df.sd[npts]<0){
      mtext(paste('Decreasing recent trend', sep=''), side=3, line=3)
    }else{
      if(df[npts]>0){
        mtext(paste('Plateauing recent trend, may be increasing', sep=''), side=3, line=3)
      }else{
        mtext(paste('Plateauing recent trend, may be decreasing', sep=''), side=3, line=3)
      }
    }
  }
  
  res
}

## function for running GAM and creating plots
DoubleTimeRate <- function(dat, timev=seq(1,length(dat),1), aggregate = 2, npts=200, plt=FALSE, figtitle=""){
  kval <- floor(length(dat)/20)
  res<- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts))
  Tv <- timev[1:(length(timev))] 
  datfull <- dat
  dat <- dat[1:(length(dat))]
  
  MGAM <- gam(dat~s(Tv, bs='gp', k=kval))
  

  xv<-seq(min(Tv),max(Tv), length=npts)
  newd <- data.frame(Tv=xv) # data.frame(Tv=xv, DW=DW)
  X0 <- predict(MGAM, newd,type="lpmatrix")
#  p <- predict(MGAM,newd,  type = "link", se.fit = TRUE)
 # meanspline <- (exp(p$fit))
  
  eps <- 1e-7 ## finite difference interval
  xv<-seq(min(Tv),max(Tv), length=npts)+eps
  newd <- data.frame(Tv=xv) # data.frame(Tv=xv, DW=DW)
  X1 <- predict(MGAM, newd,type="lpmatrix")
  Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
  
  Xi <- Xp*0 
  Xi[,1:(kval-1)+1] <- Xp[,1:(kval-1)+1] ## Xi%*%coef(MGAM) = smooth deriv i
  df <- Xi%*%coef(MGAM)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%MGAM$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  res$sdt <- df
  res$sdtup <- df+2*df.sd
  res$sdtlow <- df-2*df.sd
  
  #  MGLM <- glm(dat~(Tv), family=quasipoisson)
  #  Doubling<-c(MGLM$coefficients[2],confint(MGLM)[2,1],confint(MGLM)[2,2])
  
  if(plt==TRUE){
    par(mfrow=c(2,1))
    
    plot(as.Date(c(timev), origin = "1899-12-30"), c(datfull), xlim=as.Date(c(min(timev), max(timev)+2), origin = "1899-12-30"), pch=16, col=2, main='', ylab='Number', xlab='Time', xaxt='n')
    axis.Date(1, at=as.Date(c('2020-03-01','2020-04-01','2020-05-01','2020-06-01','2020-07-01','2020-08-01','2020-09-01','2020-10-01','2020-11-01','2020-12-01','2021-01-01','2021-02-01','2021-03-01','2021-04-01','2021-05-01','2021-06-01','2021-07-01','2021-08-01','2021-09-01','2021-10-01','2021-11-01','2021-12-01','2022-01-01','2022-02-01','2022-03-01','2022-04-01','2022-05-01','2022-06-01','2022-07-01','2022-08-01','2022-09-01','2022-10-01','2022-11-01','2022-12-01')), format='%m-%y')
    
    xv1<-c(timev, max(timev)+1:2)#timev #seq(min(Tv),max(Tv), length=npts)+eps
    #   dow <- weekdays(as.Date(xv1, origin = "1899-12-30"))
    newd <- data.frame(Tv=xv1) # data.frame(Tv=xv, DW=DW)
    p <- predict(MGAM, newd, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGAM$family$linkinv(upr)
    lwr <- MGAM$family$linkinv(lwr)
  #  abline(h=6, col='grey')
    abline(h=c(3,6,9), col='grey')
    lines(as.Date(xv1, origin = "1899-12-30"), (p$fit))
    lines(as.Date(xv1, origin = "1899-12-30"), upr, col=1, lty=2)
    lines(as.Date(xv1, origin = "1899-12-30"), lwr, col=1, lty=2)
    polygon(x = c(as.Date(xv1, origin = "1899-12-30"), rev(as.Date(xv1, origin = "1899-12-30"))),
            y = c(upr, 
                  rev(lwr)),
            col =  adjustcolor("blue", alpha.f = 0.10), border = NA)  
    #   text(as.Date(min(xv1)+40, origin = "1899-12-30"), 3/4*max(datfull), pos=4,cex=0.75, paste('Projected new deaths \n in next 2 weeks:\n',round(sum(exp(p$fit)[length(p$fit)-13:0])),' (',
    #                                                                qnbinom(c(0.025*14), mu=sum(lwr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)),',',
    #                                                               qnbinom(c(1-0.025*14), mu=sum(upr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)),')', sep=''))
    #    res <- c(round(sum(exp(p$fit)[length(p$fit)-13:0])), 
    #            qnbinom(c(0.025*14), mu=sum(lwr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)),
    #           qnbinom(c(1-0.025*14), mu=sum(upr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)))
#    res <- cbind((rnbinom(1000, mu=sum(exp(p$fit)[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE))),
 #                (rnbinom(1000, mu=sum(lwr[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE))),
  #               (rnbinom(1000, mu=sum(upr[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE))))
    #print((is.finite(res[,3])))
   # resprint <- c(quantile(res[,1], probs = 0.5),quantile(res[,2], probs = 0.025),
     #             ifelse(sum(is.finite(res[,3]))==0,NA,quantile(res[,3], probs = 0.975)))
    ##    text(as.Date(min(xv1)+70, origin = "1899-12-30"), 3/4*max(datfull), pos=4,cex=0.75, paste('Projected new deaths \n in next 2 weeks:\n',round(resprint[1],0),' (',round(resprint[2],0),',',round(resprint[3],0),')', sep=''))
    #                                                                qnbinom(c(0.025*14), mu=sum(lwr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)),',',
    #                                                               qnbinom(c(1-0.025*14), mu=sum(upr[length(p$fit)-13:0]), size=MGAM$family$getTheta(TRUE)),')', sep=''))
    res <- list(splinederiv=res, splinedaymean=p$fit )
    
    plot(as.Date(xv, origin = "1899-12-30"),df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)), ylab='Instantaneous growth rate', xlab='Time', main='', xaxt='n')
    axis.Date(1, at=as.Date(c('2020-03-01','2020-04-01','2020-05-01','2020-06-01','2020-07-01','2020-08-01','2020-09-01','2020-10-01','2020-11-01','2020-12-01','2021-01-01','2021-02-01','2021-03-01','2021-04-01','2021-05-01','2021-06-01','2021-07-01','2021-08-01','2021-09-01','2021-10-01','2021-11-01','2021-12-01','2022-01-01','2022-02-01','2022-03-01','2022-04-01','2022-05-01','2022-06-01','2022-07-01','2022-08-01','2022-09-01','2022-10-01','2022-11-01','2022-12-01')), format='%m-%y')
    polygon(x = c(c(as.Date("2020-03-23"), as.Date("2020-06-01")),rev(c(as.Date("2020-03-23"), as.Date("2020-06-01")))), y = c(0.5,0.5,-0.5,-0.5), col = adjustcolor(rgb(0.8,0.8,0.8), alpha.f = 0.50), border = NA)
    polygon(x = c(c(as.Date("2020-10-31"), as.Date("2020-12-02")),rev(c(as.Date("2020-10-31"), as.Date("2020-12-02")))), y = c(0.5,0.5,-0.5,-0.5), col = adjustcolor(rgb(0.8,0.8,0.8), alpha.f = 0.50), border = NA)
    polygon(x = c(c(as.Date("2021-01-04"), as.Date("2021-03-08")),rev(c(as.Date("2021-01-04"), as.Date("2021-03-08")))), y = c(0.5,0.5,-0.5,-0.5), col = adjustcolor(rgb(0.8,0.8,0.8), alpha.f = 0.50), border = NA)
    polygon(x = c(c(as.Date("2020-12-20"), as.Date("2021-02-15")),rev(c(as.Date("2020-12-20"), as.Date("2021-02-15")))), y = c(0.5,0.5,-0.5,-0.5), col = adjustcolor(rgb(1,0.9,0.9), alpha.f = 0.50), border = NA)
    polygon(x = c(c(as.Date("2021-03-16"), as.Date("2021-05-10")),rev(c(as.Date("2021-03-16"), as.Date("2021-05-10")))), y = c(0.5,0.5,-0.5,-0.5), col = adjustcolor(rgb(1,0.8,0.8), alpha.f = 0.50), border = NA)
    lines(as.Date(xv, origin = "1899-12-30"),df+2*df.sd,lty=2);
    lines(as.Date(xv, origin = "1899-12-30"),df-2*df.sd,lty=2)
    polygon(x = c(as.Date(xv, origin = "1899-12-30"), rev(as.Date(xv, origin = "1899-12-30"))),
            y = c(df+2*df.sd, 
                  rev(df-2*df.sd)),
            col =  adjustcolor("blue", alpha.f = 0.10), border = NA)  
    abline(h=0, col=4)
    axis(4,at=c(-0.1,-0.05,-0.033,-0.0247, 0,0.0247, 0.033, 0.0495, 0.099, 0.173, 0.3466), labels=c(c(-7,-14,-21,-28)/7, 'Infinite',c(28,21, 14, 7,4,2)/7))
    mtext(figtitle, outer=TRUE,  cex=1.5, line=-4)
#    mtext(paste('Projected new deaths in next ', aggregate, ' days: ',round(resprint[1],0),' (',round(resprint[2],0),',',round(resprint[3],0),')', sep=''), side=3, line=2)
    if(df[npts]-2*df.sd[npts]>0){
      mtext(paste('Increasing recent trend', sep=''), side=3, line=3)
    }else if(df[npts]+2*df.sd[npts]<0){
      mtext(paste('Decreasing recent trend', sep=''), side=3, line=3)
    }else{
      if(df[npts]>0){
        mtext(paste('Plateauing recent trend, may be increasing', sep=''), side=3, line=3)
      }else{
        mtext(paste('Plateauing recent trend, may be decreasing', sep=''), side=3, line=3)
      }
    }
  }
  
  res
}

lpbbinom <- function(para, kval, nval){
  aval <- para[1]*para[2]
  bval <- (1-para[1])*para[2]
  lnf <- lgamma(kval+aval)+lgamma(nval-kval+bval)-lgamma(nval+aval+bval)+
          lgamma(aval+bval)-lgamma(aval)-lgamma(bval)+lgamma(nval+1)-lgamma(kval+1)-lgamma(nval-kval+1)
  LL <- sum(lnf)
  LL
}

rbbinom <- function(nsim=1000, para=c(0.1,10), size = 50){
  # parameter size could be a vector of length nsim.
  s1<-para[1]*para[2]
  s2<-(1-para[1])*para[2]
#  print(c(s1,s2))
  betasamp <- rbeta(nsim, shape1=s1, shape2=s2)
  binomsamp<- rbinom(nsim, size=size, prob=betasamp)
  binomsamp
}

# load packages for later 
library(mgcv)
library(grDevices)
library(extraDistr)
library(MASS)

library(ggplot2)
library(patchwork)
library(EpiEstim)
library(incidence)

# Define functions

# gampred()
# Output: derives a prediction of central estimate from a model, 95% CrI and produces plots.
# Arguments: 
gamTVRpred <- function(data, pval=5, plt=TRUE, pltgr=FALSE, wrtfile=FALSE, type='SEIR', para=c(300,450000, 5,3), EEwindow=floor(2*(para[3]+para[4])),
                       figtitle = '', simulatorR=0,
                       keydates = data.frame(start=c("2020-03-23", "2020-11-05","2021-01-05"), end=c( "2020-05-31","2020-12-01","2021-03-08"))){
  start_time <- Sys.time()
  zv <- qnorm(1-pval/200,0,1) # CI will be asymptotic and so calcuate z-score
  npts <- length(data[,1]) 
  res<- data.frame(time=data[,1], vals=data[,2]) # set up a new data.frame for results, 
            # this assumes input is a 2 column matrix with first column time points and second column values to be fitted  
  res$Date=as.Date(res$time,origin = "1969-12-31") # One could set this up to use dates through but...
            # ... this is less fiddly with changes to date format in raw data though requires instead input to be have origin of 1st Jan 1970
  res$weekday=weekdays(res$Date)
  #print(mean(data[,2]))
  res$time=res$time-min(res$time) # set the time start at 0
  if(type=='SIS'){
    knot=floor(npts/(para[3]+para[4])) # set knots to be 1 generations
  }else if(type=='SEIR'){
    knot=floor(npts/(para[3]+para[4])/2) # set knots to be 2-3 generations
  }else if(type=='SEIRWeek'){
    knot=floor(npts/(para[3]+para[4])/2*7) # set knots to be 2-3 generations
  }else{
    knot=10
    }
  knot <- ifelse(knot<10, 10, knot) # unless time series is short, this method may not be reiable in this situation
  
  # run gam
  GAM=gam(vals~s(time, bs='ps',k=knot)+weekday, family=nb(link=log),data = res)
  # if desired print out summary values for later reporting
  if(wrtfile==TRUE){
    print(summary(GAM))
    print(exp(summary(GAM)$p.table[,1]))
    GAM_weekday= summary(GAM)$p.table
    GAM_weekday=as.data.frame(round(GAM_weekday,3))
    write.table(GAM_weekday, file = "GAM_weekday.txt",sep = "&",qmethod='double',quote = F)
  }
  
  # predict GAM on variable scale
  res$fullgam <- fitted(GAM)

  # predict GAM average over week day
  nt <- res$time
  Xp=predict(GAM,newdata=data.frame(time=nt,weekday="Friday"),type="lpmatrix")  #[1:(knot-1)+7]
  res$meanspline <- as.numeric(Xp%*%coef(GAM))+sum(coef(GAM)[2:7])/7
  # and infer cheap CI (but with unconditional variance).
  dfsd <- rowSums(Xp%*%GAM$Vc*Xp)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  res$dfse <- dfsd
  res$msup <- res$meanspline+zv*dfsd
  res$mslow <- res$meanspline-zv*dfsd
  # now predict on response scale
  res$pred <- exp(res$meanspline)
#  print(mean(res$pred))
  # CI firstly assuming linear response and conditional variance
  Xpse=predict(GAM,newdata=data.frame(time=nt,weekday="Friday"),type="response", se=TRUE)  #[1:(knot-1)+7]
  res$linearSE <- Xpse$se.fit
  res$predLup <- res$pred+zv*Xpse$se.fit
  res$predLlow <- res$pred-zv*Xpse$se.fit
  # then the CI assuming linear but with unconditional variance
  res$predULup <- exp(res$msup)
  res$predULlow <- exp(res$mslow)
  # then check for non-linearity issues (noting this is more computationally intensive)
  betasim <- mvrnorm(10000,coef(GAM),vcov(GAM, unconditional=TRUE))
  qntsim <- apply(exp(Xp %*% t(betasim)), 1, quantile, probs=c(0.025, 0.975))
  res$predUNLup <- qntsim[2,]*exp(sum(coef(GAM)[2:7])/7)
  res$predUNLlow <- qntsim[1,]*exp(sum(coef(GAM)[2:7])/7)

  # then calcuate the derivative for double time calcs
  eps <- 1e-3 ## finite difference interval
  X1 <- predict(GAM, data.frame(time=nt+eps,weekday="Friday"),type="lpmatrix")
  Xdp <- (X1-Xp)/eps ## maps coefficients to (fd approx.) derivatives

  newd <- data.frame(time=nt+2*eps,weekday="Friday") # data.frame(Tv=xv, DW=DW)
  X2 <- predict(GAM, newd,type="lpmatrix")
  Xdp2 <-  (Xp-2*X1+X2)/eps^2 ## maps coefficients to (fd approx.) derivatives
  newd <- data.frame(time=nt+3*eps,weekday="Friday") # data.frame(Tv=xv, DW=DW)
  X3 <- predict(GAM, newd,type="lpmatrix")
  Xdp3 <-  (X3-3*X2+3*X1-Xp)/eps^3 ## maps coefficients to (fd approx.) derivatives
  # Xdp2 <- ((Xdp2-Xdp)/eps)
  #  Xdp <- (X1-Xp)/eps ## maps coefficients to (fd approx.) derivatives
  Xi <- Xdp*0 
  Xi[,1:(knot-1)+7] <- Xdp[,1:(knot-1)+7] ## Xi%*%coef(MGAM) = smooth deriv i
  Xi2 <- Xdp2*0 
  Xi2[,1:(knot-1)+7] <- Xdp2[,1:(knot-1)+7] ## Xi%*%coef(MGAM) = smooth deriv i
  Xi3 <- Xdp3*0 
  Xi3[,1:(knot-1)+7] <- Xdp3[,1:(knot-1)+7] ## Xi%*%coef(MGAM) = smooth deriv i
  df <- Xi%*%coef(GAM)              ## ith smooth derivative 
  dfsd <- rowSums(Xi%*%GAM$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  df2 <- Xi2%*%coef(GAM)              ## ith smooth derivative 
  df2sd <- rowSums(Xi2%*%GAM$Vp*Xi2)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  df3 <- Xi3%*%coef(GAM)              ## ith smooth derivative 
  
#  res$sdt <- df
 # res$sdtup <- df+2*df.sd
  #res$sdtlow <- df-2*df.sd
  #res$time <- as.Date(xv, origin = "1899-12-30")
  # CI is approximate.
  
#  Xi <- Xdp*0
 # Xi[,1:(knot-1)+7] <- Xdp[,1:(knot-1)+7] ## Xi%*%coef(MGAM) = smooth deriv i
  #df <- Xi%*%coef(GAM)              ## ith smooth derivative
  
#  dfsd <- rowSums(Xi%*%GAM$Vc*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  res$sdt <- df
  res$sdtup <- df+zv*dfsd
  res$sdtlow <- df-zv*dfsd
  res$sdt2 <- df2
  res$sdt2up <- df2+2*df2sd
  res$sdt2low <- df2-2*df2sd
  res$sdt3 <- df3
  
# if plotting then this plots the data (red), the GAM with weekday terms (green), the GAM averaged over weekday (blue) 
  # and CI on this average GAM output. 
  if(plt==TRUE){
    mindate <- min(res$Date[1], min(as.Date(keydates$start)))
    maxdate <- max(res$Date[length(res$Date)], max(as.Date(keydates$end)))
    p1 <- ggplot(data=res)+
    geom_point(mapping=aes(x=Date, y=vals), color='red')+
    geom_line(mapping=aes(x=Date, y=fullgam), color='green')+
    geom_ribbon(mapping=aes(x=Date, ymax=predUNLup, ymin=predUNLlow), alpha=0.2, fill='blue')+
    geom_line(mapping=aes(x=Date, y=pred), color='blue')+
    ylab(paste(ifelse(type=='SIS', 'Outbreak', 'Case'), ' Count', sep=''))+xlab('Date')+
    xlim(c(mindate,maxdate))+# as.Date(c("2020-03-15", "2021-10-10")))+
    ylim(0,1.1*max(res$vals, res$predUNLup))+labs(title = paste(figtitle, sep=''))
  }

  # Now we set up calculation of reproduction number
  nv <- 10
  dt <- 1/nv
  I0 <- para[1]
  N <- para[2]
  gamma <- 1/para[3]
  Incid <- res$pred*exp(gamma*(res$time))
  Incidup <- res$predUNLup*exp(gamma*(res$time))
  Incidlow <- res$predUNLlow*exp(gamma*(res$time))
  # 
  ntv <- seq(nt[1], nt[length(nt)], dt)
  # predict over finer time grid than days
  xptrap=predict(GAM,newdata=data.frame(time=ntv,weekday="Friday"),type="lpmatrix")  
  # put this on response scale
  bootsimfine <- exp(xptrap %*% t(betasim))
  # but also retain on the daily timescale
  X0 <- Xp %*% t(betasim)
  bootsim <- exp(X0)
  # take quantiles for fast (but incorrect) CIs - these will be too wide
  qntsim <- apply(bootsimfine, 1, quantile, probs=c(0.025, 0.975))
  # this runs ok speed wise with all samples but taking a subset of the sample may be faster if debugging/testing
  sampboot <- sample(1:10000, 10000, replace=FALSE) #1:10000 #
  wdmean <- mean(exp(apply(betasim[sampboot, 2:7], 1, sum)/7))

  # use the trapesuim rule to approximate integrals [this is the slow bit if nv is large, or sample size in sampboot large]
    Fvalboot <- bootsim[, sampboot]*exp(gamma*(nt))
  Fvalgradboot <-  bootsimfine[,sampboot]*exp(gamma*(ntv)) #sum(coef(GAM)[2:7])/7+
  Fvalgradbootsum <- dt*apply(Fvalgradboot,2,cumsum)[1+nv*0:(length(nt)-1),]-
    dt/2*(Fvalgradboot[1+nv*0:(length(nt)-1),]+Fvalgradboot[1])
  Fintboot <- Fvalboot/(I0*exp(-sum(coef(GAM)[2:7])/7)+Fvalgradbootsum)/gamma
  qntsimfullint <- apply(Fintboot, 1, quantile, probs=c(0.025,0.5, 0.975))
  
  # redundant code but useful comparison of 'wrong' appraoch
  Fvalgrad <- exp(as.numeric(xptrap%*%coef(GAM))+sum(coef(GAM)[2:7])/7+gamma*(ntv))
  Fvalgradup <- qntsim[2,]*exp(sum(coef(GAM)[2:7])/7+gamma*(ntv))
  Fvalgradlow <- qntsim[1,]*exp(sum(coef(GAM)[2:7])/7+gamma*(ntv))
  Fvalgradsum = dt*(cumsum(Fvalgrad)) # -Fvalgrad[1]/2)
  Fvalgradsumup = dt*(cumsum(Fvalgradup)) # -Fvalgrad[1]/2)
  Fvalgradsumlow = dt*(cumsum(Fvalgradlow)) # -Fvalgrad[1]/2)
  Fint = I0+Fvalgradsum[1+nv*0:(length(Incid)-1)]-dt*(Fvalgrad[1+nv*0:(length(Incid)-1)]+Fvalgrad[1])/2
  Fintup = I0+Fvalgradsumup[1+nv*0:(length(Incid)-1)]-dt*(Fvalgradup[1+nv*0:(length(Incid)-1)]+Fvalgradup[1])/2
  Fintlow = I0+Fvalgradsumlow[1+nv*0:(length(Incid)-1)]-dt*(Fvalgradlow[1+nv*0:(length(Incid)-1)]+Fvalgradlow[1])/2
  
  if(type=='SEIR' || type=='SEIRWeek'){
    alpha <- 1/para[4]

    # now merge to estimator for RE
    alpha <- 1/para[4]
    res$REfulllow <- qntsimfullint[1,]*(1+res$sdtlow/alpha)
    res$REfullmed <- qntsimfullint[2,]*(1+res$sdt/alpha)
    res$REfullup <- qntsimfullint[3,]*(1+res$sdtup/alpha)
    res$RE <- Incid/Fint/gamma*(1+res$sdt/alpha)
    res$REup <- Incidup/Fintlow/gamma*(1+res$sdtup/alpha)
    res$RElow <- Incidlow/Fintup/gamma*(1+res$sdtlow/alpha)

    # use the trapesuim rule to approximate integrals [this is the slow bit if nv is large, or sample size in sampboot large]
    Gvalboot <- bootsim[, sampboot]
    Gvalgradboot <-  bootsimfine[,sampboot] #sum(coef(GAM)[2:7])/7+
    Gvalgradbootsum <- dt*apply(Gvalgradboot,2,cumsum)[1+nv*0:(length(nt)-1),]-
      dt/2*(Gvalgradboot[1+nv*0:(length(nt)-1),]+Gvalgradboot[1])
    Gintboot <- (I0*exp(-sum(coef(GAM)[2:7])/7)+Gvalgradbootsum)
 #   print(c(dim(Fintboot),dim(bootsim),dim(Gintboot)))
    RCqntsimfullint <- apply(Fintboot/(1-Gvalboot/N/alpha-Gintboot/N), 1, quantile, probs=c(0.025,0.5, 0.975))
    RCqntsimfullint <- ifelse(RCqntsimfullint<0, 0, RCqntsimfullint)
    # again redundant code which could be removed
    Gvalgrad <- exp(as.numeric(xptrap%*%coef(GAM))+sum(coef(GAM)[2:7])/7)
    Gvalgradup <- qntsim[2,]*exp(sum(coef(GAM)[2:7])/7)
    Gvalgradlow <- qntsim[1,]*exp(sum(coef(GAM)[2:7])/7)
    Gvalgradsum = dt*(cumsum(Gvalgrad)) 
    Gvalgradsumup = dt*(cumsum(Gvalgradup)) 
    Gvalgradsumlow = dt*(cumsum(Gvalgradlow))
    Gint = I0+Gvalgradsum[1+nv*0:(length(Incid)-1)]-dt*(Gvalgrad[1+nv*0:(length(Incid)-1)]+Gvalgrad[1])/2
    Gintup = I0+Gvalgradsumup[1+nv*0:(length(Incid)-1)]-dt*(Gvalgradup[1+nv*0:(length(Incid)-1)]+Gvalgradup[1])/2
    Gintlow = I0+Gvalgradsumlow[1+nv*0:(length(Incid)-1)]-dt*(Gvalgradlow[1+nv*0:(length(Incid)-1)]+Gvalgradlow[1])/2
    
    res$RC <- res$RE/(1-res$pred/alpha/N-Gint/N)
    res$RCup <- res$REup/(1-res$predUNLlow/alpha/N-Gintlow/N)
    res$RClow <- res$RElow/(1-res$predUNLup/alpha/N-Gintup/N)
    res$RCfulllow <- RCqntsimfullint[1,]*(1+res$sdtlow/alpha)
    res$RCfullmed <- RCqntsimfullint[2,]*(1+res$sdt/alpha)
    res$RCfullup <-RCqntsimfullint[3,]*(1+res$sdtup/alpha)
  }else if(type=='SIS'){
    
    res$REfulllow <- qntsimfullint[1,]
    res$REfullmed <- qntsimfullint[2,]
    res$REfullup <- qntsimfullint[3,]
    res$RE <- Incid/Fint/gamma
    res$REup <- Incidup/Fintlow/gamma
    res$RElow <- Incidlow/Fintup/gamma
    
    PrevEst <- Fintboot/(1-(I0*exp(-sum(coef(GAM)[2:7])/7)+Fvalgradbootsum)/N*exp(-gamma*res$time))
    RCqntsimfullint <- apply(PrevEst, 1, quantile, probs=c(0.025,0.5, 0.975))

    res$RCfulllow <- RCqntsimfullint[1,]
    res$RCfullmed <- RCqntsimfullint[2,]
    res$RCfullup <- RCqntsimfullint[3,]
    # redundant code
    res$RC <- res$RE/(1-Fint/N*exp(-gamma*res$time))
    res$RCup <- res$REup/(1-Fintlow/N*exp(-gamma*res$time))
    res$RClow <- res$RElow/(1-Fintup/N*exp(-gamma*res$time))
    
    #    rhotRC <- rhotRE/(1-Fint*exp(-nt*gamma)/N)
    
  }
  if(EEwindow>0){
#    print('Generate EpiEstim for comparison')
    ee_data=res[,c('Date','vals')]
    names(ee_data)=c('dates','I')
    
    window <- EEwindow
    start_t<-seq(2,length(res$t)-window+1)
    end_t<-start_t+window-1
    
    EEcomp=EpiEstim::estimate_R(ee_data,method="parametric_si",
                                config=EpiEstim::make_config(list(mean_si=para[3]+para[4],std_si=5.83,
                                                                  t_start =start_t,
                                                                  t_end =end_t)))
    # set up data frame to render plots - note the +para[3] to synchronise the R_E estimators
    datEpiEstRE <- data.frame(Date = EEcomp$dates[EEcomp$R$t_start+para[3]],  
                              Mean =EEcomp$R$'Mean(R)', Low = EEcomp$R$'Quantile.0.025(R)', Upp = EEcomp$R$'Quantile.0.975(R)')
    }
    if(plt==TRUE){
    mrc <- max(res$RCfullup[5:(length(res$RCfullup)-5)])
    if(EEwindow>0){
      if(simulatorR>0){
        p2 <- ggplot(data=res)+
          annotate("rect", xmin=as.Date(keydates$start),xmax=as.Date(keydates$end),
                   ymin=rep(0, length(keydates$start)),ymax=rep(mrc, length(keydates$start)),
                   alpha = .2, fill='lightgreen')+
          geom_hline(yintercept=1, color='white', lwd=2)+
          geom_hline(yintercept=simulatorR, color='white', lwd=2)+
          geom_vline(xintercept=as.Date(keydates$start[1]-(para[3]+para[4])), color='white', lwd=2)+
          geom_vline(xintercept=as.Date(keydates$start[1]+para[3]+para[4]), color='white', lwd=2)+
          geom_ribbon(mapping=aes(x=Date, ymax=REfullup, ymin=REfulllow), alpha=0.2, fill='red')+
          geom_ribbon(data=datEpiEstRE, mapping=aes(x=Date, ymax=Upp, ymin=Low), fill='purple', alpha=0.2)+
          geom_ribbon(mapping=aes(x=Date, ymax=RCfullup, ymin=RCfulllow), alpha=0.2, fill='black')+
          geom_line(mapping=aes(x=Date, y=REfullmed), color='red')+
          geom_line(mapping=aes(x=Date, y=RCfullmed), color='black')+
          geom_line(data=datEpiEstRE, mapping=aes(x=Date, y=Mean), color='purple')+
          ylab('Reproduction Number') + xlab('Date')+
          xlim(c(mindate,maxdate)) + ylim(0,mrc)
      }else{
        p2 <- ggplot(data=res)+
        annotate("rect", xmin=as.Date(keydates$start),xmax=as.Date(keydates$end),
                       ymin=rep(0, length(keydates$start)),ymax=rep(mrc, length(keydates$start)),
                      alpha = .2, fill='lightgreen')+
        geom_hline(yintercept=1, color='white', lwd=2)+
        geom_ribbon(mapping=aes(x=Date, ymax=REfullup, ymin=REfulllow), alpha=0.2, fill='red')+
          geom_ribbon(data=datEpiEstRE, mapping=aes(x=Date, ymax=Upp, ymin=Low), fill='purple', alpha=0.2)+
        geom_ribbon(mapping=aes(x=Date, ymax=RCfullup, ymin=RCfulllow), alpha=0.2, fill='black')+
          geom_line(mapping=aes(x=Date, y=REfullmed), color='red')+
        geom_line(mapping=aes(x=Date, y=RCfullmed), color='black')+
        geom_line(data=datEpiEstRE, mapping=aes(x=Date, y=Mean), color='purple')+
        ylab('Reproduction Number') + xlab('Date')+
        xlim(c(mindate,maxdate)) + ylim(0,mrc)
      }
    }else{
      if(simulatorR>0){
        p2 <- ggplot(data=res)+
          annotate("rect", xmin=as.Date(keydates$start),xmax=as.Date(keydates$end),
                   ymin=rep(0, length(keydates$start)),ymax=rep(mrc, length(keydates$start)),
                   alpha = .2, fill='lightgreen')+
          geom_hline(yintercept=1, color='white', lwd=2)+
          geom_hline(yintercept=simulatorR, color='white', lwd=2)+
          geom_vline(xintercept=as.Date(keydates$start[1]-para[3]), color='white', lwd=2)+
          geom_vline(xintercept=as.Date(keydates$start[1]+para[3]), color='white', lwd=2)+
          geom_ribbon(mapping=aes(x=Date, ymax=REfullup, ymin=REfulllow), alpha=0.2, fill='red')+
          geom_line(mapping=aes(x=Date, y=REfullmed), color='red')+
          geom_ribbon(mapping=aes(x=Date, ymax=RCfullup, ymin=RCfulllow), alpha=0.2, fill='black')+
          geom_line(mapping=aes(x=Date, y=RCfullmed), color='black')+
          ylab('Reproduction Number')+xlab('Date')+
          xlim(c(mindate,maxdate))+
          ylim(0,mrc)
      }else{
        p2 <- ggplot(data=res)+
          annotate("rect", xmin=as.Date(keydates$start),xmax=as.Date(keydates$end),
                   ymin=rep(0, length(keydates$start)),ymax=rep(mrc, length(keydates$start)),
                   alpha = .2, fill='lightgreen')+
      geom_hline(yintercept=1, color='white', lwd=2)+
      geom_ribbon(mapping=aes(x=Date, ymax=REfullup, ymin=REfulllow), alpha=0.2, fill='red')+
      geom_line(mapping=aes(x=Date, y=REfullmed), color='red')+
      geom_ribbon(mapping=aes(x=Date, ymax=RCfullup, ymin=RCfulllow), alpha=0.2, fill='black')+
      geom_line(mapping=aes(x=Date, y=RCfullmed), color='black')+
      ylab('Reproduction Number')+xlab('Date')+
          xlim(c(mindate,maxdate))+
      ylim(0,mrc)
      }
    }
        if(pltgr==TRUE){
    if(res$sdtlow[npts]>0){
      resstring <- paste('Increasing recent trend', sep='')
    }else if(res$sdtup[npts]<0){
      resstring <- paste('Decreasing recent trend', sep='')
    }else{
      if(res$sdt[npts]>0){
        resstring <- paste('Plateauing recent trend, may be increasing with probability ',round(1-pnorm(0, res$sdt[npts],dfsd[npts]), 3), sep='')
      }else{
        resstring <- paste('Plateauing recent trend, may be decreasing with probability ',round(pnorm(0, res$sdt[npts],dfsd[npts]),3), sep='')
      }
    }
    p3 <- ggplot(data=res, aes(x=time, y=sdt))+
      geom_ribbon(aes(ymin=sdtlow, ymax=sdtup), fill='lightblue')+#+ 
      geom_line(color='darkblue') +
      labs(x='Date', title='Derivative of spline arising from GAM', caption=resstring)+
      scale_y_continuous(
        name = expression("Instantaneous growth rate"), 
        sec.axis = sec_axis(~., name = "Doubling Time", 
                            breaks = c(-0.1,-0.05,-0.033,-0.0247, 0,0.0247, 0.033, 0.0495, 0.099, 0.173, 0.3466), 
                            labels=c(-7,-14,-21,-28, 'Infinite',28,21, 14, 7,4,2)), 
        limits = c(min(res$sdtlow), max(res$sdtup)))   # theme(
    print(p1+ p2+p3+plot_layout(ncol=1))
    }else{
      print(p1+p2+plot_layout(ncol=1))
    }
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  
      res
}




# simulates stochastic version of SEIR model
seir_tau_tv=function(dt=0.15,MaxTime=1000000,N0=10000,E0=100,I0=0,repro=2*log(2),sig=1/1.3,gam=1/5,itp=0.5,frac=0.5,plt1=T,plt2=T)
{
  t=dt*(0:MaxTime)
  S = N0-E0-I0
  E = E0
  I = I0
  R = 0
  Ihat=0
  
  for (i in 2:(MaxTime+1)){ 
    Rt=ifelse(i>MaxTime*itp,frac,1)*repro*gam*S[i-1]*I[i-1]/N0   #transmission rate
    Ri=sig*E[i-1]                  #becoming infectious rate
    Rr=gam*I[i-1]          #recover rate
    
    dmt=ifelse(Rt>0,rpois(1,Rt*dt),0)
    dmi=ifelse(Ri>0,rpois(1,Ri*dt),0)
    dmr=ifelse(Rr>0,rpois(1,Rr*dt),0)
    
    if(dmt>S[i-1]) dmt=S[i-1]
    if(dmi>E[i-1]) dmi=E[i-1]
    if(dmr>I[i-1]) dmr=I[i-1]
    
    S[i]=S[i-1]-dmt
    E[i]=E[i-1]+dmt-dmi
    I[i]=I[i-1]+dmi-dmr
    R[i]=R[i-1]+dmr
    Ihat=append(Ihat,dmi)
  }
  
  if(plt1==T){
    plot(S~t,type="l",ylim=c(0,N0),ylab="Number",xlab="Time",main="SEIR")
    lines(E~t,col="gold2")
    lines(I~t,col="firebrick")
    lines(R~t,col="blue3")
    legend("right", legend=c("S","E","I","R"),col=c(1,"gold2","firebrick","blue3"), lwd=c(1,1), cex=0.5)
  }
  
  if(plt2==T){
    plot(I~t,type="l",ylim=c(0,max(I)),ylab="Number",col="firebrick",xlab="Time"#,main="SEIR"
    )
    lines(E~t,col="gold2")
    legend("right", legend=c("E","I"),col=c("gold2","firebrick"), lwd=c(1,1), , cex=1.5)
  }
  data<- data.frame(time=t,S=S,E=E,I=I,R=R,Ihat=Ihat)
  data
}

#####tau leap with intervention
sis_tau_tv=function(dt=0.1,MaxTime=1000000,N0=10000,I0=10,trans=2,recrate=1/20,itp=0.5,frac=0.25,plt=T)
{
  t=dt*(0:MaxTime)
  S=N0-I0
  I=I0
  Ihat=0
  
  for (i in 2:(MaxTime+1)) 
  {
    R1=ifelse(i>MaxTime*itp,frac,1)*trans*(N0-I[i-1])/N0  #transmission rate
    R2=recrate*I[i-1]  # recovery rate
    
    dmt=rpois(1,R2*dt*R1)
    dmr=rpois(1,R2*dt)
    
    if(dmt>N0-I[i-1]) dmt=N0-I[i-1]
    if(dmr>I[i-1]) dmr=I[i-1]
    
#    S[i]=S[i-1]-dmt+dmr
    I[i]=I[i-1]+dmt-dmr
    Ihat=append(Ihat,dmt)
  }
  S = N0-I
  
  if(plt==T){
    #par(mfrow=c(1,1))
    plot(S~t,type="l",ylim=c(0,N0),ylab="Number",xlab="Time",main="SIS")
    lines(I~t,col="firebrick")
    legend("topright", legend=c("S","I"),col=c(1,"firebrick"), lwd=c(1,1), cex=0.5)
  }
  data<- data.frame(time=t,S=S,I=I,Ihat=Ihat)
  data
}

# Function that wraps together simulator of outbreak and R_C estimation for SIS
SimulatorTVRSIS <- function(epipara = c(14500, 26, 2, 145, 2), contpara=c(0.5,0.5), plotgraph=FALSE){
  repTS <- 10
  offset <- 18353
  duration <- 365*epipara[3]
  startdate <- as.Date(offset) # 
  interventdate <- as.Date(floor(offset+duration*contpara[1]))
  enddate <-  as.Date(offset+duration)
  
  Sim <- sis_tau_tv(MaxTime=duration*repTS,dt=1/repTS, trans=epipara[5], recrate=1/epipara[2], 
                    N0=epipara[1], I0=epipara[4], itp=contpara[1], frac=contpara[2], plt = F )
  Simdata <- cbind(offset-1+1:(floor(dim(Sim)[1])/repTS), diff(cumsum(Sim$Ihat)[1+repTS*(0:floor(dim(Sim)[1]/repTS))]))
  tvrout <- gamTVRpred(data=Simdata, type='SIS', para=c(epipara[4], epipara[1], epipara[2], 0), EEwindow=0, 
                     keydates = data.frame(start=interventdate, end="2022-04-01"), simulatorR = epipara[5],  plt=plotgraph)

  output <- cbind(tvrout$Date, tvrout$REfullmed,tvrout$REfullup,tvrout$REfulllow,tvrout$RCfullmed,tvrout$RCfullup,tvrout$RCfulllow)
  output
}

# Function that wraps together simulator of outbreak and R_C estimation for SEIR
SimulatorTVRSEIR <- function(epipara = c(10^5, log(4)*4/3, 5, 3, 10^2, 1), contpara=c(0.5,0.5), plotgraph=FALSE){
repTS <- 10
pop <- epipara[1]
R0 <- epipara[2] 
InfP <- epipara[3]
IncP <- epipara[4]
InitI <- epipara[5]
offset <- 18353
duration <- 365*epipara[6]
startdate <- as.Date(offset) # 
enddate <-  as.Date(offset+duration)

  interventdate <- as.Date(offset+contpara[1]*duration)
  FadeOut <- 0
  # lower initial conditions or transmission - be careful of fade out so loop to catch outbreak that doesnt fade out
  while(FadeOut == 0){
    Sim <- seir_tau_tv(MaxTime=duration*repTS,dt=1/repTS, repro=R0,sig=1/IncP, gam= 1/InfP, N0=pop, I0=InitI,E0=0, 
                       itp=contpara[1], frac=contpara[2], plt1 = F, plt2 = F )
    Simdata <- cbind(offset-1+1:(floor(dim(Sim)[1])/repTS), diff(cumsum(Sim$Ihat)[1+repTS*(0:floor(dim(Sim)[1]/repTS))]))
    if(Sim$R[dim(Sim)[1]] > 2*10^3 && sum(Simdata[100:365,2]<1)>0){
      
      FadeOut <- 1}
  }
  #
  Coreperiod <- (which(Simdata[100:365,2]<1)[1]+100)
  tvrout <- gamTVRpred(data=Simdata[1:Coreperiod,], para=c(InitI,pop, InfP,IncP), 
                       keydates = data.frame(start=interventdate, end=min(offset+Coreperiod, enddate)), simulatorR=R0,plt=plotgraph)#,
#                       figtitle = paste('Simulation with 50% reduction in R0=',round(R0,2),' on ', interventdate, sep=''))
  tid <- which(tvrout$Date<interventdate)

    output <- cbind(tvrout$Date, tvrout$REfullmed,tvrout$REfullup,tvrout$REfulllow,tvrout$RCfullmed,tvrout$RCfullup,tvrout$RCfulllow)
  
    output
}

# Function that wraps together simulator of outbreak and R_C estimation for SEIR
SimulatorTVRSEIRUnderReporting <- function(epipara = c(10^5, log(4)*4/3, 5, 3, 10^2, 1), contpara=c(0.5,0.5), plotgraph=FALSE){
  repTS <- 10
  
  SummaryOutput <- data.frame(dates = seq(18353, 18353+365/2,1)) 
  pop <- epipara[1]
  R0 <- epipara[2] 
  InfP <- epipara[3]
  IncP <- epipara[4]
  InitI <- epipara[5]
  offset <- 18353
  duration <- 365*epipara[6]
  startdate <- as.Date(offset) # 
  enddate <-  as.Date(offset+duration)
  Coreperiod <- 365 #(which(Simdata[100:365,2]<1)[1]+100)
  
  interventdate <- as.Date(offset+contpara[1]*duration)
    Sim <- seir_tau_tv(MaxTime=duration*repTS,dt=1/repTS, repro=R0,sig=1/IncP, gam= 1/InfP, N0=pop, I0=InitI,E0=0, 
                       itp=contpara[1], frac=contpara[2], plt1 = F, plt2 = F )
    Simdata <- cbind(offset-1+1:(floor(dim(Sim)[1])/repTS), diff(cumsum(Sim$Ihat)[1+repTS*(0:floor(dim(Sim)[1]/repTS))]))
  #
  tvrout <- gamTVRpred(data=Simdata[1:Coreperiod,], para=c(InitI,pop, InfP,IncP), 
                       keydates = data.frame(start=interventdate, end=min(offset+Coreperiod, enddate)), simulatorR=R0,plt=plotgraph,
                         figtitle = paste('Simulation with R0=',round(R0,2), sep=''))
  output <- cbind(tvrout$Date, tvrout$REfullmed,tvrout$REfullup,tvrout$REfulllow,tvrout$RCfullmed,tvrout$RCfullup,tvrout$RCfulllow)
  
  output <- data.frame(output)
  colnames(output) <- c('dates', paste('RE',sep=''), paste('REup',sep=''), paste('RElow',sep=''), 
                        paste('RC', sep=''), paste('RCup', sep=''), paste('RClow',sep=''))
  SummaryOutput <- Reduce(function(x, y) merge(x, y, all=TRUE), list(SummaryOutput, output))
  #  for(i in seq(1,9,2)){
 # SimdataUR <- cbind(offset-1+1:(floor(dim(Sim)[1])/repTS),rbeta(n=length(Simdata[,2]), shape1 = i, shape2 = 2)*Simdata[,2])
  #tvrout <- gamTVRpred(data=SimdataUR[1:Coreperiod,], para=c(InitI,pop, InfP,IncP), 
   #                    keydates = data.frame(start=interventdate, end=min(offset+Coreperiod, enddate)), simulatorR=R0,plt=plotgraph,
    #                   figtitle = paste('Simulation with R0=',round(R0,2),' and mean under-reporting of ', round(i/(i+2),2), sep=''))#,
  #}
  for(i in seq(1,4,1)){
    SimdataUR <- cbind(offset-1+1:(floor(dim(Sim)[1])/repTS),rbeta(n=length(Simdata[,2]), shape1 = 10^(i-2), shape2 = 2*10^(i-2))*Simdata[,2])
    tvrout <- gamTVRpred(data=SimdataUR[1:Coreperiod,], para=c(InitI,pop, InfP,IncP), 
                         keydates = data.frame(start=interventdate, end=min(offset+Coreperiod, enddate)), simulatorR=R0,plt=plotgraph,
                         figtitle = paste('Simulation with R0=',round(R0,2),' and mean under-reporting of 1/3 and overdisperision ', 10^(i-2), sep=''))#,
    output <- cbind(tvrout$Date, tvrout$REfullmed,tvrout$REfullup,tvrout$REfulllow,tvrout$RCfullmed,tvrout$RCfullup,tvrout$RCfulllow)
    
    output <- data.frame(output)
    colnames(output) <- c('dates', paste('RE', i,sep=''), paste('REup', i,sep=''), paste('RElow', i,sep=''), 
                          paste('RC', i,sep=''), paste('RCup', i,sep=''), paste('RClow', i,sep=''))
    SummaryOutput <- Reduce(function(x, y) merge(x, y, all=TRUE), list(SummaryOutput, output))
  }
  
#  output <- cbind(tvrout$Date, tvrout$REfullmed,tvrout$REfullup,tvrout$REfulllow,tvrout$RCfullmed,tvrout$RCfullup,tvrout$RCfulllow)
  
  SummaryOutput
}
