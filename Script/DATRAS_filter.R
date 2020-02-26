#load packages----
library(nonlinearTseries)
library(purrr)
library(plyr)
library(tidyr)
library(pbapply)
library(dplyr)
library(data.table)
library(ggplot2)
library(magick)
library(here) 
library(magrittr)
library(png)
library(grid)
library(jpeg)
library(gridExtra)
library(rgdal)
library(sf)
library(gratia)
library(mgcv)
library(itsadug)
library(INDperform)
library(gtools)
library(imager)
library(mgcViz)
library(pracma)
library(cdfquantreg)
library(remotes)
library(tseries)
library(NlinTS)
library(astsa)
library(TTR)
library(fpp)
library(Rsenal)
library(ggmap)
library(rbin)

st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

#function for selecting specific column to omit NA----
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

#deriv function----
tmp <- tempfile()
Deriv <- function(mod, n = 200, eps = 1e-7, newdata) {
  #if(isTRUE(inherits(mod,``list``)))
  mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  # number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  return(lD)
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term))
    term <- term.labs
  Term <- match(term, term.labs)
  ##term <- term[match(term, term.labs)]
  if(any(miss <- is.na(Term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  ## if(is.na(term))
  ##     stop("'term' not a valid model term.")
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- length(object$gamModel$y) - sum(object$gamModel$edf)
  tVal <- qt(1 - (alpha/2), residual.df)
  ## tVal <- qt(1 - (alpha/2), object$gamModel$df.residual)
  for(i in seq_along(term)) {
    upr <- object[[term[i]]]$deriv + tVal * object[[term[i]]]$se.deriv
    lwr <- object[[term[i]]]$deriv - tVal * object[[term[i]]]$se.deriv
    res[[term[i]]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

plot.Deriv <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term, eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term))
    term <- term.labs
  Term <- match(term, term.labs)
  if(any(miss <- is.na(Term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(is.na(Term)))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  ## tVal <- qt(1 - (alpha/2), x$gamModel$df.residual)
  residual.df <- length(x$gamModel$y) - sum(x$gamModel$edf)
  tVal <- qt(1 - (alpha/2), residual.df)
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")[Term]
    names(xlab) <- xlab
  }
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  CI <- confint(x, term = term, alpha = alpha)
  for(i in seq_along(term)) {
    ## for(i in seq_len(l)) {
    upr <- CI[[term[i]]]$upper
    lwr <- CI[[term[i]]]$lower
    ylim <- range(upr, lwr)
    plot(x$eval[,term[i]], x[[term[i]]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[term[i]], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,term[i]], rev(x$eval[,term[i]])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,term[i]], upr, lty = "dashed")
      lines(x$eval[,term[i]], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 1)
      S <- signifD(x[[term[i]]]$deriv, x[[term[i]]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,term[i]], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,term[i]], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
}

#function to generate random values from multivariate normal dist.----
rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}







#set path----
path=("C:/Users/mrademaker/Documents/Research projects/STCNWS/DATRAS/cpue_length_hour/")

#file with species info----
species_info=read.table(paste(path,"Data/spec_info.txt",sep=""),sep="\t",header=TRUE)
species_info$file_name =gsub(" ","_",species_info$Scientific.name)

#open seascapes map and convert to wgs84 long,lat----
ssc_map=readOGR("C:/Users/mrademaker/Documents/Research projects/STCNWS/Seascapes/To_Mark_ArcGIS.",layer="seascapes")
ssc_map_wgs84 <- spTransform(ssc_map, CRS("+proj=longlat +datum=WGS84")) 

#filter and adjust data_sets----
for (i in 3:3){#1:nrow(species_info)){
  file_name=paste(path,"Data/",species_info$file_name[i],sep="")
  print(file_name)
  #Read in data
  data=read.csv(paste(file_name,".csv",sep=""))
  
  #Subset with same fishing gear
  data=subset(data, Gear == 'GOV')
  
  #Subset only during day
  #data=subset(data, DayNight == "D")
  
  #Bin depth in 20m classes
  #print(max(data$Depth))
  #data$Depth=abs(data$Depth)
  #data$Depth_bin=factor(cut(data$Depth,breaks=c(0,50,100,150,200,250,300,350,400)))
  # 
  #Length in cm
  data$LngtClass_cm=data$LngtClass/10

  #Add info on weight conversion formula
  data$WeightConv=species_info$LW.conversion[i]

  #Calculate weight for corresponding Length Class
  data$Weight=species_info$a[i]*(data$LngtClass_cm^species_info$b[i])

  # Multiply weight by CPUE number per hour for indication of total weight(biomass) per length class
  data$Total_wgt=data$CPUE_number_per_hour*data$Weight

  #Determine location (Seascape)
  coordinates(data)= ~ ShootLong + ShootLat
  proj4string(data) <- proj4string(ssc_map_wgs84)
  data_in=as.data.frame(over(data, ssc_map_wgs84))
  data_ssc=as.data.frame(cbind(as.data.frame(data),data_in))

  # # # check on map
  # # shapefile_df <- fortify(ssc_map_wgs84)
  # #
  # # map <- ggplot() +
  # #   geom_path(data = shapefile_df,
  #             aes(x = long, y = lat, group = group),
  #             color = 'black', size = 1)+
  #   ggtitle("Seascapes plot and NS-IBTS survey data (Limanda limanda)")
  # print(map)
  #
  # p2 = map+geom_point(data=as.data.frame(data),aes(x=ShootLong,y=ShootLat),color="black",fill="red",size=1,shape=21)
  # print(p2)
  #
  #save filtered file
  file_name=sprintf("C:/Users/mrademaker/Documents/Research projects/STCNWS/DATRAS/cpue_length_hour/Data/Filtered_%s.csv",species_info$file_name[i])
  write.csv(data_ssc,file_name)
  #
}


################################################################################
####GAM Anscombe transformed (sqrt(x+3/8)) and scaled by max (BETA errror)   ###-----
################################################################################
for (j in c(1,13)){
  data=read.csv(paste(path,sprintf("Data/Filtered_%s.csv",species_info$file_name[j]),sep=""))
  names(data)[names(data) == "id"] <- "Seascapenr"
  print(species_info$file_name[j])
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  ssc_data=completeFun(data,"Seascapenr")
  ssc_data=subset(ssc_data,Quarter==1)
  data_agg=ssc_data
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species+Ship+DayNight, data=ssc_data, sum)
  #print(nrow(data_agg))
  #positive_data_agg=subset(data_agg, Total_wgt!= 0)
  #print(nrow(positive_data_agg))
  
  # # # # #set outlier weight limit (99th percentile)
  #   nf_p=quantile(data_agg$Total_wgt,0.95)
  #   if (nf_p > 0){
  #    #subset of data excluding outlier weights
  #     data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
  #   }
  # # 
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000

  #separate subset per seascape, anscombe transform, max scale it and finally rbind back into single dataset----
  #1----
  data_agg1=subset(data_agg,Seascapenr==1)
  data_agg1$ansc_biomass= 2*nthroot(data_agg1$biomass_kg+(3/8),2)
  data_agg1$scaled_ansc_biomass=scaleTR(data_agg1$ansc_biomass)
  data_agg1$scaled_ansc_biomass[is.na(data_agg1$scaled_ansc_biomass)] <- 0.0001

  #2-----
  data_agg2=subset(data_agg,Seascapenr==2)
  data_agg2$ansc_biomass= 2*nthroot(data_agg2$biomass_kg+(3/8),2)
  data_agg2$scaled_ansc_biomass=scaleTR(data_agg2$ansc_biomass)
  data_agg2$scaled_ansc_biomass[is.na(data_agg2$scaled_ansc_biomass)] <- 0.0001

  #3----
  data_agg3=subset(data_agg,Seascapenr==3)
  data_agg3$ansc_biomass= 2*nthroot(data_agg3$biomass_kg+(3/8),2)
  data_agg3$scaled_ansc_biomass=scaleTR(data_agg3$ansc_biomass)
  data_agg3$scaled_ansc_biomass[is.na(data_agg3$scaled_ansc_biomass)] <- 0.0001

  #4----
  data_agg4=subset(data_agg,Seascapenr==4)
  data_agg4$ansc_biomass= 2*nthroot(data_agg4$biomass_kg+(3/8),2)
  data_agg4$scaled_ansc_biomass=scaleTR(data_agg4$ansc_biomass)
  data_agg4$scaled_ansc_biomass[is.na(data_agg4$scaled_ansc_biomass)] <- 0.0001

  #5----
  data_agg5=subset(data_agg,Seascapenr==5)
  data_agg5$ansc_biomass= 2*nthroot(data_agg5$biomass_kg+(3/8),2)
  data_agg5$scaled_ansc_biomass=scaleTR(data_agg5$ansc_biomass)
  data_agg5$scaled_ansc_biomass[is.na(data_agg5$scaled_ansc_biomass)] <- 0.0001


  #6----
  data_agg6=subset(data_agg,Seascapenr==6)
  data_agg6$ansc_biomass= 2*nthroot(data_agg6$biomass_kg+(3/8),2)
  data_agg6$scaled_ansc_biomass=scaleTR(data_agg6$ansc_biomass)
  data_agg6$scaled_ansc_biomass[is.na(data_agg6$scaled_ansc_biomass)] <- 0.0001

  #7----
  data_agg7=subset(data_agg,Seascapenr==7)
  data_agg7$ansc_biomass= 2*nthroot(data_agg7$biomass_kg+(3/8),2)
  data_agg7$scaled_ansc_biomass=scaleTR(data_agg7$ansc_biomass)
  data_agg7$scaled_ansc_biomass[is.na(data_agg7$scaled_ansc_biomass)] <- 0.0001

  #8----
  data_agg8=subset(data_agg,Seascapenr==8)
  data_agg8$ansc_biomass= 2*nthroot(data_agg8$biomass_kg+(3/8),2)
  data_agg8$scaled_ansc_biomass=scaleTR(data_agg8$ansc_biomass)
  data_agg8$scaled_ansc_biomass[is.na(data_agg8$scaled_ansc_biomass)] <- 0.0001

  #9----
  data_agg9=subset(data_agg,Seascapenr==9)
  data_agg9$ansc_biomass= 2*nthroot(data_agg9$biomass_kg+(3/8),2)
  data_agg9$scaled_ansc_biomass=scaleTR(data_agg9$ansc_biomass)
  data_agg9$scaled_ansc_biomass[is.na(data_agg9$scaled_ansc_biomass)] <- 0.0001

  #10----
  data_agg10=subset(data_agg,Seascapenr==10)
  data_agg10$ansc_biomass= 2*nthroot(data_agg10$biomass_kg+(3/8),2)
  data_agg10$scaled_ansc_biomass=scaleTR(data_agg10$ansc_biomass)
  data_agg10$scaled_ansc_biomass[is.na(data_agg10$scaled_ansc_biomass)] <- 0.0001


  #rbind----
  data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9,data_agg10)

  ####################################################
  #GAM model with separate smoother per seascape
  ###################################################
  print("gam")
  data_agg_ansc$Seascapenr=as.factor(data_agg_ansc$Seascapenr) #Set as factor
  data_agg_ansc = start_event(data_agg_ansc, column="Year",event="Seascapenr") #Sep time series per factor level


  mod=bam(scaled_ansc_biomass~Seascapenr+s(Year, by=Seascapenr)+s(Ship,bs="re")+s(DayNight,bs="re"),data=data_agg_ansc,method="REML",family=betar)
#   #s(Ship,bs="re"),+s(Ship,Year,bs="re")
#   print("gam fitted")
#   print(AIC(mod))
#   par(mfrow=c(2,2))
#   gam.check(mod)
#   #observed_fitted_plot(mod)
#   #par(mfrow=c(2,2))
# 
#   dev.copy(png,file=paste(path,sprintf("Seascapes/Q1/%s/GAM/GAM_ANSC_diagnostics.png",species_info$file_name[j]),sep=""))
#   dev.off()
#   
#   
#   # ################################################################
#   # # Plot predicted smooths with all data                         #----
#   # ################################################################
# 
  # newd=expand.grid("Year"=seq(1977,2019,length = 126),
  #            "Seascapenr" = c(1,2,3,4,5,6,7,8,9,10),
  #            "Ship" = "WAH3", "DayNight" ="D")
  # 
  #  pred <- predict(mod, newd, se.fit = TRUE,type="response",exclude=c("s(Ship)","s(DayNight)"))
  #  se.fit <- pred$se.fit
  #  pred <- transform(cbind(data.frame(pred), newd),
  #                    uprP = fit + (2 * se.fit),
  #                    lwrP = fit - (2 * se.fit))
  # #
  # #separate plots per area
  # S1=pred[which(pred$Seascapenr==1),]
  # S2=pred[which(pred$Seascapenr==2),]
  # S3=pred[which(pred$Seascapenr==3),]
  # S4=pred[which(pred$Seascapenr==4),]
  # S5=pred[which(pred$Seascapenr==5),]
  # S6=pred[which(pred$Seascapenr==6),]
  # S7=pred[which(pred$Seascapenr==7),]
  # S8=pred[which(pred$Seascapenr==8),]
  # S9=pred[which(pred$Seascapenr==9),]
  # S10=pred[which(pred$Seascapenr==10),]
  # 
   img=readPNG(source=paste(path,sprintf("images/spec_%s.png",species_info$file_name[j]),sep=""))
   g <- rasterGrob(img, interpolate=TRUE)
  # # # #
  # # # # ### plot smooths with all data -----
  # s1plot=ggplot(data=S1,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg1)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S1")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # #Add species image
  # s1plot=s1plot+annotation_custom(g,xmin = min(data_agg1$Year), xmax = min(data_agg1$Year)+9, ymin = max(data_agg1$scaled_ansc_biomass)-0.15*(max(data_agg1$scaled_ansc_biomass)), ymax = max(data_agg1$scaled_ansc_biomass))
  # 
  # s2plot=ggplot(data=S2,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S2")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # #Add species image
  # #s2plot=s2plot+annotation_custom(g,xmin = min(data_agg2$Year), xmax = min(data_agg2$Year)+8, ymin = max(data_agg2$scaled_ansc_biomass)-0.15*(max(data_agg2$scaled_ansc_biomass)), ymax = max(data_agg2$scaled_ansc_biomass))
  # 
  # 
  # s3plot=ggplot(data=S3,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg3)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S3")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # #Add species image
  # #s3plot=s3plot+annotation_custom(g,xmin = min(data_agg3$Year), xmax = min(data_agg3$Year)+8, ymin = max(data_agg3$scaled_ansc_biomass)-0.15*(max(data_agg3$scaled_ansc_biomass)), ymax = max(data_agg3$scaled_ansc_biomass))
  # 
  # s4plot=ggplot(data=S4,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg4)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S4")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # #Add species image
  # #s4plot=s4plot+annotation_custom(g,xmin = min(data_agg4$Year), xmax = min(data_agg4$Year)+8, ymin = max(data_agg4$scaled_ansc_biomass)-0.15*(max(data_agg4$scaled_ansc_biomass)), ymax = max(data_agg4$scaled_ansc_biomass))
  # 
  # 
  # s5plot=ggplot(data=S5,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg5)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S5")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # #Add species image
  # #s5plot=s5plot+annotation_custom(g,xmin = min(data_agg5$Year), xmax = min(data_agg5$Year)+8, ymin = max(data_agg5$scaled_ansc_biomass)-0.15*(max(data_agg5$scaled_ansc_biomass)), ymax = max(data_agg5$scaled_ansc_biomass))
  # 
  # 
  # s6plot=ggplot(data=S6,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg6)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S6")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # #Add species image
  # #s6plot=s6plot+annotation_custom(g,xmin = min(data_agg6$Year), xmax = min(data_agg6$Year)+8, ymin = max(data_agg6$scaled_ansc_biomass)-0.15*(max(data_agg6$scaled_ansc_biomass)), ymax = max(data_agg6$scaled_ansc_biomass))
  # 
  # 
  # s7plot=ggplot(data=S7,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg7)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S7")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # #Add species image
  # #s7plot=s7plot+annotation_custom(g,xmin = min(data_agg7$Year), xmax = min(data_agg7$Year)+8, ymin = max(data_agg7$scaled_ansc_biomass)-0.15*(max(data_agg7$scaled_ansc_biomass)), ymax = max(data_agg7$scaled_ansc_biomass))
  # 
  # s8plot=ggplot(data=S8,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg8)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S8")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # #Add species image
  # #s8plot=s8plot+annotation_custom(g,xmin = min(data_agg8$Year), xmax = min(data_agg8$Year)+8, ymin = max(data_agg8$scaled_ansc_biomass)-0.15*(max(data_agg8$scaled_ansc_biomass)), ymax = max(data_agg8$scaled_ansc_biomass))
  # 
  # 
  # s9plot=ggplot(data=S9,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg9)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S9")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # #Add species image
  # #s9plot=s9plot+annotation_custom(g,xmin = min(data_agg9$Year), xmax = min(data_agg9$Year)+8, ymin = max(data_agg9$scaled_ansc_biomass)-0.15*(max(data_agg9$scaled_ansc_biomass)), ymax = max(data_agg9$scaled_ansc_biomass))
  # 
  # s10plot=ggplot(data=S10,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg10)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S10")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # #Add species image
  # #s10plot=s10plot+annotation_custom(g,xmin = min(data_agg10$Year), xmax = min(data_agg10$Year)+8, ymin = max(data_agg10$scaled_ansc_biomass)-0.15*(max(data_agg10$scaled_ansc_biomass)), ymax = max(data_agg10$scaled_ansc_biomass))
  # 
  # 
  # cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,s10plot,ncol=4,nrow=3)
  # ggsave(cplot,file=paste(path,sprintf("Seascapes/Q1/%s/GAM/GAM_ANSC_data_smooths.png",species_info$file_name[j]),sep=""),height=8, width=10,dpi = 600)#

# 
  # ################################################################
  # # Plot predicted smooths at smooth scale and simultaneous CI   #----
  # ################################################################
  #
  newd=expand.grid("Year"=seq(1977,2019,length = 126),
                   "Seascapenr" = c(1,2,3,4,5,6,7,8,9,10),
                   "Ship" = "WAH3", "DayNight"="D")

  pred <- predict(mod, newd, se.fit = TRUE, exclude=c("s(Ship)","s(DayNight)"))
  se.fit <- pred$se.fit
  pred <- transform(cbind(data.frame(pred), newd),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit))

  #################
  Vb <- vcov(mod)

  pred <- predict(mod, newd, se.fit = TRUE, exclude=c("s(Ship)","s(DayNight)"))
  se.fit <- pred$se.fit

  set.seed(42)
  N <- 10000

  BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  #
  # #xp matrix where basisfunctions of the model have been evaluated at 400 time point values per area
  Cg <- predict(mod, newd, type = "lpmatrix", exclude = "s(Ship)")
  simDev <- Cg %*% t(BUdiff)

  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  masd <- apply(absDev, 2L, max)
  crit <- quantile(masd, prob = 0.95, type = 8)
  pred <- transform(cbind(data.frame(pred), newd),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit),
                    uprS = fit + (crit * se.fit),
                    lwrS = fit - (crit * se.fit),
                    ci_width = (fit+(crit*se.fit))-(fit-(crit*se.fit)))
  write.csv(pred,file = paste(path,sprintf("Seascapes/Q1/%s/GAM/fitted_trends_ci.csv",species_info$file_name[j]),sep=""))
  #
  #
  # # #separate plots per area
  S1=pred[which(pred$Seascapenr==1),]
  S2=pred[which(pred$Seascapenr==2),]
  S3=pred[which(pred$Seascapenr==3),]
  S4=pred[which(pred$Seascapenr==4),]
  S5=pred[which(pred$Seascapenr==5),]
  S6=pred[which(pred$Seascapenr==6),]
  S7=pred[which(pred$Seascapenr==7),]
  S8=pred[which(pred$Seascapenr==8),]
  S9=pred[which(pred$Seascapenr==9),]
  S10=pred[which(pred$Seascapenr==10),]

  ### plot smooths with all data -----
  s1plot=ggplot(data=S1,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S1")+
    labs(x="Year",y="ANSC fish biomass trend")
  s1plot=s1plot+annotation_custom(g,xmin = min(data_agg1$Year), xmax = min(data_agg1$Year)+9, ymin = max(S1$fit)-0.2*(max(S1$fit)), ymax = max(S1$fit))


  s2plot=ggplot(data=S2,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S2")+
    labs(x="Year",y="ANSC fish biomass trend")


  s3plot=ggplot(data=S3,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S3")+
    labs(x="Year",y="ANSC fish biomass trend")


  s4plot=ggplot(data=S4,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S4")+
    labs(x="Year",y="ANSC fish biomass trend")

  s5plot=ggplot(data=S5,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S5")+
    labs(x="Year",y="ANSC fish biomass trend")

  s6plot=ggplot(data=S6,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S6")+
    labs(x="Year",y="ANSC fish biomass trend")

  s7plot=ggplot(data=S7,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S7")+
    labs(x="Year",y="ANSC fish biomass trend")
  s8plot=ggplot(data=S8,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S8")+
    labs(x="Year",y="ANSC fish biomass trend")

  s9plot=ggplot(data=S9,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S9")+
    labs(x="Year",y="ANSC fish biomass trend")

  s10plot=ggplot(data=S10,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S10")+
    labs(x="Year",y="ANSC fish biomass trend")

  cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,s10plot,ncol=4,nrow=3)
  ggsave(cplot,file=paste(path,sprintf("Seascapes/Q1/%s/GAM/GAM_ANSC_smooths.png",species_info$file_name[j]),sep=""),height=8, width=10,dpi = 600)#
  #
  # #save residual autocorrelation plot
  # par(mfrow=c(1,2))
  # acf(residuals(mod),main="ACF(mod)")
  # pacf(residuals(mod),main="partial ACF(mod)")
  # dev.copy(png,file=paste(path,sprintf("Seascapes/Q1/%s/GAM/GAM_ANSC_autocorrelation_3.png",species_info$file_name[j]),sep=""))
  # dev.off()

#   # # #####################################################
#   # # # Compute difference smooths using simultaneous CI's #  ----
#   # # #####################################################
  combi=c(1,2,3,4,5,6,7,8,9,10)
  combi_width=as.data.frame(combinations(n=10,r=2,v=combi))

  diff_df=data.frame(matrix(0, ncol = 3, nrow = 45))
  names(diff_df)=c("combi","period_diff","ci_width")

  for (i in 1:nrow(combi_width)){
    print(i)
    pdat=expand.grid("Year"=seq(1977,2019,length = 126),
                       "Seascapenr" = c(1,2,3,4,5,6,7,8,9,10),
                       "Ship" = "WAH3","DayNight"="D")

    xp <- predict(mod, pdat, type = 'lpmatrix',exclude=c("s(Ship)","s(DayNight)"))
    col1 <- grepl(sprintf("Seascapenr%s",combi_width$V1[i]), colnames(xp))
    col2 <- grepl(sprintf("Seascapenr%s",combi_width$V2[i]), colnames(xp))
    row1 <- pdat[["Seascapenr"]] == combi_width$V1[i]
    row2 <- pdat[["Seascapenr"]] == combi_width$V2[i]

    ## difference rows of pred for data from comparison
    X <- xp[row1, ] - xp[row2, ]
    ## zero out cols of X related to splines for other lochs
    X[, ! (col1 | col2)] <- 0
    ## zero out the parametric cols
    X[, !grepl('^s\\(', colnames(xp))] <- 0

    #predicted values for differences
    Vb=vcov(mod)
    dif <- X %*% coef(mod)
    se <- sqrt(rowSums((X %*% Vb) * X))
    pred=cbind(dif,se)

    ## CI of the estimates
    set.seed(42)
    N <- 10000

    #predi <- predict(m1, pdat, se.fit = TRUE)
    #se.fit <- predi$se.fit

    BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
    Cg <- predict(mod, pdat, type = "lpmatrix",exclude=c("s(Ship)","s(DayNight)"))
    simDev <- Cg %*% t(BUdiff)

    absDev <- abs(sweep(simDev, 1, se, FUN = "/"))
    masd <- apply(absDev, 2L, max)
    crit <- quantile(masd, prob = 0.95, type = 8)

    newd <- expand.grid(Year = seq(1977, 2019, length = 126))

    S12pred <- transform(cbind(data.frame(pred), newd),
                         uprP = V1 + (2 * se),
                         lwrP = V1 - (2 * se),
                         uprS = V1 + (crit * se),
                         lwrS = V1 - (crit * se))
    S12pred$sigdif <-ifelse(S12pred$lwrS<0 & S12pred$uprS<0,"lower",ifelse(S12pred$lwrS>0 & S12pred$uprS>0,"higher","no"))
    S12higher=S12pred
    S12higher["V1"][S12pred["sigdif"]!="higher"]=NA
    S12lower=S12pred
    S12lower["V1"][S12pred["sigdif"]!="lower"]=NA

    #nrow(S12pred)
    #length(which(S12pred$sigdif=="higher"))
    #length(which(S12pred$sigdif=="lower"))

    per_diff=(length(which(S12pred$sigdif=="higher"))+length(which(S12pred$sigdif=="lower")))/nrow(S12pred)
    S12pred$ci_width=S12pred$uprS-S12pred$lwrS
    ci_width=mean(S12pred$ci_width)

    diff_df$combi[i]=  sprintf("%s - %s",combi_width$V1[i],combi_width$V2[i])
    diff_df$period_diff[i]=per_diff
    diff_df$ci_width[i]=ci_width


    DifS12plot=
      ggplot(S12pred, aes(x = Year, y = V1)) +
      theme_bw()+
      #geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
      geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
      geom_line(data=S12pred,aes(x=Year,y=V1),col="black",size=0.5) +
      geom_line(data=S12higher,aes(x=Year,y=V1),col="seagreen1",size=1)+
      geom_line(data=S12lower,aes(x=Year,y=V1),col="red",size=1)+
      ggtitle(sprintf("%s Seascape %s - %s",species_info$Scientific.name[j],combi_width$V1[i],combi_width$V2[i]))+
      labs(x = "Year", y = 'Difference Fish biomass trend')+
      theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
      theme(axis.text.x = element_text(face="bold",size=9))+
      theme(axis.text.y = element_text(face="bold",size=9))+
      theme(axis.title.x = element_text(face="bold",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
      theme(axis.title.y = element_text(face="bold",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
      theme(axis.line = element_line(size=1))
      ggsave(DifS12plot,file = paste(path,sprintf("Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooth_%s_%s.png",species_info$file_name[j],combi_width$V1[i],combi_width$V2[i]),sep=""), height=3.5, width=4.5,dpi = 600)
  #DifS12plot
#   #
   }
   write.csv(diff_df,file = paste(path,sprintf("Seascapes/Q1/%s/GAM/Diff_smooth/num_diff.csv",species_info$file_name[j]),sep=""))
}

#save different smooths in single large image            
for (j in 1:16){
  print(j)
  path
  p=paste(path,sprintf("Seascapes/Q1/%s/GAM/Diff_smooth",species_info$file_name[j]),sep="")
  #p
  filenames=list.files(path=p,pattern="*.png",full.names=TRUE)
  #filenames
  rl <- lapply(filenames, png::readPNG)
  gl <- lapply(rl, grid::rasterGrob)
  cp=grid.arrange(gl[[1]],gl[[2]],gl[[3]],gl[[4]],gl[[5]],gl[[6]],
                  gl[[7]],gl[[8]],gl[[9]],gl[[10]],gl[[11]],gl[[12]],
                  gl[[13]],gl[[14]],gl[[15]],gl[[16]],gl[[17]],gl[[18]],
                  gl[[19]],gl[[20]],gl[[21]],gl[[22]],gl[[23]],gl[[24]],
                  gl[[25]],gl[[26]],gl[[27]],gl[[28]],gl[[29]],gl[[30]],
                  gl[[31]],gl[[32]],gl[[33]],gl[[34]],gl[[35]],gl[[36]],
                  gl[[37]],gl[[38]],gl[[39]],gl[[40]],gl[[41]],gl[[42]],
                  gl[[43]],gl[[44]],gl[[45]],ncol=9,nrow=5)
  ggsave(cp,file=paste(path,sprintf("Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooths.png",species_info$file_name[j]),sep=""),height=15, width=25,dpi = 300)
}


###########################################
#### compute temporal correlation file #### ----
###########################################
library(MFDFA)
library(mpmi)
#selected species for temporal 
for (k in 1:1){# in c(1,2,3,4,5,7,12,13,14)){
  data=read.csv(paste(path,sprintf("Data/Filtered_%s.csv",species_info$file_name[k]),sep=""))
  names(data)[names(data) == "id"] <- "Seascapenr"
  print(species_info$file_name[k])
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  ssc_data=completeFun(data,"Seascapenr")
  ssc_data=subset(ssc_data,Quarter==1)
  data_agg=ssc_data
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species+Ship+DayNight, data=ssc_data, sum)
  #print(nrow(data_agg))
  #positive_data_agg=subset(data_agg, Total_wgt!= 0)
  #print(nrow(positive_data_agg))
  
  # # # # #set outlier weight limit (99th percentile)
  #   nf_p=quantile(data_agg$Total_wgt,0.95)
  #   if (nf_p > 0){
  #    #subset of data excluding outlier weights
  #     data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
  #   }
  # # 
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000
  
  #separate subset per seascape, anscombe transform, max scale it and finally rbind back into single dataset----
  #1----
  data_agg1=subset(data_agg,Seascapenr==1)
  data_agg1$ansc_biomass= 2*nthroot(data_agg1$biomass_kg+(3/8),2)
  data_agg1$scaled_ansc_biomass=scaleTR(data_agg1$ansc_biomass)
  data_agg1$scaled_ansc_biomass[is.na(data_agg1$scaled_ansc_biomass)] <- 0.0001
  
  #2-----
  data_agg2=subset(data_agg,Seascapenr==2)
  data_agg2$ansc_biomass= 2*nthroot(data_agg2$biomass_kg+(3/8),2)
  data_agg2$scaled_ansc_biomass=scaleTR(data_agg2$ansc_biomass)
  data_agg2$scaled_ansc_biomass[is.na(data_agg2$scaled_ansc_biomass)] <- 0.0001
  
  #3----
  data_agg3=subset(data_agg,Seascapenr==3)
  data_agg3$ansc_biomass= 2*nthroot(data_agg3$biomass_kg+(3/8),2)
  data_agg3$scaled_ansc_biomass=scaleTR(data_agg3$ansc_biomass)
  data_agg3$scaled_ansc_biomass[is.na(data_agg3$scaled_ansc_biomass)] <- 0.0001
  
  #4----
  data_agg4=subset(data_agg,Seascapenr==4)
  data_agg4$ansc_biomass= 2*nthroot(data_agg4$biomass_kg+(3/8),2)
  data_agg4$scaled_ansc_biomass=scaleTR(data_agg4$ansc_biomass)
  data_agg4$scaled_ansc_biomass[is.na(data_agg4$scaled_ansc_biomass)] <- 0.0001
  
  #5----
  data_agg5=subset(data_agg,Seascapenr==5)
  data_agg5$ansc_biomass= 2*nthroot(data_agg5$biomass_kg+(3/8),2)
  data_agg5$scaled_ansc_biomass=scaleTR(data_agg5$ansc_biomass)
  data_agg5$scaled_ansc_biomass[is.na(data_agg5$scaled_ansc_biomass)] <- 0.0001
  
  
  #6----
  data_agg6=subset(data_agg,Seascapenr==6)
  data_agg6$ansc_biomass= 2*nthroot(data_agg6$biomass_kg+(3/8),2)
  data_agg6$scaled_ansc_biomass=scaleTR(data_agg6$ansc_biomass)
  data_agg6$scaled_ansc_biomass[is.na(data_agg6$scaled_ansc_biomass)] <- 0.0001
  
  #7----
  data_agg7=subset(data_agg,Seascapenr==7)
  data_agg7$ansc_biomass= 2*nthroot(data_agg7$biomass_kg+(3/8),2)
  data_agg7$scaled_ansc_biomass=scaleTR(data_agg7$ansc_biomass)
  data_agg7$scaled_ansc_biomass[is.na(data_agg7$scaled_ansc_biomass)] <- 0.0001
  
  #8----
  data_agg8=subset(data_agg,Seascapenr==8)
  data_agg8$ansc_biomass= 2*nthroot(data_agg8$biomass_kg+(3/8),2)
  data_agg8$scaled_ansc_biomass=scaleTR(data_agg8$ansc_biomass)
  data_agg8$scaled_ansc_biomass[is.na(data_agg8$scaled_ansc_biomass)] <- 0.0001
  
  #9----
  data_agg9=subset(data_agg,Seascapenr==9)
  data_agg9$ansc_biomass= 2*nthroot(data_agg9$biomass_kg+(3/8),2)
  data_agg9$scaled_ansc_biomass=scaleTR(data_agg9$ansc_biomass)
  data_agg9$scaled_ansc_biomass[is.na(data_agg9$scaled_ansc_biomass)] <- 0.0001
  
  #10----
  data_agg10=subset(data_agg,Seascapenr==10)
  data_agg10$ansc_biomass= 2*nthroot(data_agg10$biomass_kg+(3/8),2)
  data_agg10$scaled_ansc_biomass=scaleTR(data_agg10$ansc_biomass)
  data_agg10$scaled_ansc_biomass[is.na(data_agg10$scaled_ansc_biomass)] <- 0.0001

  ###############################################################################
  data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9,data_agg10)
  Mean_data=aggregate(biomass_kg~Year+Seascapenr,data=data_agg_ansc,mean)
  combi=c(1,2,3,4,5,6,7,8,9,10)
  combi_width=as.data.frame(combinations(n=10,r=2,v=combi))

  lagged_MI=data.frame(matrix(ncol=6,nrow=0))
  colnames(lagged_MI)=c("Seascapes","lag","si","cmi","bcmi","z")
  cnames=c("ts1","ts2","lag1","lag2","lag3","lag4","lag5","lag6","lag7","lag8","lag9","lag10")
  lag_list=c(0,1,2,3,4,5,6,7,8,9,10)
  
  #TRY FOR FITTED POSTERIOR SMOOTH TRENDS
  
  
  for (j in 1:nrow(combi_width)){
    V1=combi_width$V1[j]
    V2=combi_width$V2[j]
    ts1=subset(Mean_data,Seascapenr==V1)$biomass_kg
    ts2=subset(Mean_data,Seascapenr==V2)$biomass_kg

    dat=data.matrix(cbind(ts1,ts2,lead(ts2,1),lead(ts2,2),lead(ts2,3),lead(ts2,4),lead(ts2,5),
                          lead(ts2,6),lead(ts2,7),lead(ts2,8),lead(ts2,9),lead(ts2,10)))

    colnames(dat)=cnames
    c_var=paste(V1,"&",V2,sep=" ")
    print(c_var)

    
    for (i in lag_list){
      data=dat[,c(1,(2+i))]
      #dmx=discretize2d(data[,1],data[,2],numBins1=1000,numBins2=1000,r1=c(0,max(ts1)),r2=c(0,max(ts2)))#,r1=c(0,1),r2=c(0,1))
      #H1=entropy(rowSums(dmx))
      #H2=entropy(colSums(dmx))
      #H12=entropy(dmx)
      #mi=H1+H2-H12
      #miq=mi/H1
      mi=cmi.pw(data[,1],data[,2])$mi
      si=cmi.pw(data[,1],data[,1])$mi
      miq=mi/si
      bcmi=cmi.pw(data[,1],data[,2])$bcmi
      miz=cmi.pw(data[,1],data[,2])$z
      df=data.frame(c_var,i,si,mi,bcmi,miz)
      colnames(df)=c("Seascapes","lag","si","cmi","bcmi","z")
      lagged_MI=rbind(lagged_MI,df)
     }
  }
  write.csv(lagged_MI,paste(path,sprintf("Seascapes/Q1/%s/GAM/MI/mutinf_lag.csv",species_info$file_name[k]),sep=""))
}

#   ####################################################
#   #GAM model with separate smoother per seascape
# # ###################################################
#     print("gam")
#     data_agg_ansc$Seascapenr=as.factor(data_agg_ansc$Seascapenr) #Set as factor
#     data_agg_ansc = start_event(data_agg_ansc, column="Year",event="Seascapenr") #Sep time series per factor level
#   
#   
#     mod=bam(scaled_ansc_biomass~Seascapenr+s(Year, by=Seascapenr)+s(Ship,bs="re")+s(DayNight,bs="re"),data=data_agg_ansc,method="REML",family=betar)
#     print("gam fitted")
#   
#     combi=c(1,2,3,4,5,6,7,8,9,10)
#     combi_width=as.data.frame(combinations(n=10,r=2,v=combi))
#     combi_width2=as.matrix(combinations(n=10,r=2,v=combi))
#   
#     #newdat = as.data.frame(1977:2019)
#     #colnames(newdat)="Year"
#     #fitted_samples()
#     #100 simulations from model posterior distribution including sampling variation
#     sim_dat=data.frame(smooth_samples(mod,n=1000,term="Seascapenr",n_vals=43))#scale="response"))
#     sim_dat$Year=sim_dat$.x1
#     #sim_dat$Year=rep(data_agg_ansc$Year,1000)
#     #sim_dat$Seascapenr=rep(data_agg_ansc$Seascapenr,1000)
#     #sim_dat2=aggregate(response~Year+draw+Seascapenr,sim_dat,mean)

#   
#   
#   #difference 1ce and look at ccf
#   autocor=apply(combi_width2,MARGIN=1,
#                 FUN=function(x){
#                   V1=x[1]
#                   #print(V1)
#                   V2=x[2]
#                   #print(V2)
#                   c_var=paste(V1,"&",V2,sep=" ")
#                   print(c_var)
# 
#                   #sv1=sim_dat2[sim_dat2[,"Seascapenr"]==V1,]
#                   sv1=subset(sim_dat2,Seascapenr==V1)
#                   #sv2=sim_dat2[sim_dat2[,"Seascapenr"]==V2,]
#                   sv2=subset(sim_dat2,Seascapenr==V2)
# 
#                   simulations=c(seq(1:1000))
# 
#                   sapply(simulations,FUN=function(y){
# 
#                     draw_v1=subset(sv1,draw==y[1])
#                     draw_v2=subset(sv2,draw==y[1])
# 
#                     dv1=draw_v1$response
#                     dv2=draw_v2$response
# 
#                     sim=as.numeric(y[1])
# 
# 
#                     b=ccf(dv1,dv2,lag.max=10,plot=FALSE) #temporal correlation
#                     acf=b$acf
#                     lag=b$lag
#                     Seascapes=rep(c_var,length(acf))
#                     simu=rep(sim,length(acf))
#                     dat=data.frame(cbind(acf,lag,Seascapes,simu))
#                     return(dat)
#                     })
# 
#                 })
#   # 
#   # 
#   #sim_dat=as.data.frame(smooth_samples(mod,n=100,term="Seascapenr",n_vals=43))
#   #sim_dat2=data.matrix(smooth_samples(mod,n=100,term="Seascapenr",n_vals=43))
#   #sim_dat2=data.matrix(smooth_samples(mod,n=100,t,n_vals=43))
#   # 
#   # Year=floor(as.numeric(sim_dat2[,".x1"]))
#   # sim_dat2=cbind(sim_dat2,Year)
#   # #plot(x=sim_dat$Year,y=sim_dat$value)
#   # #sv1=sim_dat2[sim_dat2[,"Seascapenr"]==1,]
#   # #draw_v1=sv1[sv1[,"draw"]==1,]
#   # #go iteratively through smooth pair combinations and simulate autocorrelation
#   # sim_autocor=matrix(ncol=4,nrow=0)
#   # #draw_v1[,8]
#   # 
#   # 
#   # 
#   # 
#   # #difference 1ce and look at ccf
#   # autocor=apply(combi_width2,MARGIN=1,
#   #               FUN=function(x){
#   #                 V1=x[1]
#   #                 #print(V1)
#   #                 V2=x[2]
#   #                 #print(V2)
#   #                 c_var=paste(V1,"&",V2,sep=" ")
#   #                 print(c_var)
#   #                 
#   #                 #sv1=sim_dat2[sim_dat2[,"Seascapenr"]==V1,]
#   #                 sv1=subset(sim_dat2,Seascapenr==V1)
#   #                 #sv2=sim_dat2[sim_dat2[,"Seascapenr"]==V2,]
#   #                 sv2=subset(sim_dat2,Seascapenr==V2)
#   #                 
#   #                 simulations=c(seq(1:1000))
#   #                 
#   #                 sapply(simulations,FUN=function(y){
#   #                   
#   #                   draw_v1=subset(sv1,draw==y[1])
#   #                   draw_v2=subset(sv2,draw==y[1])
#   #                   
#   #                   dv1=draw_v1$response
#   #                   dv2=draw_v2$response
#   #                   
#   #                   sim=as.numeric(y[1])
#   #                   
#   #                   
#   #                   dat=q_L_AXHA(dv1,dv2,lags=c(0,1,2,3,4,5,6,7,8,9,10),L=5,q=2) #temporal correlation
#   #                   dat=cbind(dat,c(0,1,2,3,4,5,6,7,8,9,10)) # lags
#   #                   dat=cbind(dat,c(sim,sim,sim,sim,sim,sim,sim,sim,sim,sim,sim)) # draw
#   #                   dat=cbind(dat,c(c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var)) #seascapes
#   # # 
#   #                   return(dat)
#   #            })
#    
#     
#   #prepare as csv
#   df  <-  as.data.frame(autocor)
#   df=data.frame(t(df))
#   df=as.data.frame(lapply(df,unlist))
#   colnames(df)=c("Acf","lag","Seascapes","sim")
#   #df$lag_dup=NULL
#   write.csv(df,paste(path,sprintf("Seascapes/Q1/%s/GAM/CCF/autocor_lag_undif.csv",species_info$file_name[j]),sep=""))
 }


#######################################
### Create autocor graphs from file ###----
#######################################
library(tidyverse)
library(RColorBrewer)
library(colorRamps)
library(stringr)

#seascape combinations to recall
combi=c(1,2,3,4,5,6,7,8,9,10)
combi_width=as.data.frame(combinations(n=10,r=2,v=combi))


cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)

for (j in c(1,2,3,4,5,7,12,13,14)){
   data=read.csv(paste(path,sprintf("Seascapes/Q1/%s/GAM/MI/mutinf_lag.csv",species_info$file_name[j]),sep=""))
   #names(data)[1:5]=c("ID","Acf","lag","Seascapes","sim")
   print(species_info$file_name[j])
   #cols=rev("RdYlBu")
   
    #cols <- brewer.pal(11, 'Spectral')
    #cols
    #combi_data=aggregate(Acf~lag+Seascapes,data=data,FUN=mean)
    #add_cols=str_split_fixed(combi_data$Seascapes," & ",2)
    #combi_data$S1=as.numeric(add_cols[,1])
    #combi_data$S2=as.numeric(add_cols[,2])
    #combi_data$Seascapes2=as.numeric(gsub(" & ","",combi_data$Seascapes))
    # combi_data$Seascapes=gsub("&","-",combi_data$Seascapes)
    # combi_data$Seascapes=factor(combi_data$Seascapes, levels = c("1 - 2", "1 - 3", "1 - 4", "1 - 5", "1 - 6", "1 - 7", "1 - 8", "1 - 9", "1 - 10",
    #                                                              "2 - 3", "2 - 4", "2 - 5", "2 - 6", "2 - 7", "2 - 8", "2 - 9", "2 - 10", "3 - 4",
    #                                                              "3 - 5", "3 - 6", "3 - 7", "3 - 8", "3 - 9", "3 - 10", "4 - 5", "4 - 6", "4 - 7",
    #                                                              "4 - 8", "4 - 9", "4 - 10", "5 - 6","5 - 7", "5 - 8", "5 - 9", "5 - 10", "6 - 7",
    #                                                              "6 - 8", "6 - 9", "6 - 10", "7 - 8", "7 - 9", "7 - 10", "8 - 9", "8 - 10", "9 - 10"))
    # #
   data$Seascapes=gsub("&","-",data$Seascapes) 
   data$Seascapes=factor(data$Seascapes, levels = c("1 - 2", "1 - 3", "1 - 4", "1 - 5", "1 - 6", "1 - 7", "1 - 8", "1 - 9", "1 - 10",
                                                                "2 - 3", "2 - 4", "2 - 5", "2 - 6", "2 - 7", "2 - 8", "2 - 9", "2 - 10", "3 - 4",
                                                                "3 - 5", "3 - 6", "3 - 7", "3 - 8", "3 - 9", "3 - 10", "4 - 5", "4 - 6", "4 - 7",
                                                                "4 - 8", "4 - 9", "4 - 10", "5 - 6","5 - 7", "5 - 8", "5 - 9", "5 - 10", "6 - 7",
                                                                "6 - 8", "6 - 9", "6 - 10", "7 - 8", "7 - 9", "7 - 10", "8 - 9", "8 - 10", "9 - 10"))

# 
#    pal=rev(brewer.pal(n=5,name="RdBu"))
#    pal
   p=ggplot(data = data, aes(x=Seascapes, y=lag, fill = bcmi))+
      geom_tile(color = "white")+
      scale_fill_gradientn(colors=mypalette,limits=c(-0.1,1))+#(colours=pal,limits=c(-0.1,1)) + #rainbow(21)
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size=9))+
      theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1,size=9))+
      scale_y_continuous(breaks = round(seq(min(combi_data2$lag), max(combi_data2$lag), by = 1),1))
    p
   #  #p
   #  
   #   # q=ggplot(data = combi_data, aes(x=Seascapes, y=lag))+#, fill = Accf))+
   #   #  #geom_tile(color = "white")+
   #   #  geom_point(aes(colour=Acf,size=Acf))+
   #   #  scale_color_gradientn(colors=rainbow(11),limits=c(-1,1)) + #rainbow(21)
   #   #  theme_minimal()+
   #   #  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size=9))+
   #   #  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1,size=9))+
   #   #  scale_y_continuous(breaks = round(seq(min(combi_data2$lag), max(combi_data2$lag), by = 1),1))
   #   # #q
     ggsave(p,file = paste(path,sprintf("Seascapes/Q1/%s/GAM/MI/Tiled_bcmi_%s.png",species_info$file_name[j],species_info$file_name[j]),sep=""), height=2.5, width=8,dpi = 600)
    #ggsave(q,file = paste(path,sprintf("Seascapes/Q1/%s/GAM/CCF/diff1_Circled_ccf_%s.png",species_info$file_name[j],species_info$file_name[j]),sep=""), height=6, width=8,dpi = 600)
#   

   # for (i in 1:nrow(combi_width)){
   #  sc_combi=sprintf("%s - %s",combi_width$V1[i],combi_width$V2[i])
   # 
   #  #subset data
   #  combi_data=subset(data, Seascapes==sc_combi)
   #  #combi_data=select(combi_data,-ID,sim)#drop unneccesary columns
   # 
   #  # aggregate autocorrelation per lag (mean per lag) for subset
   #  #combi_data2=do.call(data.frame, aggregate(.~lag+Seascapes,data=combi_data,function(x)c(mean=mean(x),se=std_err(x))))
   #  #combi_data2$upr_conf=combi_data2$Acf.mean+(2*combi_data2$Acf.se)
   #  #combi_data2$lwr_conf=combi_data2$Acf.mean-(2*combi_data2$Acf.se)
   # 
   #  #create graph
   #  autocor_plot=ggplot(combi_data, aes(x = lag, y = bcmi))+
   #               theme_classic()+
   #               geom_bar(stat="identity",col="black",fill="gray",alpha=0.5,width=0.5)+#col="black",fill="gray",width=0.5)+
   #               #geom_errorbar(aes(ymin=lwr_conf,ymax=upr_conf),width=0.2)+
   # 
   #               #geom_hline(yintercept=0.3, linetype="dashed", color = "red",size=1,alpha=0.5)+
   #               #geom_hline(yintercept=-0.3, linetype="dashed", color = "red",size=1,alpha=0.5)+
   #               geom_hline(yintercept=0, color = "black",size=0.5)+
   # 
   # 
   #               ggtitle(sprintf("%s Mutual information S%s & S%s",species_info$Scientific.name[j],combi_width$V1[i],combi_width$V2[i]))+
   #               labs(x = "Lag (Year)", y = 'BCMI')+
   #               scale_x_continuous(breaks = round(seq(min(combi_data$lag), max(combi_data$lag), by = 1),1))+
   #               #scale_y_continuous(breaks = round(seq(min(combi_data$cor), max(combi_data$cor), by = 0.05),1))+
   # 
   #                theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
   #                theme(axis.text.x = element_text(face="bold",size=9))+
   #                theme(axis.text.y = element_text(face="bold",size=9))+
   #                theme(axis.title.x = element_text(face="bold",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
   #                theme(axis.title.y = element_text(face="bold",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
   #                ylim(-0.1,1)
   #  #autocor_plot
   #  ggsave(autocor_plot,file = paste(path,sprintf("Seascapes/Q1/%s/GAM/CCF/CCF_smooth_%s_%s_ANSC.png",species_info$file_name[j],combi_width$V1[i],combi_width$V2[i]),sep=""), height=3.5, width=4.5,dpi = 600)

  }
}


####################################################################
###figure to illustrate fitting of gam and posterior simulations####
####################################################################
#dataprep----
data=read.csv(paste(path,sprintf("Data/Filtered_%s.csv",species_info$file_name[1]),sep=""))
names(data)[names(data) == "id"] <- "Seascapenr"
print(species_info$file_name[1])
# Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
data$ID = cumsum(!duplicated(data[3:7]))

ssc_data=completeFun(data,"Seascapenr")
ssc_data=subset(ssc_data,Quarter==1)
data_agg=ssc_data
data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species+Ship+DayNight, data=ssc_data, sum)
#print(nrow(data_agg))
#positive_data_agg=subset(data_agg, Total_wgt!= 0)
#print(nrow(positive_data_agg))

# # # # #set outlier weight limit (99th percentile)
#   nf_p=quantile(data_agg$Total_wgt,0.95)
#   if (nf_p > 0){
#    #subset of data excluding outlier weights
#     data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
#   }
# # 

# convert weights from grams to kg
data_agg$biomass_kg=data_agg$Total_wgt/1000

#separate subset per seascape, anscombe transform, max scale it and finally rbind back into single dataset----
#1----
data_agg1=subset(data_agg,Seascapenr==1)
data_agg1$ansc_biomass= 2*nthroot(data_agg1$biomass_kg+(3/8),2)
data_agg1$scaled_ansc_biomass=scaleTR(data_agg1$ansc_biomass)
data_agg1$scaled_ansc_biomass[is.na(data_agg1$scaled_ansc_biomass)] <- 0.0001

#2-----
data_agg2=subset(data_agg,Seascapenr==2)
data_agg2$ansc_biomass= 2*nthroot(data_agg2$biomass_kg+(3/8),2)
data_agg2$scaled_ansc_biomass=scaleTR(data_agg2$ansc_biomass)
data_agg2$scaled_ansc_biomass[is.na(data_agg2$scaled_ansc_biomass)] <- 0.0001

#3----
data_agg3=subset(data_agg,Seascapenr==3)
data_agg3$ansc_biomass= 2*nthroot(data_agg3$biomass_kg+(3/8),2)
data_agg3$scaled_ansc_biomass=scaleTR(data_agg3$ansc_biomass)
data_agg3$scaled_ansc_biomass[is.na(data_agg3$scaled_ansc_biomass)] <- 0.0001

#4----
data_agg4=subset(data_agg,Seascapenr==4)
data_agg4$ansc_biomass= 2*nthroot(data_agg4$biomass_kg+(3/8),2)
data_agg4$scaled_ansc_biomass=scaleTR(data_agg4$ansc_biomass)
data_agg4$scaled_ansc_biomass[is.na(data_agg4$scaled_ansc_biomass)] <- 0.0001

#5----
data_agg5=subset(data_agg,Seascapenr==5)
data_agg5$ansc_biomass= 2*nthroot(data_agg5$biomass_kg+(3/8),2)
data_agg5$scaled_ansc_biomass=scaleTR(data_agg5$ansc_biomass)
data_agg5$scaled_ansc_biomass[is.na(data_agg5$scaled_ansc_biomass)] <- 0.0001


#6----
data_agg6=subset(data_agg,Seascapenr==6)
data_agg6$ansc_biomass= 2*nthroot(data_agg6$biomass_kg+(3/8),2)
data_agg6$scaled_ansc_biomass=scaleTR(data_agg6$ansc_biomass)
data_agg6$scaled_ansc_biomass[is.na(data_agg6$scaled_ansc_biomass)] <- 0.0001

#7----
data_agg7=subset(data_agg,Seascapenr==7)
data_agg7$ansc_biomass= 2*nthroot(data_agg7$biomass_kg+(3/8),2)
data_agg7$scaled_ansc_biomass=scaleTR(data_agg7$ansc_biomass)
data_agg7$scaled_ansc_biomass[is.na(data_agg7$scaled_ansc_biomass)] <- 0.0001

#8----
data_agg8=subset(data_agg,Seascapenr==8)
data_agg8$ansc_biomass= 2*nthroot(data_agg8$biomass_kg+(3/8),2)
data_agg8$scaled_ansc_biomass=scaleTR(data_agg8$ansc_biomass)
data_agg8$scaled_ansc_biomass[is.na(data_agg8$scaled_ansc_biomass)] <- 0.0001

#9----
data_agg9=subset(data_agg,Seascapenr==9)
data_agg9$ansc_biomass= 2*nthroot(data_agg9$biomass_kg+(3/8),2)
data_agg9$scaled_ansc_biomass=scaleTR(data_agg9$ansc_biomass)
data_agg9$scaled_ansc_biomass[is.na(data_agg9$scaled_ansc_biomass)] <- 0.0001

#10----
data_agg10=subset(data_agg,Seascapenr==10)
data_agg10$ansc_biomass= 2*nthroot(data_agg10$biomass_kg+(3/8),2)
data_agg10$scaled_ansc_biomass=scaleTR(data_agg10$ansc_biomass)
data_agg10$scaled_ansc_biomass[is.na(data_agg10$scaled_ansc_biomass)] <- 0.0001


#rbind----
data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9,data_agg10)

####################################################
#GAM model with separate smoother per seascape
###################################################
print("gam")
data_agg_ansc$Seascapenr=as.factor(data_agg_ansc$Seascapenr) #Set as factor
data_agg_ansc = start_event(data_agg_ansc, column="Year",event="Seascapenr") #Sep time series per factor level


mod=bam(scaled_ansc_biomass~Seascapenr+s(Year, by=Seascapenr)+s(Ship,bs="re")+s(DayNight,bs="re"),data=data_agg_ansc,method="REML",family=betar)

sim_dat=as.data.frame(predicted_samples(mod,n=1000,scale="response"))
sim_dat$Seascapenr=rep(data_agg_ansc$Seascapenr,1000)
sim_dat$Year=rep(data_agg_ansc$Year,1000)
sim_dat2=aggregate(response~Year+draw+Seascapenr,sim_dat,mean)

###
Vb <- vcov(mod)
newd=expand.grid("Year"=seq(1977,2019,length = 43),
                 "Seascapenr" = c(1,2,3,4,5,6,7,8,9,10),
                 "Ship" = "WAH3",
                 "DayNight" ="D")
pred <- predict(mod, newd, se.fit = TRUE,exclude=c("s(Ship)","s(DayNight)"))
se.fit <- pred$se.fit

set.seed(42)
N <- 10000

BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

#xp matrix where basisfunctions of the model have been evaluated at 400 time point values per area
Cg <- predict(mod, newd, type = "lpmatrix",exclude="s(Ship)","s(DayNight)")
simDev <- Cg %*% t(BUdiff)

absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
masd <- apply(absDev, 2L, max)
crit <- quantile(masd, prob = 0.95, type = 8)
pred <- transform(cbind(data.frame(pred), newd),
                  uprP = fit + (2 * se.fit),
                  lwrP = fit - (2 * se.fit),
                  uprS = fit + (crit * se.fit),
                  lwrS = fit - (crit * se.fit))
#separate plots per area
S1=pred[which(pred$Seascapenr==1),]
S2=pred[which(pred$Seascapenr==2),]


#subset data for seascape 1 and 2
sim_dat_S1=subset(sim_dat2,Seascapenr==1)
sim_dat_S1$uprS=rep(S1$uprS,1000)#+2.8,100)
sim_dat_S1$lwrS=rep(S1$lwrS,1000)#+2.8,100)
sim_dat_S2=subset(sim_dat2,Seascapenr==2)

s1d1=subset(sim_dat_S1,draw==2)
s2d1=subset(sim_dat_S2,draw==2)

b=ccf(diff(s1d1$response,1),diff(s2d1$response,1),plot=FALSE)
acf=b$acf
lag=b$lag
Seascapes=rep(c_var,length(acf))
dat=cbind(acf,lag,Seascapes)
  
#GGPLOTS---
#sim1
s1_plot=ggplot(sim_dat_S1, aes(x = Year, y = response, group=draw))+
  theme_classic()+
  #geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.005,fill="blue")+
  geom_line(alpha=0.025)+
  ggtitle(sprintf("Posterior simulations seascape %s",sim_dat_S1$Seascapenr[1]))+
  labs(x = "Year", y = 'Fish biomass trend')+
  #scale_x_continuous(breaks = round(seq(min(combi_data2$lag), max(combi_data2$lag), by = 1),1))+
  theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
  theme(axis.text.x = element_text(face="bold",size=12))+
  theme(axis.text.y = element_text(face="bold",size=12))+
  theme(axis.title.x = element_text(face="bold",size=12,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  theme(axis.title.y = element_text(face="bold",size=12,margin = margin(t = 0, r = 10, b = 0, l = 0)))
s1_plot
ggsave(s1_plot,file = paste(path,sprintf("images/s1simulated_means1000.png"),sep=""), height=3.5, width=4.5,dpi = 600)

#sim2
s2_plot=ggplot(sim_dat_S2, aes(x = Year, y = response, group=draw))+
  theme_classic()+
  geom_line(alpha=0.025)+
  ggtitle(sprintf("Posterior simulations seascape %s",sim_dat_S2$Seascapenr[1]))+
  labs(x = "Year", y = 'Fish biomass trend')+
  #scale_x_continuous(breaks = round(seq(min(combi_data2$lag), max(combi_data2$lag), by = 1),1))+
  theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
  theme(axis.text.x = element_text(face="bold",size=12))+
  theme(axis.text.y = element_text(face="bold",size=12))+
  theme(axis.title.x = element_text(face="bold",size=12,margin = margin(t = 10, r = 0, b = 0, l = 5)))+
  theme(axis.title.y = element_text(face="bold",size=12,margin = margin(t = 0, r = 10, b = 0, l = 0)))
s2_plot
ggsave(s2_plot,file = paste(path,sprintf("images/s2simulated_means1000.png"),sep=""), height=3.5, width=4.5,dpi = 600)

?smooth_samples

getAnywhere(smooth_samples())
smooth_samples
showMethods(smooth_samples)
getMethod(smooth_samples)
