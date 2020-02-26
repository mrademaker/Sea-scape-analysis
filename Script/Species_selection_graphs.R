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
for (j in c(2:15)){
  print(j)
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
  
  #2-----
  data_agg2=subset(data_agg,Seascapenr==2)
  
  #3----
  data_agg3=subset(data_agg,Seascapenr==3)
  
  #4----
  data_agg4=subset(data_agg,Seascapenr==4)
  
  #5----
  data_agg5=subset(data_agg,Seascapenr==5)
  
  #6----
  data_agg6=subset(data_agg,Seascapenr==6)
  
  #7----
  data_agg7=subset(data_agg,Seascapenr==7)

  
  #8----
  data_agg8=subset(data_agg,Seascapenr==8)
  
  #9----
  data_agg9=subset(data_agg,Seascapenr==9)

  #10----
  data_agg10=subset(data_agg,Seascapenr==10)

  
  
  #rbind----
  data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9,data_agg10)
  
  ####################################################
  #GAM model with separate smoother per seascape
  ###################################################
  print("gam")
  data_agg_ansc$Seascapenr=as.factor(data_agg_ansc$Seascapenr) #Set as factor
  data_agg_ansc = start_event(data_agg_ansc, column="Year",event="Seascapenr") #Sep time series per factor level
  
  
  mod=bam(biomass_kg~Seascapenr+s(Year, by=Seascapenr)+s(Ship,bs="re")+s(DayNight,bs="re"),data=data_agg_ansc,method="REML")
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
  # ################################################################
  # Plot predicted smooths with all data                         #----
  # ################################################################
  # 
   newd=expand.grid("Year"=seq(1977,2019,length = 126),
              "Seascapenr" = c(1,2,3,4,5,6,7,8,9,10),
              "Ship" = "WAH3", "DayNight" ="D")
   
    pred <- predict(mod, newd, se.fit = TRUE,type="response",exclude=c("s(Ship)","s(DayNight)"))
    se.fit <- pred$se.fit
    pred <- transform(cbind(data.frame(pred), newd),
                      uprP = fit + (2 * se.fit),
                      lwrP = fit - (2 * se.fit))
  
  #separate plots per area
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
  # 
  img=readPNG(source=paste(path,sprintf("images/spec_%s.png",species_info$file_name[j]),sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  # # # #
  # # # # ### plot smooths with all data -----
   s1plot=ggplot(data=S1,aes(x=Year))+
     theme_bw()+
     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
     geom_point(aes(x=Year,y=biomass_kg),alpha=0.2,colour="blue",data=data_agg1)+
     geom_line(aes(y=fit),col="black",size=1)+
     ggtitle("S1")+
     labs(x="Year",y="Fish biomass trend")
   #Add species image
   s1plot=s1plot+annotation_custom(g,xmin = min(data_agg1$Year), xmax = min(data_agg1$Year)+9, ymin = max(data_agg1$biomass_kg)-0.15*(max(data_agg1$biomass_kg)), ymax = max(data_agg1$biomass_kg))
  # 
  s2plot=ggplot(data=S2,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=biomass_kg),alpha=0.2,colour="blue",data=data_agg2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S2")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s2plot=s2plot+annotation_custom(g,xmin = min(data_agg2$Year), xmax = min(data_agg2$Year)+8, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg))


  s3plot=ggplot(data=S3,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=biomass_kg),alpha=0.2,colour="blue",data=data_agg3)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S3")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s3plot=s3plot+annotation_custom(g,xmin = min(data_agg3$Year), xmax = min(data_agg3$Year)+8, ymin = max(data_agg3$biomass_kg)-0.15*(max(data_agg3$biomass_kg)), ymax = max(data_agg3$biomass_kg))

  s4plot=ggplot(data=S4,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=biomass_kg),alpha=0.2,colour="blue",data=data_agg4)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S4")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s4plot=s4plot+annotation_custom(g,xmin = min(data_agg4$Year), xmax = min(data_agg4$Year)+8, ymin = max(data_agg4$biomass_kg)-0.15*(max(data_agg4$biomass_kg)), ymax = max(data_agg4$biomass_kg))


  s5plot=ggplot(data=S5,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=biomass_kg),alpha=0.2,colour="blue",data=data_agg5)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S5")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s5plot=s5plot+annotation_custom(g,xmin = min(data_agg5$Year), xmax = min(data_agg5$Year)+8, ymin = max(data_agg5$biomass_kg)-0.15*(max(data_agg5$biomass_kg)), ymax = max(data_agg5$biomass_kg))


  s6plot=ggplot(data=S6,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=biomass_kg),alpha=0.2,colour="blue",data=data_agg6)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S6")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s6plot=s6plot+annotation_custom(g,xmin = min(data_agg6$Year), xmax = min(data_agg6$Year)+8, ymin = max(data_agg6$biomass_kg)-0.15*(max(data_agg6$biomass_kg)), ymax = max(data_agg6$biomass_kg))


  s7plot=ggplot(data=S7,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=biomass_kg),alpha=0.2,colour="blue",data=data_agg7)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S7")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s7plot=s7plot+annotation_custom(g,xmin = min(data_agg7$Year), xmax = min(data_agg7$Year)+8, ymin = max(data_agg7$biomass_kg)-0.15*(max(data_agg7$biomass_kg)), ymax = max(data_agg7$biomass_kg))

  s8plot=ggplot(data=S8,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=biomass_kg),alpha=0.2,colour="blue",data=data_agg8)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S8")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s8plot=s8plot+annotation_custom(g,xmin = min(data_agg8$Year), xmax = min(data_agg8$Year)+8, ymin = max(data_agg8$biomass_kg)-0.15*(max(data_agg8$biomass_kg)), ymax = max(data_agg8$biomass_kg))


  s9plot=ggplot(data=S9,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=biomass_kg),alpha=0.2,colour="blue",data=data_agg9)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S9")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s9plot=s9plot+annotation_custom(g,xmin = min(data_agg9$Year), xmax = min(data_agg9$Year)+8, ymin = max(data_agg9$biomass_kg)-0.15*(max(data_agg9$biomass_kg)), ymax = max(data_agg9$biomass_kg))

  s10plot=ggplot(data=S10,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=biomass_kg),alpha=0.2,colour="blue",data=data_agg10)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S10")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s10plot=s10plot+annotation_custom(g,xmin = min(data_agg10$Year), xmax = min(data_agg10$Year)+8, ymin = max(data_agg10$biomass_kg)-0.15*(max(data_agg10$biomass_kg)), ymax = max(data_agg10$biomass_kg))


  cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,s10plot,ncol=4,nrow=3)
  ggsave(cplot,file=paste(path,sprintf("Seascapes/Q1/%s/GAM/GAM_ANSC_data_smooths.png",species_info$file_name[j]),sep=""),height=8, width=10,dpi = 600)#
}
