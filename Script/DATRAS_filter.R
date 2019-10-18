#load packages----
library(purrr)
library(plyr)
library(tidyr)
library(pbapply)
library(dplyr)
library(data.table)
library(ggplot2)
library(magick)
library(here) # For making the script run without a wd
library(magrittr) # For piping the logo
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

st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

lagged_DCCA_CC=function(x,y,k,lags){
  
  ## calculate cumulated deviation time series Xt
  xx<- cumsum(x - mean(x))
  yy<- cumsum(y - mean(y))
  
  
  data=data.frame(time=seq(1:length(x)))
  data$time=as.numeric(seq(1:length(x)))
  
  data$xx=xx
  data$yy=yy
  
  store_f2_x=list()
  store_f2_y=list()
  store_f2_xy=list()
  
  lag_list=lags
  dcca_list=list()
  
  for (tau in lag_list){
    #### equation 2. fit linear model on sliding window
    for (i in 1:length(xx)){
      # so sliding window step size 1 untill i + window size reaches end of time series
      if (i + k + tau < length(xx)){
        #print(i)
        dat_x=data[i:(i+k),] #subset of x data in window
        dat_y=data[(i+tau):(i+k+tau),] #subset of y data in window
        
        #for x
        x_act=dat_x$xx                  #data values in window
        mod_x=lm(xx~time,data=dat_x)    #fit linear model
        x_hat = mod_x$fitted.values     #fitted values x in window
        f2DCCA_x = sum(((x_act-x_hat)^2))/(k+1)    #equation 2 
        store_f2_x = as.numeric(append(store_f2_x,f2DCCA_x))#add to list
        
        #for y
        y_act=dat_y$yy
        mod_y=lm(yy~time,data=dat_y)
        y_hat = mod_y$fitted.values
        f2DCCA_y = sum(((y_act-y_hat)^2))/(k+1)
        store_f2_y = as.numeric(append(store_f2_y,f2DCCA_y))
        
        #for xy
        f2DCCA_xy =sum(((x_act-x_hat)*(y_act-y_hat)))/(k+1) #equation 4
        store_f2_xy = as.numeric(append(store_f2_xy,f2DCCA_xy))
      }
    }
    F_DCCA_x=sum(store_f2_x)/(length(store_f2_x))#-k)       #(equation3)
    F_DCCA_y=sum(store_f2_y)/(length(store_f2_y))#-k) 
    #F_DFA_y=mean(store_f2_y)
    F_DCCA_xy=sum(store_f2_xy)/(length(store_f2_xy))#-k)   #equation 5.
    #F_DCCA_xy=mean(store_f2_xy)
    ptauDCCA = F_DCCA_xy/(sqrt(F_DCCA_x)*sqrt(F_DCCA_y)) #equation 1.
    dcca_list=as.numeric(append(dcca_list,ptauDCCA))
  }
  
  ptauDCCA=matrix(dcca_list)#as.data.frame(dcca_list)
  #ptauDCCA$lag=(0:(length(dcca_list)-1))
  #names(ptauDCCA)[names(ptauDCCA)=="dcca_list"]="corr"
  return(ptauDCCA)
}

#Time lagged detrended cross-correlation coefficient
q_L_AXHA=function(x,y,L,q,lags){

  ## calculate cumulated deviation time series Xt
  xx<- cumsum(x - mean(x))  ## Equation 2
  yy<- cumsum(y - mean(y))  ## Equation 2

  data=data.frame(time=seq(1:length(x)))
  data$time=as.numeric(seq(1:length(x)))
  data$xx=xx
  data$yy=yy

  store_fqXY=list()
  store_fqXX=list()
  store_fqYY=list()
  lag_list=lags
  AXHA_list=list()

  for (tau in lag_list){
    for (i in 1:length(xx)){
      if (i + L + tau < length(xx)){
        #dat_x=data[i:(i+L),]
        #dat_y=data[(i+tau):(i+L+tau),]

        #equation 2
        XY=(xx[i]-xx[i+L])*(yy[i+tau]-yy[i+tau+L])

        #equation 3
        fqXY=sign(XY)*(abs(XY)^(q/2))

        fqXX = (xx[i]-xx[i+L])^q
        fqYY = (yy[i+tau]-yy[i+L+tau])^q

        store_fqXY=as.numeric(append(store_fqXY,fqXY))
        store_fqXX=as.numeric(append(store_fqXX,fqXX))
        store_fqYY=as.numeric(append(store_fqYY,fqYY))
      }
    }
    FqXY=mean(store_fqXY)
    FqXX=mean(store_fqXX)
    FqYY=mean(store_fqYY)

    FqXX*FqYY
    sqrt(FqXX*FqYY)
    pqLAXHA=FqXY/sqrt(FqXX*FqYY)
    AXHA_list=as.numeric(append(AXHA_list,pqLAXHA))
  }
  pqLAXHA=as.data.frame(AXHA_list)
  pqLAXHA$lag=(0:(length(AXHA_list)-1))
  names(pqLAXHA)[names(pqLAXHA)=="AXHA_list"]="corr"
  return(pqLAXHA)
}


q_L_AXHA=function(x,y,L,q,lags){
  
  ## calculate cumulated deviation time series Xt
  xx<- cumsum(x - mean(x))  ## Equation 2
  yy<- cumsum(y - mean(y))  ## Equation 2
  
  data=data.frame(time=seq(1:length(x)))
  data$time=as.numeric(seq(1:length(x)))
  data$xx=xx
  data$yy=yy
  
  store_fqXY=list()
  store_fqXX=list()
  store_fqYY=list()
  lag_list=lags
  AXHA_list=list()
  
  for (tau in lag_list){  
    for (i in 1:length(xx)){
      if (i + L + tau < length(xx)){
        #dat_x=data[i:(i+L),]
        #dat_y=data[(i+tau):(i+L+tau),]
        
        #equation 2
        XY=(xx[i]-xx[i+L])*(yy[i+tau]-yy[i+tau+L])
        
        #equation 3
        fqXY=sign(XY)*((abs(XY))^(q/2))
        
        fqXX = (xx[i]-xx[i+L])^q
        fqYY = (yy[i]-yy[i+L+tau])^q
        
        store_fqXY=as.numeric(append(store_fqXY,fqXY))
        store_fqXX=as.numeric(append(store_fqXX,fqXX))
        store_fqYY=as.numeric(append(store_fqYY,fqYY))
      }
    }
    FqXY=mean(store_fqXY)
    FqXX=mean(store_fqXX)
    FqYY=mean(store_fqYY)
    
    FqXX*FqYY
    sqrt(FqXX*FqYY)
    pqLAXHA=FqXY/sqrt(FqXX*FqYY)
    AXHA_list=as.numeric(append(AXHA_list,pqLAXHA))
  }
  pqLAXHA=as.data.frame(AXHA_list)
  pqLAXHA$lag=(0:(length(AXHA_list)-1))
  names(pqLAXHA)[names(pqLAXHA)=="AXHA_list"]="corr"
  return(pqLAXHA)
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
for (i in 1:nrow(species_info)){
  file_name=paste(path,"Data/",species_info$file_name[i],sep="")
  print(file_name)
  #Read in data
  data=read.csv(paste(file_name,".csv",sep=""))
  
  #Subset with same fishing gear
  data=subset(data, Gear == 'GOV')
  
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

  # # check on map
  # shapefile_df <- fortify(ssc_map_wgs84)
  # 
  # map <- ggplot() +
  #   geom_path(data = shapefile_df,
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



# ####For loop - Plot trend in biomass (kg) ---- 
# for (i in 14:16){
#   print(paste(i,species_info$file_name[i],sep=","))
#   # Read in data
#   data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[i]),sep=""))
#   
#   # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
#   data$ID = cumsum(!duplicated(data[3:7]))
#   
#   #Aggregate weight over each haul (sum)
#   data_agg=aggregate(Total_wgt ~ ID+Year, data=data, sum)
#   
#   #set outlier weight limit (95th percentile)
#   nf_p=quantile(data_agg$Total_wgt,0.95) 
#   
#   #subset of data excluding outlier weights
#   data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
#   
#   # convert weights from grams to kg
#   data_agg$biomass_kg=data_agg$Total_wgt/1000 
#   
#   #mean per year to get better overview
#   data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
#   
#   #Plot
#   p=ggplot(data_agg2,aes(Year,biomass_kg))+
#     geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
#     geom_point(size=3, alpha=0.5)+
#     theme_classic()+
#     ggtitle(species_info$Scientific.name[i])+
#     labs(x = "Year", y = 'Mean CPUE (kg)')+
#     theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
#     theme(axis.text.x = element_text(face="bold",size=18))+
#     theme(axis.text.y = element_text(face="bold",size=18))+
#     theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
#     theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
#     theme(axis.line = element_line(size=1))
#   p
#   
#   #Add species image
#   img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
#   g <- rasterGrob(img, interpolate=TRUE)
#   
#   p2=p+annotation_custom(g,xmin = 1977, xmax = 1986, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
#   ggsave(p2, file = paste(path,sprintf("graphs/Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
# }
# 
# 
# ####Separate plots for all data (mix Q1,Q3), clean Q1: 1977-2018 & clean Q1-Q3: 1991-2018 + set up indexes
# for (i in 1:12){
#   print(paste(i,species_info$file_name[i],sep=","))
#   # Read in data
#   data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[i]),sep=""))
#   
#   # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
#   data$ID = cumsum(!duplicated(data[3:7]))
#  
#   #All data----
#   #Aggregate weight over each haul (sum)
#   data_agg=aggregate(Total_wgt ~ ID+Year, data=data, sum)
#   
#   #set outlier weight limit (95th percentile)
#   nf_p=quantile(data_agg$Total_wgt,0.95) 
#   
#   #subset of data excluding outlier weights
#   data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
#   
#   # convert weights from grams to kg
#   data_agg$biomass_kg=data_agg$Total_wgt/1000 
#   
#   #mean per year to get better overview
#   data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
#   
#   
#   #Plot actual biomass
#   p=ggplot(data_agg2,aes(Year,biomass_kg))+
#     geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
#     geom_point(size=3, alpha=0.5)+
#     theme_classic()+
#     ggtitle(species_info$Scientific.name[i])+
#     labs(x = "Year", y = 'Mean CPUE (kg)')+
#     theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
#     theme(axis.text.x = element_text(face="bold",size=18))+
#     theme(axis.text.y = element_text(face="bold",size=18))+
#     theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
#     theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
#     theme(axis.line = element_line(size=1))
#   p
#   
#   #Add species image
#   img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
#   g <- rasterGrob(img, interpolate=TRUE)
#   
#   p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$biomass_kg)-0.15*max(data_agg2$biomass_kg),ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
#   ggsave(p2, file = paste(path,sprintf("graphs/All/Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
#   
#   #index
#   #add 1 to each biomass value to exclude potential problems with zero starting values
#   data_agg2$biomass_kg=data_agg2$biomass_kg+1
#   base = data_agg2$biomass_kg[1]
#   print(data_agg2$biomass_kg[1])
#   print(base)
#   
#   data_agg2$indexed_biomass = (data_agg2$biomass_kg/base)*100
#   
#   #Plot indexed biomass
#   p=ggplot(data_agg2,aes(Year,indexed_biomass))+
#     geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
#     geom_point(size=3, alpha=0.5)+
#     theme_classic()+
#     ggtitle(species_info$Scientific.name[i])+
#     labs(x = "Year", y = 'Indexed biomass')+
#     theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
#     theme(axis.text.x = element_text(face="bold",size=18))+
#     theme(axis.text.y = element_text(face="bold",size=18))+
#     theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
#     theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
#     theme(axis.line = element_line(size=1))
#   p
#   
#   #Add species image
#   img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
#   g <- rasterGrob(img, interpolate=TRUE)
#   
#   #different image size for callionymus lyra or will not display proportionally
#   if(species_info$file_name[i]!= "Callionymus_lyra"){
#   p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.15*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
#   ggsave(p2, file = paste(path,sprintf("graphs/All/Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
#   } else{
#     p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.25*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
#     ggsave(p2, file = paste(path,sprintf("graphs/All/Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
#     
#   }
#   
#   #Q1 (1977-2018) ---- 
#   Q1_data = subset(data,Quarter==1)
#   write.csv(Q1_data,paste(path,sprintf("Q1_filtered_%s.csv",species_info$file_name[i]),sep=""))
#   #Aggregate weight over each haul (sum)
#   data_agg=aggregate(Total_wgt ~ ID+Year, data=Q1_data, sum)
#   
#   #set outlier weight limit (95th percentile)
#   nf_p=quantile(data_agg$Total_wgt,0.95) 
#   
#   #subset of data excluding outlier weights
#   data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
#   
#   # convert weights from grams to kg
#   data_agg$biomass_kg=data_agg$Total_wgt/1000 
#   
#   #mean per year to get better overview
#   data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
#   
#   #Plot
#   p=ggplot(data_agg2,aes(Year,biomass_kg))+
#     geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
#     geom_point(size=3, alpha=0.5)+
#     theme_classic()+
#     ggtitle(species_info$Scientific.name[i])+
#     labs(x = "Year", y = 'Mean CPUE (kg)')+
#     theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
#     theme(axis.text.x = element_text(face="bold",size=18))+
#     theme(axis.text.y = element_text(face="bold",size=18))+
#     theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
#     theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
#     theme(axis.line = element_line(size=1))
#   p
#   
#   #Add species image
#   img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
#   g <- rasterGrob(img, interpolate=TRUE)
#   
#   p2=p+annotation_custom(g,xmin = 1977, xmax = 1986, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
#   ggsave(p2, file = paste(path,sprintf("graphs/Q1/Q1_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
#   
#   #index
#   #add 1 to each biomass value to exclude potential problems with zero starting values
#   data_agg2$biomass_kg=data_agg2$biomass_kg+1
#   base = data_agg2$biomass_kg[1]
#   print(data_agg2$biomass_kg[1])
#   print(base)
#   
#   data_agg2$indexed_biomass = (data_agg2$biomass_kg/base)*100
#   
#   #Plot indexed biomass
#   p=ggplot(data_agg2,aes(Year,indexed_biomass))+
#     geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
#     geom_point(size=3, alpha=0.5)+
#     theme_classic()+
#     ggtitle(species_info$Scientific.name[i])+
#     labs(x = "Year", y = 'Indexed biomass')+
#     theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
#     theme(axis.text.x = element_text(face="bold",size=18))+
#     theme(axis.text.y = element_text(face="bold",size=18))+
#     theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
#     theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
#     theme(axis.line = element_line(size=1))
#   p
#   
#   #Add species image
#   img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
#   g <- rasterGrob(img, interpolate=TRUE)
#   
#   #different image size for callionymus lyra or will not display proportionally
#   if(species_info$file_name[i]!= "Callionymus_lyra"){
#     p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.15*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
#     ggsave(p2, file = paste(path,sprintf("graphs/Q1/Q1_Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
#   } else{
#     p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.25*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
#     ggsave(p2, file = paste(path,sprintf("graphs/Q1/Q1_Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
#     
#   }
#   
#   
#   #Q1-Q3 (1991-2018)----
#   Q1_Q3_data = subset(data,Year>=1991)
#   write.csv(Q1_Q3_data,paste(path,sprintf("Q1_Q3_filtered_%s.csv",species_info$file_name[i]),sep=""))
#   Q1_Q3_data$QY = cumsum(!duplicated(Q1_Q3_data[3:4]))
#   
#   #Aggregate weight over each haul (sum)
#   data_agg=aggregate(Total_wgt ~ ID+Year+QY, data=Q1_Q3_data, sum)
#   
#   #set outlier weight limit (95th percentile)
#   nf_p=quantile(data_agg$Total_wgt,0.95) 
#   
#   #subset of data excluding outlier weights
#   data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
#   
#   # convert weights from grams to kg
#   data_agg$biomass_kg=data_agg$Total_wgt/1000 
#   
#   #mean per year to get better overview
#   data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
#   
#   #Plot
#   p=ggplot(data_agg2,aes(Year,biomass_kg))+
#     geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
#     geom_point(size=3, alpha=0.5)+
#     theme_classic()+
#     ggtitle(species_info$Scientific.name[i])+
#     labs(x = "Year", y = 'Mean CPUE (kg)')+
#     theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
#     theme(axis.text.x = element_text(face="bold",size=18))+
#     theme(axis.text.y = element_text(face="bold",size=18))+
#     theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
#     theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
#     theme(axis.line = element_line(size=1))
#   p
#   
#   #Add species image
#   img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
#   g <- rasterGrob(img, interpolate=TRUE)
#   
#   p2=p+annotation_custom(g,xmin = 1977, xmax = 1986, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
#   ggsave(p2, file = paste(path,sprintf("graphs/Q1_Q3/Q1_Q3_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
# 
#   #index
#   #add 1 to each biomass value to exclude potential problems with zero starting values
#   data_agg2$biomass_kg=data_agg2$biomass_kg+1
#   base = data_agg2$biomass_kg[1]
#   print(data_agg2$biomass_kg[1])
#   print(base)
#   
#   data_agg2$indexed_biomass = (data_agg2$biomass_kg/base)*100
#   
#   #Plot indexed biomass
#   p=ggplot(data_agg2,aes(Year,indexed_biomass))+
#     geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
#     geom_point(size=3, alpha=0.5)+
#     theme_classic()+
#     ggtitle(species_info$Scientific.name[i])+
#     labs(x = "Year", y = 'Indexed biomass')+
#     theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
#     theme(axis.text.x = element_text(face="bold",size=18))+
#     theme(axis.text.y = element_text(face="bold",size=18))+
#     theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
#     theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
#     theme(axis.line = element_line(size=1))
#   p
#   
#   #Add species image
#   img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
#   g <- rasterGrob(img, interpolate=TRUE)
#   
#   #different image size for callionymus lyra or will not display proportionally
#   if(species_info$file_name[i]!= "Callionymus_lyra"){
#     p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.15*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
#     ggsave(p2, file = paste(path,sprintf("graphs/Q1_Q3/Q1_Q3_Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
#   } else{
#     p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.25*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
#     ggsave(p2, file = paste(path,sprintf("graphs/Q1_Q3/Q1_Q3_Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
#   }
#   
# }
# 
# 
# 
# 
# ####Separate plots for all 9 seascapes all data ----
# for (i in 16:16){
#   print(paste(i,species_info$file_name[i],sep=","))
#   # Read in data
#   data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[i]),sep=""))
# 
#   # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
#   data$ID = cumsum(!duplicated(data[3:7]))
#   
#   for(j in 1:9){
#       print(paste("Seascape",j,sep=" "))
#       ssc_data=subset(data,Seascapenr==j)
#     
#       data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr, data=ssc_data, sum)
#       
#       #set outlier weight limit (95th percentile)
#       nf_p=quantile(data_agg$Total_wgt,0.95) 
#       
#       #subset of data excluding outlier weights
#       #data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
#       
#       if (nf_p > 0){
#         #subset of data excluding outlier weights
#         data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
#       }
#       
#       
#       # convert weights from grams to kg
#       data_agg$biomass_kg=data_agg$Total_wgt/1000 
#       
#       #mean per year to get better overview
#       data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
#       
#       #Plot
#       p=ggplot(data_agg2,aes(Year,biomass_kg))+
#         geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
#         geom_point(size=3, alpha=0.5)+
#         theme_classic()+
#         ggtitle(sprintf("%s Seascape %s",species_info$Scientific.name[i],j))+
#         labs(x = "Year", y = 'Mean CPUE (kg)')+
#         theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
#         theme(axis.text.x = element_text(face="bold",size=18))+
#         theme(axis.text.y = element_text(face="bold",size=18))+
#         theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
#         theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
#         theme(axis.line = element_line(size=1))
#       
#       #Add species image
#       img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
#       g <- rasterGrob(img, interpolate=TRUE)
#       
#       p2=p+annotation_custom(g,xmin = min(data_agg2$Year), xmax = min(data_agg2$Year)+9, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
#       ggsave(p2, file = paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],j,species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
#     
#     }
#   
#     #stack in single image
#     plot1 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],1,species_info$file_name[i]),sep=""))
#     plot2 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],2,species_info$file_name[i]),sep=""))
#     plot3 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],3,species_info$file_name[i]),sep=""))
#     plot4 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],4,species_info$file_name[i]),sep=""))
#     plot5 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],5,species_info$file_name[i]),sep=""))
#     plot6 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],6,species_info$file_name[i]),sep=""))
#     plot7 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],7,species_info$file_name[i]),sep=""))
#     plot8 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],8,species_info$file_name[i]),sep=""))
#     plot9 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],9,species_info$file_name[i]),sep=""))
#     
#     cplot=grid.arrange(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot5),
#                  rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),rasterGrob(plot9),ncol=3)
#     ggsave(cplot, file = paste(path,sprintf("graphs/Seascapes/All/%s/compiled_Trend_%s.png",species_info$file_name[i],species_info$file_name[i]),sep=""), height=21, width=27,dpi = 600) # adjust dpi )
# 
# }
# 
# 
# 
# 
# #############GAM INDEX VERSION (single averaged datapoint/year)############## -----
# #############################################################################
# for (j in 16:16){
#   data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[j]),sep=""))
#   print(species_info$file_name[j])
#   
#   # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
#   data$ID = cumsum(!duplicated(data[3:7]))
#   
#   ssc_data=completeFun(data,"Seascapenr")
#   ssc_data=subset(ssc_data,Quarter==1)
#   data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species, data=ssc_data, sum)
#   
#   #set outlier weight limit (95th percentile)
#   nf_p=quantile(data_agg$Total_wgt,0.95) 
#   if (nf_p > 0){
#     #subset of data excluding outlier weights
#     data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
#   }
#   
#   # convert weights from grams to kg
#   data_agg$biomass_kg=data_agg$Total_wgt/1000 
#   
#   #mean per year and seascape to get better overview
#   data_agg_year=aggregate(biomass_kg~Year+Seascapenr+Species,data=data_agg,mean)
#   
#   #separate subset per seascape, index it and finally rbind back into single dataset----
#   #1----
#   data_agg1=subset(data_agg_year,Seascapenr==1)
#   data_agg1$biomass_kg=data_agg1$biomass_kg+100 #scale up for index
#   base = data_agg1$biomass_kg[1]
#   print(data_agg1$biomass_kg[1])
#   print(base)
#   data_agg1$indexed_biomass = (data_agg1$biomass_kg/base)*100
#   data_agg1$biomass_kg=data_agg1$biomass_kg-100 #return back for actual value
#   
#   #2-----
#   data_agg2=subset(data_agg_year,Seascapenr==2)
#   data_agg2$biomass_kg=data_agg2$biomass_kg+100
#   base = data_agg2$biomass_kg[1]
#   print(data_agg2$biomass_kg[1])
#   print(base)
#   data_agg2$indexed_biomass = (data_agg2$biomass_kg/base)*100
#   data_agg2$biomass_kg=data_agg2$biomass_kg-100 #return back for actual value
#   
#   #3----
#   data_agg3=subset(data_agg_year,Seascapenr==3)
#   data_agg3$biomass_kg=data_agg3$biomass_kg+100
#   base = data_agg3$biomass_kg[1]
#   print(data_agg3$biomass_kg[1])
#   print(base)
#   data_agg3$indexed_biomass = (data_agg3$biomass_kg/base)*100
#   data_agg3$biomass_kg=data_agg3$biomass_kg-100 #return back for actual value
#   
#   #4----
#   data_agg4=subset(data_agg_year,Seascapenr==4)
#   data_agg4$biomass_kg=data_agg4$biomass_kg+100
#   base = data_agg4$biomass_kg[1]
#   print(data_agg4$biomass_kg[1])
#   print(base)
#   data_agg4$indexed_biomass = (data_agg4$biomass_kg/base)*100
#   data_agg4$biomass_kg=data_agg4$biomass_kg-100 #return back for actual value
#   
#   #5----
#   data_agg5=subset(data_agg_year,Seascapenr==5)
#   data_agg5$biomass_kg=data_agg5$biomass_kg+100
#   base = data_agg5$biomass_kg[1]
#   print(data_agg5$biomass_kg[1])
#   print(base)
#   data_agg5$indexed_biomass = (data_agg5$biomass_kg/base)*100
#   data_agg5$biomass_kg=data_agg5$biomass_kg-100 #return back for actual value
#   
#   #6----
#   data_agg6=subset(data_agg_year,Seascapenr==6)
#   data_agg6$biomass_kg=data_agg6$biomass_kg+100
#   base = data_agg6$biomass_kg[1]
#   print(data_agg6$biomass_kg[1])
#   print(base)
#   data_agg6$indexed_biomass = (data_agg6$biomass_kg/base)*100
#   data_agg6$biomass_kg=data_agg6$biomass_kg-100 #return back for actual value
#   
#   #7----
#   data_agg7=subset(data_agg_year,Seascapenr==7)
#   data_agg7$biomass_kg=data_agg7$biomass_kg+100
#   base = data_agg7$biomass_kg[1]
#   print(data_agg7$biomass_kg[1])
#   print(base)
#   data_agg7$indexed_biomass = (data_agg7$biomass_kg/base)*100
#   data_agg7$biomass_kg=data_agg7$biomass_kg-100 #return back for actual value
#   
#   #8----
#   data_agg8=subset(data_agg_year,Seascapenr==8)
#   data_agg8$biomass_kg=data_agg8$biomass_kg+100
#   base = data_agg8$biomass_kg[1]
#   print(data_agg8$biomass_kg[1])
#   data_agg8$indexed_biomass = (data_agg8$biomass_kg/base)*100
#   data_agg8$biomass_kg=data_agg8$biomass_kg-100 #return back for actual value
#   
#   #9----
#   data_agg9=subset(data_agg_year,Seascapenr==9)
#   data_agg9$biomass_kg=data_agg9$biomass_kg+100
#   base = data_agg9$biomass_kg[1]
#   print(data_agg9$biomass_kg[1])
#   print(base)
#   data_agg9$indexed_biomass = (data_agg9$biomass_kg/base)*100
#   data_agg9$biomass_kg=data_agg9$biomass_kg-100 #return back for actual value
#   
#   #rbind----
#   data_agg_indexed=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9)
#   
#   
#   ####################################################
#   #GAM model with separate smoother per seascape
#   ###################################################
#   
#   data_agg_indexed$Seascapenr=as.factor(data_agg_indexed$Seascapenr) #Set as factor
#   data_agg_indexed = start_event(data_agg_indexed, column="Year",event="Seascapenr") #Sep time series per factor level
#   
#   #m=bam(biomass_kg~Seascapenr+s(Year, by=Seascapenr),data=data_agg_indexed,method="REML")
#   #(valRho <- acf(resid(m), plot=FALSE)$acf[2])
#   #AIC(m)
#   
#   mi=bam(indexed_biomass~Seascapenr+s(Year, by=Seascapenr),data=data_agg_indexed,method="REML")
#   (valRho <- acf(resid(mi), plot=FALSE)$acf[2])
#   AIC(mi)
#   
#   
#   mi1=bam(indexed_biomass~Seascapenr+s(Year,by=Seascapenr),data=data_agg_indexed,method="REML",AR.start=data_agg_indexed$start.event,rho=valRho)
#   AIC(mi1)
#   
#   par(mfrow=c(2,2))
#   gam.check(mi1)
#   #savePlot(acf,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ACF.png",species_info$file_name[j]),sep=""))
#   dev.copy(png,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_diagnostics.png",species_info$file_name[j]),sep=""))
#   dev.off()
# 
#   
#   
#   ##############################################################  #commented out/obsolete
#   #Manual Trial difference smooths Seascape 1 - Seascape 2----
#     # pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
#     # #Control - Nurse plants
#     # xp <- predict(mi1, pdat, type = 'lpmatrix')
#     # col1 <- grepl("Seascapenr4", colnames(xp))
#     # col2 <- grepl("Seascapenr7", colnames(xp))
#     # row1 <- pdat[["Seascapenr"]] == 4
#     # row2 <- pdat[["Seascapenr"]] == 7
#     # 
#     # 
#     # ## difference rows of pred for data from comparison
#     # X <- xp[row1, ] - xp[row2, ]
#     # ## zero out cols of X related to splines for other lochs
#     # X[, ! (col1 | col2)] <- 0
#     # ## zero out the parametric cols
#     # X[, !grepl('^s\\(', colnames(xp))] <- 0
#     # 
#     # #predicted values for differences
#     # Vb=vcov(mi1)
#     # dif <- X %*% coef(mi1)
#     # se <- sqrt(rowSums((X %*% Vb) * X))
#     # pred=cbind(dif,se)
#     # 
#     # ## CI of the estimates
#     # set.seed(42)
#     # N <- 10000
#     # 
#     # #predi <- predict(m1, pdat, se.fit = TRUE)
#     # #se.fit <- predi$se.fit
#     # 
#     # BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)    
#     # Cg <- predict(mi1, pdat, type = "lpmatrix")
#     # simDev <- Cg %*% t(BUdiff)
#     # 
#     # absDev <- abs(sweep(simDev, 1, se, FUN = "/"))
#     # masd <- apply(absDev, 2L, max)
#     # crit <- quantile(masd, prob = 0.95, type = 8)  
#     # 
#     # newd <- expand.grid(Year = seq(1977, 2019, length = 126))
#     # 
#     # S12pred <- transform(cbind(data.frame(pred), newd),
#     #                      uprP = V1 + (2 * se),
#     #                      lwrP = V1 - (2 * se),
#     #                      uprS = V1 + (crit * se),
#     #                      lwrS = V1 - (crit * se))  
#     # S12pred$sigdif <-ifelse(S12pred$lwrS<0 & S12pred$uprS<0,"lower",ifelse(S12pred$lwrS>0 & S12pred$uprS>0,"higher","no"))
#     # S12higher=S12pred
#     # S12higher["V1"][S12pred["sigdif"]!="higher"]=NA
#     # S12lower=S12pred
#     # S12lower["V1"][S12pred["sigdif"]!="lower"]=NA
#     # 
#     # #S12higher[S12higher$sigdif!="higher",]=NA
#     # 
#     # #S12higher=S12pred[which(S12pred$sigdif=="higher"),]
#     # #S12pred$higher=S12pred$sigdif
#     # #S12pred$higher[S12pred$higher!="higher"]=NA
#     # #S12pred$lower=S12pred$sigdif
#     # #S12pred$lower[S12pred$lower!="lower"]=NA
#     # 
#     # #mean(S12higher$V1)
#     # #sd(S12higher$V1)
#     # #S12lower=as.data.frame(S12pred$lower)
#     # #S12higher=as.data.frame(S12pred$higher)
#     # #mean(HCNlower$V1)
#     # #sd(HCNlower$V1)
#     # S12pred$sigdif=factor(levels="lower","higher","no")
#     # 
#     # DifS12plot=
#     #   ggplot(S12pred, aes(x = Year, y = V1)) +
#     #   theme_classic()+
#     #   #geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
#     #   geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.1) +
#     #   geom_line(data=S12pred,aes(x=Year,y=V1),col="black",size=1) +
#     #   geom_line(data=S12higher,aes(x=Year,y=V1),col="springgreen1",size=2)+
#     #   geom_line(data=S12lower,aes(x=Year,y=V1),col="red",size=2)+
#     #   ggtitle("Seascape 1 - 2")+
#     #   labs(x = "Year", y = 'Difference Fish biomass index')+
#     #   theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
#     #   theme(axis.text.x = element_text(face="bold",size=18))+
#     #   theme(axis.text.y = element_text(face="bold",size=18))+
#     #   theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
#     #   theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
#     #   theme(axis.line = element_line(size=1))
#     # 
#     # DifS12plot
#   
#   
#   ##############################################################
#   
#    Vb <- vcov(mi1)
#    newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
#    pred <- predict(mi1, newd, se.fit = TRUE)
#    se.fit <- pred$se.fit
#   # 
#    set.seed(42)
#    N <- 10000
#   # 
#    BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
#   # 
#   # xp matrix where basisfunctions of the model have been evaluated at 126 time point values per area
#    Cg <- predict(mi1, newd, type = "lpmatrix")
#    simDev <- Cg %*% t(BUdiff)
#   # # 
#    absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
#    masd <- apply(absDev, 2L, max)
#    crit <- quantile(masd, prob = 0.95, type = 8)
#    pred <- transform(cbind(data.frame(pred), newd),
#                      uprP = fit + (2 * se.fit),
#                      lwrP = fit - (2 * se.fit),
#                      uprS = fit + (crit * se.fit),
#                      lwrS = fit - (crit * se.fit))
#   # # 
#   # # 
#   #separate plots per area
#   S1=pred[which(pred$Seascapenr==1),]
#   S2=pred[which(pred$Seascapenr==2),]
#   S3=pred[which(pred$Seascapenr==3),]
#   S4=pred[which(pred$Seascapenr==4),]
#   S5=pred[which(pred$Seascapenr==5),]
#   S6=pred[which(pred$Seascapenr==6),]
#   S7=pred[which(pred$Seascapenr==7),]
#   S8=pred[which(pred$Seascapenr==8),]
#   S9=pred[which(pred$Seascapenr==9),]
#   
#   # # # 
#    s1plot=ggplot(data=S1,aes(x=Year))+
#      geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#      geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#      geom_line(aes(y=fit),col="black",size=1)+
#      ggtitle("S1")+
#      labs(x="Year",y="Indexed fish biomass trend")
#    s2plot=ggplot(data=S2,aes(x=Year))+
#      geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#      geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#      geom_line(aes(y=fit),col="black",size=1)+
#      ggtitle("S2")+
#      labs(x="Year",y="Indexed fish biomass trend")
#    s3plot=ggplot(data=S3,aes(x=Year))+
#      geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#      geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#      geom_line(aes(y=fit),col="black",size=1)+
#      ggtitle("S3")+
#      labs(x="Year",y="Indexed fish biomass trend")
#    s4plot=ggplot(data=S4,aes(x=Year))+
#      geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#      geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#      geom_line(aes(y=fit),col="black",size=1)+
#      ggtitle("S4")+
#      labs(x="Year",y="Indexed fish biomass trend")
#    s5plot=ggplot(data=S5,aes(x=Year))+
#      geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#      geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#      geom_line(aes(y=fit),col="black",size=1)+
#      ggtitle("S5")+
#      labs(x="Year",y="Indexed fish biomass trend")
#    s6plot=ggplot(data=S6,aes(x=Year))+
#      geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#      geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#      geom_line(aes(y=fit),col="black",size=1)+
#      ggtitle("S6")+
#      labs(x="Year",y="Indexed fish biomass trend")
#    s7plot=ggplot(data=S7,aes(x=Year))+
#      geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#      geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#      geom_line(aes(y=fit),col="black",size=1)+
#      ggtitle("S7")+
#      labs(x="Year",y="Indexed fish biomass trend")
#    s8plot=ggplot(data=S8,aes(x=Year))+
#      geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#      geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#      geom_line(aes(y=fit),col="black",size=1)+
#      ggtitle("S8")+
#      labs(x="Year",y="Indexed fish biomass trend")
#    s9plot=ggplot(data=S9,aes(x=Year))+
#      geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#      geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#      geom_line(aes(y=fit),col="black",size=1)+
#      ggtitle("S9")+
#      labs(x="Year",y="Indexed fish biomass trend")
#  
#   cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,ncol=3,nrow=3)
#   ggsave(cplot,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_smooths.png",species_info$file_name[j]),sep=""),height=7, width=9,dpi = 600)# 
#   
#   #save residual autocorrelation plot
#   par(mfrow=c(1,1))
#   acf_resid(mi1,split_pred="Seascapenr",main="ACF resid(m1)")
#   dev.copy(png,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_autocorrelation.png",species_info$file_name[j]),sep=""))
#   dev.off()
#   
#   ######################################
#   #Compute difference smooths ----
#   ######################################
#   combi=c(1,2,3,4,5,6,7,8,9)
#   combi_width=as.data.frame(combinations(n=9,r=2,v=combi))
#   pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
#   xp <- predict(mi1, pdat, type = 'lpmatrix')
#   
#   sprintf("Seascapenr%s",combi_width$V1[2])
#   sprintf("Seascapenr%s",combi_width$V2[2])
#   
#   #for (i in 1:nrow(combi_width)){
#   for (i in 1:nrow(combi_width)){
#     print(i)
#     pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
#     xp <- predict(mi1, pdat, type = 'lpmatrix')
#     col1 <- grepl(sprintf("Seascapenr%s",combi_width$V1[i]), colnames(xp))
#     col2 <- grepl(sprintf("Seascapenr%s",combi_width$V2[i]), colnames(xp))
#     row1 <- pdat[["Seascapenr"]] == combi_width$V1[i]
#     row2 <- pdat[["Seascapenr"]] == combi_width$V2[i] 
#     
#     ## difference rows of pred for data from comparison
#     X <- xp[row1, ] - xp[row2, ]
#     ## zero out cols of X related to splines for other lochs
#     X[, ! (col1 | col2)] <- 0
#     ## zero out the parametric cols
#     X[, !grepl('^s\\(', colnames(xp))] <- 0
#     
#     #predicted values for differences
#     Vb=vcov(mi1)
#     dif <- X %*% coef(mi1)
#     se <- sqrt(rowSums((X %*% Vb) * X))
#     pred=cbind(dif,se)
#     
#     ## CI of the estimates
#     set.seed(42)
#     N <- 10000
#     
#     #predi <- predict(m1, pdat, se.fit = TRUE)
#     #se.fit <- predi$se.fit
#     
#     BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)    
#     Cg <- predict(mi1, pdat, type = "lpmatrix")
#     simDev <- Cg %*% t(BUdiff)
#     
#     absDev <- abs(sweep(simDev, 1, se, FUN = "/"))
#     masd <- apply(absDev, 2L, max)
#     crit <- quantile(masd, prob = 0.95, type = 8)  
#     
#     newd <- expand.grid(Year = seq(1977, 2019, length = 126))
#     
#     S12pred <- transform(cbind(data.frame(pred), newd),
#                          uprP = V1 + (2 * se),
#                          lwrP = V1 - (2 * se),
#                          uprS = V1 + (crit * se),
#                          lwrS = V1 - (crit * se))  
#     S12pred$sigdif <-ifelse(S12pred$lwrS<0 & S12pred$uprS<0,"lower",ifelse(S12pred$lwrS>0 & S12pred$uprS>0,"higher","no"))
#     S12higher=S12pred
#     S12higher["V1"][S12pred["sigdif"]!="higher"]=NA
#     S12lower=S12pred
#     S12lower["V1"][S12pred["sigdif"]!="lower"]=NA
#     
#     DifS12plot=
#       ggplot(S12pred, aes(x = Year, y = V1)) +
#       theme_bw()+
#       #geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
#       geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
#       geom_line(data=S12pred,aes(x=Year,y=V1),col="black",size=0.5) +
#       geom_line(data=S12higher,aes(x=Year,y=V1),col="seagreen1",size=1)+
#       geom_line(data=S12lower,aes(x=Year,y=V1),col="red",size=1)+
#       ggtitle(sprintf("%s Seascape %s - %s",species_info$Scientific.name[j],combi_width$V1[i],combi_width$V2[i]))+
#       labs(x = "Year", y = 'Difference Fish biomass index')+
#       theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
#       theme(axis.text.x = element_text(face="bold",size=9))+
#       theme(axis.text.y = element_text(face="bold",size=9))+
#       theme(axis.title.x = element_text(face="bold",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
#       theme(axis.title.y = element_text(face="bold",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
#       theme(axis.line = element_line(size=1))
#     ggsave(DifS12plot,file = paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooth_%s_%s.png",species_info$file_name[j],combi_width$V1[i],combi_width$V2[i]),sep=""), height=3.5, width=4.5,dpi = 600)
# 
#   }
# }  
# 
# #save different smooths in single large image            
# for (j in 16:16){
#   print(j)
#   p=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth",species_info$file_name[j]),sep="")
#   filenames=list.files(path=p,pattern="*.png",full.names=TRUE)
#   
#   rl <- lapply(filenames, png::readPNG)
#   gl <- lapply(rl, grid::rasterGrob)
#   cp=grid.arrange(gl[[1]],gl[[2]],gl[[3]],gl[[4]],gl[[5]],gl[[6]],
#                   gl[[7]],gl[[8]],gl[[9]],gl[[10]],gl[[11]],gl[[12]],
#                   gl[[13]],gl[[14]],gl[[15]],gl[[16]],gl[[17]],gl[[18]],
#                   gl[[19]],gl[[20]],gl[[21]],gl[[22]],gl[[23]],gl[[24]],
#                   gl[[25]],gl[[26]],gl[[27]],gl[[28]],gl[[29]],gl[[30]],
#                   gl[[31]],gl[[32]],gl[[33]],gl[[34]],gl[[35]],gl[[36]],ncol=6,nrow=6)
#   ggsave(cp,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooths.png",species_info$file_name[j]),sep=""),height=21, width=27,dpi = 300)
# }
# 
# 
# 
# 
# ################################################################################
# ####GAM Anscombe transformed (sqrt(x+3/8)) and scaled by max (BETA errror)   ###-----
# ################################################################################
# list_AIC_beta_2= vector("list",16)
# 
# 
# for (j in 1:16){
#   data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[j]),sep=""))
#   print(species_info$file_name[j])
#   
#   # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
#   data$ID = cumsum(!duplicated(data[3:7]))
#   
#   ssc_data=completeFun(data,"Seascapenr")
#   ssc_data=subset(ssc_data,Quarter==1)
#   data_agg=ssc_data
#   data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species, data=ssc_data, sum)
#   # 
#    # #set outlier weight limit (95th percentile)
#    #  nf_p=quantile(data_agg$Total_wgt,0.95)
#    #  if (nf_p > 0){
#    #   #subset of data excluding outlier weights
#    #    data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
#    #  }
#   
#   # convert weights from grams to kg
#   data_agg$biomass_kg=data_agg$Total_wgt/1000
#   
#   #separate subset per seascape, anscombe transform, max scale it and finally rbind back into single dataset----
#   #1----
#   data_agg1=subset(data_agg,Seascapenr==1)
#   data_agg1$ansc_biomass= nthroot(data_agg1$biomass_kg+(3/8),2)
#   data_agg1$scaled_ansc_biomass=data_agg1$ansc_biomass/max(data_agg1$ansc_biomass)
# 
#   #2-----
#   data_agg2=subset(data_agg,Seascapenr==2)
#   data_agg2$ansc_biomass= nthroot(data_agg2$biomass_kg+(3/8),2)
#   data_agg2$scaled_ansc_biomass=data_agg2$ansc_biomass/max(data_agg2$ansc_biomass)
#   
#   #3----
#   data_agg3=subset(data_agg,Seascapenr==3)
#   data_agg3$ansc_biomass= nthroot(data_agg3$biomass_kg+(3/8),2)
#   data_agg3$scaled_ansc_biomass=data_agg3$ansc_biomass/max(data_agg3$ansc_biomass)
#  
#   #4----
#   data_agg4=subset(data_agg,Seascapenr==4)
#   data_agg4$ansc_biomass= nthroot(data_agg4$biomass_kg+(3/8),2)
#   data_agg4$scaled_ansc_biomass=data_agg4$ansc_biomass/max(data_agg4$ansc_biomass)
# 
#   #5----
#   data_agg5=subset(data_agg,Seascapenr==5)
#   data_agg5$ansc_biomass= nthroot(data_agg5$biomass_kg+(3/8),2)
#   data_agg5$scaled_ansc_biomass=data_agg5$ansc_biomass/max(data_agg5$ansc_biomass)
#   
#   #6----
#   data_agg6=subset(data_agg,Seascapenr==6)
#   data_agg6$ansc_biomass= nthroot(data_agg6$biomass_kg+(3/8),2)
#   data_agg6$scaled_ansc_biomass=data_agg6$ansc_biomass/max(data_agg6$ansc_biomass)
#   
#   #7----
#   data_agg7=subset(data_agg,Seascapenr==7)
#   data_agg7$ansc_biomass= nthroot(data_agg7$biomass_kg+(3/8),2)
#   data_agg7$scaled_ansc_biomass=data_agg7$ansc_biomass/max(data_agg7$ansc_biomass)
#   
#   #8----
#   data_agg8=subset(data_agg,Seascapenr==8)
#   data_agg8$ansc_biomass= nthroot(data_agg8$biomass_kg+(3/8),2)
#   data_agg8$scaled_ansc_biomass=data_agg8$ansc_biomass/max(data_agg8$ansc_biomass)
#   
#   
#   #9----
#   data_agg9=subset(data_agg,Seascapenr==9)
#   data_agg9$ansc_biomass= nthroot(data_agg9$biomass_kg+(3/8),2)
#   data_agg9$scaled_ansc_biomass=data_agg9$ansc_biomass/max(data_agg9$ansc_biomass)
#   
#   #rbind----
#   data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9)
#   
#   ####################################################
#   #GAM model with separate smoother per seascape
#   ###################################################
#   print("gam")
#   data_agg_ansc$Seascapenr=as.factor(data_agg_ansc$Seascapenr) #Set as factor
#   data_agg_ansc = start_event(data_agg_ansc, column="Year",event="Seascapenr") #Sep time series per factor level
#   
#   
#   mod=bam(scaled_ansc_biomass~Seascapenr+s(Year, by=Seascapenr),data=data_agg_ansc,method="REML",family=betar)
#   print("gam fitted")
#   print(AIC(mod))
#   par(mfrow=c(2,2))
#   gam.check(mod)
#   #list_AIC_beta[[j]]=AIC(mod)
#   
#   dev.copy(png,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_diagnostics_2.png",species_info$file_name[j]),sep=""))
#   dev.off()
#   
#   
#   # ################################################################
#   # # Plot predicted smooths with all data                         #----
#   # ################################################################
#   newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
#   pred <- predict(mod, newd, se.fit = TRUE,type="response")
#   se.fit <- pred$se.fit
#   pred <- transform(cbind(data.frame(pred), newd),
#                     uprP = fit + (2 * se.fit),
#                     lwrP = fit - (2 * se.fit))
# 
#   #separate plots per area
#   S1=pred[which(pred$Seascapenr==1),]
#   S2=pred[which(pred$Seascapenr==2),]
#   S3=pred[which(pred$Seascapenr==3),]
#   S4=pred[which(pred$Seascapenr==4),]
#   S5=pred[which(pred$Seascapenr==5),]
#   S6=pred[which(pred$Seascapenr==6),]
#   S7=pred[which(pred$Seascapenr==7),]
#   S8=pred[which(pred$Seascapenr==8),]
#   S9=pred[which(pred$Seascapenr==9),]
# 
# 
#   img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[j]),sep=""))
#   g <- rasterGrob(img, interpolate=TRUE)
#   
#   ### plot smooths with all data -----  
#   s1plot=ggplot(data=S1,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg1)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S1")+
#     labs(x="Year",y="ANSC fish biomass trend")
#   #Add species image
#   s1plot=s1plot+annotation_custom(g,xmin = min(data_agg1$Year), xmax = min(data_agg1$Year)+9, ymin = max(data_agg1$scaled_ansc_biomass)-0.15*(max(data_agg1$scaled_ansc_biomass)), ymax = max(data_agg1$scaled_ansc_biomass))
# 
#   
#   s2plot=ggplot(data=S2,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg2)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S2")+
#     labs(x="Year",y="ANSC fish biomass trend")
#   #Add species image
#   #s2plot=s2plot+annotation_custom(g,xmin = min(data_agg2$Year), xmax = min(data_agg2$Year)+8, ymin = max(data_agg2$scaled_ansc_biomass)-0.15*(max(data_agg2$scaled_ansc_biomass)), ymax = max(data_agg2$scaled_ansc_biomass))
#   
#   
#   s3plot=ggplot(data=S3,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg3)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S3")+
#     labs(x="Year",y="ANSC fish biomass trend")
#   #Add species image
#   #s3plot=s3plot+annotation_custom(g,xmin = min(data_agg3$Year), xmax = min(data_agg3$Year)+8, ymin = max(data_agg3$scaled_ansc_biomass)-0.15*(max(data_agg3$scaled_ansc_biomass)), ymax = max(data_agg3$scaled_ansc_biomass))
#   
#   s4plot=ggplot(data=S4,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg4)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S4")+
#     labs(x="Year",y="ANSC fish biomass trend")
#   #Add species image
#   #s4plot=s4plot+annotation_custom(g,xmin = min(data_agg4$Year), xmax = min(data_agg4$Year)+8, ymin = max(data_agg4$scaled_ansc_biomass)-0.15*(max(data_agg4$scaled_ansc_biomass)), ymax = max(data_agg4$scaled_ansc_biomass))
#   
#   
#   s5plot=ggplot(data=S5,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg5)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S5")+
#     labs(x="Year",y="ANSC fish biomass trend")
#   #Add species image
#   #s5plot=s5plot+annotation_custom(g,xmin = min(data_agg5$Year), xmax = min(data_agg5$Year)+8, ymin = max(data_agg5$scaled_ansc_biomass)-0.15*(max(data_agg5$scaled_ansc_biomass)), ymax = max(data_agg5$scaled_ansc_biomass))
# 
#   
#   s6plot=ggplot(data=S6,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg6)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S6")+
#     labs(x="Year",y="ANSC fish biomass trend")
#   #Add species image
#   #s6plot=s6plot+annotation_custom(g,xmin = min(data_agg6$Year), xmax = min(data_agg6$Year)+8, ymin = max(data_agg6$scaled_ansc_biomass)-0.15*(max(data_agg6$scaled_ansc_biomass)), ymax = max(data_agg6$scaled_ansc_biomass))
#   
#   
#   s7plot=ggplot(data=S7,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg7)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S7")+
#     labs(x="Year",y="ANSC fish biomass trend")
#   #Add species image
#   #s7plot=s7plot+annotation_custom(g,xmin = min(data_agg7$Year), xmax = min(data_agg7$Year)+8, ymin = max(data_agg7$scaled_ansc_biomass)-0.15*(max(data_agg7$scaled_ansc_biomass)), ymax = max(data_agg7$scaled_ansc_biomass))
#   
#   s8plot=ggplot(data=S8,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg8)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S8")+
#     labs(x="Year",y="ANSC fish biomass trend")
#   #Add species image
#   #s8plot=s8plot+annotation_custom(g,xmin = min(data_agg8$Year), xmax = min(data_agg8$Year)+8, ymin = max(data_agg8$scaled_ansc_biomass)-0.15*(max(data_agg8$scaled_ansc_biomass)), ymax = max(data_agg8$scaled_ansc_biomass))
#   
#   
#   s9plot=ggplot(data=S9,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg9)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S9")+
#     labs(x="Year",y="ANSC fish biomass trend")
#   #Add species image
#   #s9plot=s9plot+annotation_custom(g,xmin = min(data_agg9$Year), xmax = min(data_agg9$Year)+8, ymin = max(data_agg9$scaled_ansc_biomass)-0.15*(max(data_agg9$scaled_ansc_biomass)), ymax = max(data_agg9$scaled_ansc_biomass))
#   
#   
#   cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,ncol=3,nrow=3)
#   ggsave(cplot,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_data_smooths.png",species_info$file_name[j]),sep=""),height=7, width=9,dpi = 600)#
# 
# 
#   ################################################################
#   # Plot predicted smooths at smooth scale and simultaneous CI   #----
#   ################################################################
#   newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
#   pred <- predict(mod, newd, se.fit = TRUE,type="response")
#   se.fit <- pred$se.fit
#   pred <- transform(cbind(data.frame(pred), newd),
#                     uprP = fit + (2 * se.fit),
#                     lwrP = fit - (2 * se.fit))
#   ###
#   Vb <- vcov(mod)
#   newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
#   pred <- predict(mod, newd, se.fit = TRUE)
#   se.fit <- pred$se.fit
# 
#   set.seed(42)
#   N <- 10000
# 
#   BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
# 
#   #xp matrix where basisfunctions of the model have been evaluated at 400 time point values per area
#   Cg <- predict(mod, newd, type = "lpmatrix")
#   simDev <- Cg %*% t(BUdiff)
# 
#   absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
#   masd <- apply(absDev, 2L, max)
#   crit <- quantile(masd, prob = 0.95, type = 8)
#   pred <- transform(cbind(data.frame(pred), newd),
#                     uprP = fit + (2 * se.fit),
#                     lwrP = fit - (2 * se.fit),
#                     uprS = fit + (crit * se.fit),
#                     lwrS = fit - (crit * se.fit))
#   #separate plots per area
#   S1=pred[which(pred$Seascapenr==1),]
#   S2=pred[which(pred$Seascapenr==2),]
#   S3=pred[which(pred$Seascapenr==3),]
#   S4=pred[which(pred$Seascapenr==4),]
#   S5=pred[which(pred$Seascapenr==5),]
#   S6=pred[which(pred$Seascapenr==6),]
#   S7=pred[which(pred$Seascapenr==7),]
#   S8=pred[which(pred$Seascapenr==8),]
#   S9=pred[which(pred$Seascapenr==9),]
# 
#   ### plot smooths with all data -----
#   s1plot=ggplot(data=S1,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S1")+
#     labs(x="Year",y="ANSC fish biomass trend")
# 
#   s2plot=ggplot(data=S2,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S2")+
#     labs(x="Year",y="ANSC fish biomass trend")
# 
#   s3plot=ggplot(data=S3,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S3")+
#     labs(x="Year",y="ANSC fish biomass trend")
# 
#   s4plot=ggplot(data=S4,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S4")+
#     labs(x="Year",y="ANSC fish biomass trend")
# 
#   s5plot=ggplot(data=S5,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S5")+
#     labs(x="Year",y="ANSC fish biomass trend")
# 
#   s6plot=ggplot(data=S6,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S6")+
#     labs(x="Year",y="ANSC fish biomass trend")
# 
#   s7plot=ggplot(data=S7,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S7")+
#     labs(x="Year",y="ANSC fish biomass trend")
#   s8plot=ggplot(data=S8,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S8")+
#     labs(x="Year",y="ANSC fish biomass trend")
# 
#   s9plot=ggplot(data=S9,aes(x=Year))+
#     theme_bw()+
#     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
#     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
#     geom_line(aes(y=fit),col="black",size=1)+
#     ggtitle("S9")+
#     labs(x="Year",y="ANSC fish biomass trend")
# 
#   cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,ncol=3,nrow=3)
#   ggsave(cplot,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_smooths.png",species_info$file_name[j]),sep=""),height=7, width=9,dpi = 600)#
# 
#   #save residual autocorrelation plot
#   par(mfrow=c(1,2))
#   acf(residuals(mod),main="ACF(mod)")
#   pacf(residuals(mod),main="partial ACF(mod)")
#   dev.copy(png,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_autocorrelation.png",species_info$file_name[j]),sep=""))
#   dev.off()
# 
#   #####################################################
#   # Compute difference smooths using simultaneous CI's #  ----
#   #####################################################
#   combi=c(1,2,3,4,5,6,7,8,9)
#   combi_width=as.data.frame(combinations(n=9,r=2,v=combi))
#   pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
#   xp <- predict(mod, pdat, type = 'lpmatrix')
# 
#   sprintf("Seascapenr%s",combi_width$V1[2])
#   sprintf("Seascapenr%s",combi_width$V2[2])
# 
#   for (i in 1:nrow(combi_width)){
#     print(i)
#     pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
#     xp <- predict(mod, pdat, type = 'lpmatrix')
#     col1 <- grepl(sprintf("Seascapenr%s",combi_width$V1[i]), colnames(xp))
#     col2 <- grepl(sprintf("Seascapenr%s",combi_width$V2[i]), colnames(xp))
#     row1 <- pdat[["Seascapenr"]] == combi_width$V1[i]
#     row2 <- pdat[["Seascapenr"]] == combi_width$V2[i]
# 
#     ## difference rows of pred for data from comparison
#     X <- xp[row1, ] - xp[row2, ]
#     ## zero out cols of X related to splines for other lochs
#     X[, ! (col1 | col2)] <- 0
#     ## zero out the parametric cols
#     X[, !grepl('^s\\(', colnames(xp))] <- 0
# 
#     #predicted values for differences
#     Vb=vcov(mod)
#     dif <- X %*% coef(mod)
#     se <- sqrt(rowSums((X %*% Vb) * X))
#     pred=cbind(dif,se)
# 
#     ## CI of the estimates
#     set.seed(42)
#     N <- 10000
# 
#     #predi <- predict(m1, pdat, se.fit = TRUE)
#     #se.fit <- predi$se.fit
# 
#     BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
#     Cg <- predict(mod, pdat, type = "lpmatrix")
#     simDev <- Cg %*% t(BUdiff)
# 
#     absDev <- abs(sweep(simDev, 1, se, FUN = "/"))
#     masd <- apply(absDev, 2L, max)
#     crit <- quantile(masd, prob = 0.95, type = 8)
# 
#     newd <- expand.grid(Year = seq(1977, 2019, length = 126))
# 
#     S12pred <- transform(cbind(data.frame(pred), newd),
#                          uprP = V1 + (2 * se),
#                          lwrP = V1 - (2 * se),
#                          uprS = V1 + (crit * se),
#                          lwrS = V1 - (crit * se))
#     S12pred$sigdif <-ifelse(S12pred$lwrS<0 & S12pred$uprS<0,"lower",ifelse(S12pred$lwrS>0 & S12pred$uprS>0,"higher","no"))
#     S12higher=S12pred
#     S12higher["V1"][S12pred["sigdif"]!="higher"]=NA
#     S12lower=S12pred
#     S12lower["V1"][S12pred["sigdif"]!="lower"]=NA
# 
#     DifS12plot=
#       ggplot(S12pred, aes(x = Year, y = V1)) +
#       theme_bw()+
#       #geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
#       geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
#       geom_line(data=S12pred,aes(x=Year,y=V1),col="black",size=0.5) +
#       geom_line(data=S12higher,aes(x=Year,y=V1),col="seagreen1",size=1)+
#       geom_line(data=S12lower,aes(x=Year,y=V1),col="red",size=1)+
#       ggtitle(sprintf("%s Seascape %s - %s",species_info$Scientific.name[j],combi_width$V1[i],combi_width$V2[i]))+
#       labs(x = "Year", y = 'Difference Fish biomass trend')+
#       theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
#       theme(axis.text.x = element_text(face="bold",size=9))+
#       theme(axis.text.y = element_text(face="bold",size=9))+
#       theme(axis.title.x = element_text(face="bold",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
#       theme(axis.title.y = element_text(face="bold",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
#       theme(axis.line = element_line(size=1))
#     ggsave(DifS12plot,file = paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooth_%s_%s_ANSC.png",species_info$file_name[j],combi_width$V1[i],combi_width$V2[i]),sep=""), height=3.5, width=4.5,dpi = 600)
# 
#   }
# }
# 
# 
# #save different smooths in single large image            
# for (j in 1:1){
#   print(j)
#   p=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth",species_info$file_name[j]),sep="")
#   filenames=list.files(path=p,pattern="*ANSC.png",full.names=TRUE)
#   
#   rl <- lapply(filenames, png::readPNG)
#   gl <- lapply(rl, grid::rasterGrob)
#   cp=grid.arrange(gl[[1]],gl[[2]],gl[[3]],gl[[4]],gl[[5]],gl[[6]],
#                   gl[[7]],gl[[8]],gl[[9]],gl[[10]],gl[[11]],gl[[12]],
#                   gl[[13]],gl[[14]],gl[[15]],gl[[16]],gl[[17]],gl[[18]],
#                   gl[[19]],gl[[20]],gl[[21]],gl[[22]],gl[[23]],gl[[24]],
#                   gl[[25]],gl[[26]],gl[[27]],gl[[28]],gl[[29]],gl[[30]],
#                   gl[[31]],gl[[32]],gl[[33]],gl[[34]],gl[[35]],gl[[36]],ncol=6,nrow=6)
#   ggsave(cp,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooths_ANSC.png",species_info$file_name[j]),sep=""),height=21, width=27,dpi = 300)
# }




################################################################################
####GAM Anscombe transformed (sqrt(x+3/8)) and scaled by max (BETA errror)   ###-----
################################################################################
list_AIC_beta_2= vector("list",16)


for (j in 12:16){
  data=read.csv(paste(path,sprintf("Data/Filtered_%s.csv",species_info$file_name[j]),sep=""))
  names(data)[names(data) == "id"] <- "Seascapenr"
  print(species_info$file_name[j])
  
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  ssc_data=completeFun(data,"Seascapenr")
  ssc_data=subset(ssc_data,Quarter==1)
  data_agg=ssc_data
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species, data=ssc_data, sum)

  # # #set outlier weight limit (99th percentile)
  #   nf_p=quantile(data_agg$Total_wgt,0.99)
  #   if (nf_p > 0){
  #    #subset of data excluding outlier weights
  #     data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
  #   }
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000
  
  #separate subset per seascape, anscombe transform, max scale it and finally rbind back into single dataset----
  #1----
  data_agg1=subset(data_agg,Seascapenr==1)
  data_agg1$ansc_biomass= nthroot(data_agg1$biomass_kg+(3/8),2)
  data_agg1$scaled_ansc_biomass=data_agg1$ansc_biomass/max(data_agg1$ansc_biomass)
  data_agg1$scaled_ansc_biomass=scaleTR(data_agg1$ansc_biomass)
  data_agg1$scaled_ansc_biomass[is.na(data_agg1$scaled_ansc_biomass)] <- 0.0001
  
  #2-----
  data_agg2=subset(data_agg,Seascapenr==2)
  data_agg2$ansc_biomass= nthroot(data_agg2$biomass_kg+(3/8),2)
  data_agg2$scaled_ansc_biomass=data_agg2$ansc_biomass/max(data_agg2$ansc_biomass)
  data_agg2$scaled_ansc_biomass=scaleTR(data_agg2$ansc_biomass)
  data_agg2$scaled_ansc_biomass[is.na(data_agg2$scaled_ansc_biomass)] <- 0.0001
  
  #3----
  data_agg3=subset(data_agg,Seascapenr==3)
  data_agg3$ansc_biomass= nthroot(data_agg3$biomass_kg+(3/8),2)
  data_agg3$scaled_ansc_biomass=data_agg3$ansc_biomass/max(data_agg3$ansc_biomass)
  data_agg3$scaled_ansc_biomass=scaleTR(data_agg3$ansc_biomass)
  data_agg3$scaled_ansc_biomass[is.na(data_agg3$scaled_ansc_biomass)] <- 0.0001
  
  #4----
  data_agg4=subset(data_agg,Seascapenr==4)
  data_agg4$ansc_biomass= nthroot(data_agg4$biomass_kg+(3/8),2)
  data_agg4$scaled_ansc_biomass=data_agg4$ansc_biomass/max(data_agg4$ansc_biomass)
  data_agg4$scaled_ansc_biomass=scaleTR(data_agg4$ansc_biomass)
  data_agg4$scaled_ansc_biomass[is.na(data_agg4$scaled_ansc_biomass)] <- 0.0001
  
  #5----
  data_agg5=subset(data_agg,Seascapenr==5)
  data_agg5$ansc_biomass= nthroot(data_agg5$biomass_kg+(3/8),2)
  data_agg5$scaled_ansc_biomass=data_agg5$ansc_biomass/max(data_agg5$ansc_biomass)
  data_agg5$scaled_ansc_biomass=scaleTR(data_agg5$ansc_biomass)
  data_agg5$scaled_ansc_biomass[is.na(data_agg5$scaled_ansc_biomass)] <- 0.0001

  
  #6----
  data_agg6=subset(data_agg,Seascapenr==6)
  data_agg6$ansc_biomass= nthroot(data_agg6$biomass_kg+(3/8),2)
  data_agg6$scaled_ansc_biomass=data_agg6$ansc_biomass/max(data_agg6$ansc_biomass)
  data_agg6$scaled_ansc_biomass=scaleTR(data_agg6$ansc_biomass)
  data_agg6$scaled_ansc_biomass[is.na(data_agg6$scaled_ansc_biomass)] <- 0.0001
  
  #7----
  data_agg7=subset(data_agg,Seascapenr==7)
  data_agg7$ansc_biomass= nthroot(data_agg7$biomass_kg+(3/8),2)
  data_agg7$scaled_ansc_biomass=data_agg7$ansc_biomass/max(data_agg7$ansc_biomass)
  data_agg7$scaled_ansc_biomass=scaleTR(data_agg7$ansc_biomass)
  data_agg7$scaled_ansc_biomass[is.na(data_agg7$scaled_ansc_biomass)] <- 0.0001
  
  #8----
  data_agg8=subset(data_agg,Seascapenr==8)
  data_agg8$ansc_biomass= nthroot(data_agg8$biomass_kg+(3/8),2)
  data_agg8$scaled_ansc_biomass=data_agg8$ansc_biomass/max(data_agg8$ansc_biomass)
  data_agg8$scaled_ansc_biomass=scaleTR(data_agg8$ansc_biomass)
  data_agg8$scaled_ansc_biomass[is.na(data_agg8$scaled_ansc_biomass)] <- 0.0001
  
  #9----
  data_agg9=subset(data_agg,Seascapenr==9)
  data_agg9$ansc_biomass= nthroot(data_agg9$biomass_kg+(3/8),2)
  data_agg9$scaled_ansc_biomass=data_agg9$ansc_biomass/max(data_agg9$ansc_biomass)
  data_agg9$scaled_ansc_biomass=scaleTR(data_agg9$ansc_biomass)
  data_agg9$scaled_ansc_biomass[is.na(data_agg9$scaled_ansc_biomass)] <- 0.0001
  
  #10----
  data_agg10=subset(data_agg,Seascapenr==10)
  data_agg10$ansc_biomass= nthroot(data_agg10$biomass_kg+(3/8),2)
  data_agg10$scaled_ansc_biomass=data_agg10$ansc_biomass/max(data_agg10$ansc_biomass)
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
  
  
  mod=bam(scaled_ansc_biomass~Seascapenr+s(Year, by=Seascapenr),data=data_agg_ansc,method="REML",family=betar)
  print("gam fitted")
  #print(AIC(mod))
  #par(mfrow=c(2,2))
  #gam.check(mod)
  #list_AIC_beta[[j]]=AIC(mod)
  
  #dev.copy(png,file=paste(path,sprintf("Seascapes/Q1/%s/GAM/GAM_ANSC_diagnostics_3.png",species_info$file_name[j]),sep=""))
  #dev.off()
  
  
  # # ################################################################
  # # # Plot predicted smooths with all data                         #----
  # # ################################################################
  # newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9,10))
  # pred <- predict(mod, newd, se.fit = TRUE,type="response")
  # se.fit <- pred$se.fit
  # pred <- transform(cbind(data.frame(pred), newd),
  #                   uprP = fit + (2 * se.fit),
  #                   lwrP = fit - (2 * se.fit))
  # 
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
  # 
  # img=readPNG(source=paste(path,sprintf("images/spec_%s.png",species_info$file_name[j]),sep=""))
  # g <- rasterGrob(img, interpolate=TRUE)
  # #
  # # ### plot smooths with all data -----
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
  # ggsave(cplot,file=paste(path,sprintf("Seascapes/Q1/%s/GAM/GAM_ANSC_data_smooths_5.png",species_info$file_name[j]),sep=""),height=8, width=10,dpi = 600)#
  # 

  # ################################################################
  # # Plot predicted smooths at smooth scale and simultaneous CI   #----
  # ################################################################
  newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9,10))
  pred <- predict(mod, newd, se.fit = TRUE)
  se.fit <- pred$se.fit
  pred <- transform(cbind(data.frame(pred), newd),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit))
  ###
  Vb <- vcov(mod)
  newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9,10))
  pred <- predict(mod, newd, se.fit = TRUE)
  se.fit <- pred$se.fit

  set.seed(42)
  N <- 10000

  BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

  #xp matrix where basisfunctions of the model have been evaluated at 400 time point values per area
  Cg <- predict(mod, newd, type = "lpmatrix")
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
  # ### plot smooths with all data -----
  # s1plot=ggplot(data=S1,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S1")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # s1plot=s1plot+annotation_custom(g,xmin = min(data_agg1$Year), xmax = min(data_agg1$Year)+9, ymin = max(S1$fit)-0.2*(max(S1$fit)), ymax = max(S1$fit))
  # 
  # 
  # s2plot=ggplot(data=S2,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S2")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # 
  # 
  # s3plot=ggplot(data=S3,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S3")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # 
  # 
  # s4plot=ggplot(data=S4,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S4")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # 
  # s5plot=ggplot(data=S5,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S5")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # 
  # s6plot=ggplot(data=S6,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S6")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # 
  # s7plot=ggplot(data=S7,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S7")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # s8plot=ggplot(data=S8,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S8")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # 
  # s9plot=ggplot(data=S9,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S9")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # 
  # s10plot=ggplot(data=S10,aes(x=Year))+
  #   theme_bw()+
  #   geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
  #   geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
  #   geom_line(aes(y=fit),col="black",size=1)+
  #   ggtitle("S10")+
  #   labs(x="Year",y="ANSC fish biomass trend")
  # 
  # cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,s10plot,ncol=5,nrow=2)
  # ggsave(cplot,file=paste(path,sprintf("Seascapes/Q1/%s/GAM/GAM_ANSC_smooths_3.png",species_info$file_name[j]),sep=""),height=5, width=9,dpi = 600)#
  # 
  # #save residual autocorrelation plot
  # par(mfrow=c(1,2))
  # acf(residuals(mod),main="ACF(mod)")
  # pacf(residuals(mod),main="partial ACF(mod)")
  # dev.copy(png,file=paste(path,sprintf("Seascapes/Q1/%s/GAM/GAM_ANSC_autocorrelation_3.png",species_info$file_name[j]),sep=""))
  # dev.off()

  # # #####################################################
  # # # Compute difference smooths using simultaneous CI's #  ----
  # # #####################################################
  combi=c(1,2,3,4,5,6,7,8,9,10)
  combi_width=as.data.frame(combinations(n=10,r=2,v=combi))
  pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9,10))
  xp <- predict(mod, pdat, type = 'lpmatrix')

  sprintf("Seascapenr%s",combi_width$V1[2])
  sprintf("Seascapenr%s",combi_width$V2[2])
  

  diff_df=data.frame(matrix(0, ncol = 3, nrow = 45))
  names(diff_df)=c("combi","period_diff","ci_width")
  for (i in 1:nrow(combi_width)){
    print(i)
    pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9,10))
    xp <- predict(mod, pdat, type = 'lpmatrix')
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
    Cg <- predict(mod, pdat, type = "lpmatrix")
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
        
     
    # DifS12plot=
    #   ggplot(S12pred, aes(x = Year, y = V1)) +
    #   theme_bw()+
    #   #geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
    #   geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
    #   geom_line(data=S12pred,aes(x=Year,y=V1),col="black",size=0.5) +
    #   geom_line(data=S12higher,aes(x=Year,y=V1),col="seagreen1",size=1)+
    #   geom_line(data=S12lower,aes(x=Year,y=V1),col="red",size=1)+
    #   ggtitle(sprintf("%s Seascape %s - %s",species_info$Scientific.name[j],combi_width$V1[i],combi_width$V2[i]))+
    #   labs(x = "Year", y = 'Difference Fish biomass trend')+
    #   theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
    #   theme(axis.text.x = element_text(face="bold",size=9))+
    #   theme(axis.text.y = element_text(face="bold",size=9))+
    #   theme(axis.title.x = element_text(face="bold",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
    #   theme(axis.title.y = element_text(face="bold",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
    #   theme(axis.line = element_line(size=1))
  #   ggsave(DifS12plot,file = paste(path,sprintf("Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooth_%s_%s_ANSC.png",species_info$file_name[j],combi_width$V1[i],combi_width$V2[i]),sep=""), height=3.5, width=4.5,dpi = 600)
  #DifS12plot
  # 
  }
  write.csv(diff_df,file = paste(path,sprintf("Seascapes/Q1/%s/GAM/Diff_smooth/num_diff.csv",species_info$file_name[j]),sep=""))
}


#save different smooths in single large image            
for (j in 1:16){
  print(j)
  p=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth",species_info$file_name[j]),sep="")
  filenames=list.files(path=p,pattern="*ANSC.png",full.names=TRUE)
  
  rl <- lapply(filenames, png::readPNG)
  gl <- lapply(rl, grid::rasterGrob)
  cp=grid.arrange(gl[[1]],gl[[2]],gl[[3]],gl[[4]],gl[[5]],gl[[6]],
                  gl[[7]],gl[[8]],gl[[9]],gl[[10]],gl[[11]],gl[[12]],
                  gl[[13]],gl[[14]],gl[[15]],gl[[16]],gl[[17]],gl[[18]],
                  gl[[19]],gl[[20]],gl[[21]],gl[[22]],gl[[23]],gl[[24]],
                  gl[[25]],gl[[26]],gl[[27]],gl[[28]],gl[[29]],gl[[30]],
                  gl[[31]],gl[[32]],gl[[33]],gl[[34]],gl[[35]],gl[[36]],ncol=6,nrow=6)
  ggsave(cp,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooths_ANSC.png",species_info$file_name[j]),sep=""),height=21, width=27,dpi = 300)
}






###########################################
#### compute temporal correlation file #### ----
###########################################
for (j in 1:16){
  data=read.csv(paste(path,sprintf("Data/Filtered_%s.csv",species_info$file_name[j]),sep=""))
  names(data)[names(data) == "id"] <- "Seascapenr"
  print(species_info$file_name[j])
  
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  ssc_data=completeFun(data,"Seascapenr")
  ssc_data=subset(ssc_data,Quarter==1)
  data_agg=ssc_data
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species, data=ssc_data, sum)
  
  # # #set outlier weight limit (99th percentile)
  #   nf_p=quantile(data_agg$Total_wgt,0.99)
  #   if (nf_p > 0){
  #    #subset of data excluding outlier weights
  #     data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
  #   }
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000
  
  #separate subset per seascape, anscombe transform, max scale it and finally rbind back into single dataset----
  #1----
  data_agg1=subset(data_agg,Seascapenr==1)
  data_agg1$ansc_biomass= nthroot(data_agg1$biomass_kg+(3/8),2)
  data_agg1$scaled_ansc_biomass=data_agg1$ansc_biomass/max(data_agg1$ansc_biomass)
  data_agg1$scaled_ansc_biomass=scaleTR(data_agg1$ansc_biomass)
  data_agg1$scaled_ansc_biomass[is.na(data_agg1$scaled_ansc_biomass)] <- 0.0002
  
  #2-----
  data_agg2=subset(data_agg,Seascapenr==2)
  data_agg2$ansc_biomass= nthroot(data_agg2$biomass_kg+(3/8),2)
  data_agg2$scaled_ansc_biomass=data_agg2$ansc_biomass/max(data_agg2$ansc_biomass)
  data_agg2$scaled_ansc_biomass=scaleTR(data_agg2$ansc_biomass)
  data_agg2$scaled_ansc_biomass[is.na(data_agg2$scaled_ansc_biomass)] <- 0.0002
  
  #3----
  data_agg3=subset(data_agg,Seascapenr==3)
  data_agg3$ansc_biomass= nthroot(data_agg3$biomass_kg+(3/8),2)
  data_agg3$scaled_ansc_biomass=data_agg3$ansc_biomass/max(data_agg3$ansc_biomass)
  data_agg3$scaled_ansc_biomass=scaleTR(data_agg3$ansc_biomass)
  data_agg3$scaled_ansc_biomass[is.na(data_agg3$scaled_ansc_biomass)] <- 0.0002
  
  #4----
  data_agg4=subset(data_agg,Seascapenr==4)
  data_agg4$ansc_biomass= nthroot(data_agg4$biomass_kg+(3/8),2)
  data_agg4$scaled_ansc_biomass=data_agg4$ansc_biomass/max(data_agg4$ansc_biomass)
  data_agg4$scaled_ansc_biomass=scaleTR(data_agg4$ansc_biomass)
  data_agg4$scaled_ansc_biomass[is.na(data_agg4$scaled_ansc_biomass)] <- 0.0002
  
  #5----
  data_agg5=subset(data_agg,Seascapenr==5)
  data_agg5$ansc_biomass= nthroot(data_agg5$biomass_kg+(3/8),2)
  data_agg5$scaled_ansc_biomass=data_agg5$ansc_biomass/max(data_agg5$ansc_biomass)
  data_agg5$scaled_ansc_biomass=scaleTR(data_agg5$ansc_biomass)
  data_agg5$scaled_ansc_biomass[is.na(data_agg5$scaled_ansc_biomass)] <- 0.0002
  
  
  #6----
  data_agg6=subset(data_agg,Seascapenr==6)
  data_agg6$ansc_biomass= nthroot(data_agg6$biomass_kg+(3/8),2)
  data_agg6$scaled_ansc_biomass=data_agg6$ansc_biomass/max(data_agg6$ansc_biomass)
  data_agg6$scaled_ansc_biomass=scaleTR(data_agg6$ansc_biomass)
  data_agg6$scaled_ansc_biomass[is.na(data_agg6$scaled_ansc_biomass)] <- 0.0002
  
  #7----
  data_agg7=subset(data_agg,Seascapenr==7)
  data_agg7$ansc_biomass= nthroot(data_agg7$biomass_kg+(3/8),2)
  data_agg7$scaled_ansc_biomass=data_agg7$ansc_biomass/max(data_agg7$ansc_biomass)
  data_agg7$scaled_ansc_biomass=scaleTR(data_agg7$ansc_biomass)
  data_agg7$scaled_ansc_biomass[is.na(data_agg7$scaled_ansc_biomass)] <- 0.0002
  
  #8----
  data_agg8=subset(data_agg,Seascapenr==8)
  data_agg8$ansc_biomass= nthroot(data_agg8$biomass_kg+(3/8),2)
  data_agg8$scaled_ansc_biomass=data_agg8$ansc_biomass/max(data_agg8$ansc_biomass)
  data_agg8$scaled_ansc_biomass=scaleTR(data_agg8$ansc_biomass)
  data_agg8$scaled_ansc_biomass[is.na(data_agg8$scaled_ansc_biomass)] <- 0.0002
  
  #9----
  data_agg9=subset(data_agg,Seascapenr==9)
  data_agg9$ansc_biomass= nthroot(data_agg9$biomass_kg+(3/8),2)
  data_agg9$scaled_ansc_biomass=data_agg9$ansc_biomass/max(data_agg9$ansc_biomass)
  data_agg9$scaled_ansc_biomass=scaleTR(data_agg9$ansc_biomass)
  data_agg9$scaled_ansc_biomass[is.na(data_agg9$scaled_ansc_biomass)] <- 0.0002
  
  #10----
  data_agg10=subset(data_agg,Seascapenr==10)
  data_agg10$ansc_biomass= nthroot(data_agg10$biomass_kg+(3/8),2)
  data_agg10$scaled_ansc_biomass=data_agg10$ansc_biomass/max(data_agg10$ansc_biomass)
  data_agg10$scaled_ansc_biomass=scaleTR(data_agg10$ansc_biomass)
  data_agg10$scaled_ansc_biomass[is.na(data_agg10$scaled_ansc_biomass)] <- 0.0002
  
  #rbind----
  data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9,data_agg10)
  
  ####################################################
  #GAM model with separate smoother per seascape
  ###################################################
  print("gam")
  data_agg_ansc$Seascapenr=as.factor(data_agg_ansc$Seascapenr) #Set as factor
  data_agg_ansc = start_event(data_agg_ansc, column="Year",event="Seascapenr") #Sep time series per factor level
  
  
  mod=bam(scaled_ansc_biomass~Seascapenr+s(Year, by=Seascapenr),data=data_agg_ansc,method="REML",family=betar)
  
  print("gam fitted")
  
  combi=c(1,2,3,4,5,6,7,8,9,10)
  combi_width=as.data.frame(combinations(n=10,r=2,v=combi))
  combi_width2=as.matrix(combinations(n=10,r=2,v=combi))
  
  newdat = as.data.frame(1977:2019)
  colnames(newdat)="Year"
  
  #100 simulations from model posterior distribution
  sim_dat=as.data.frame(smooth_samples(mod,n=100,new_data=newdat,n_vals=43))
  sim_dat2=data.matrix(smooth_samples(mod,n=100,new_data=newdat,n_vals=43))
  Year=floor(as.numeric(sim_dat2[,".x1"]))
  sim_dat2=cbind(sim_dat2,Year)
  #plot(x=sim_dat$Year,y=sim_dat$value)
  sv1=sim_dat2[sim_dat2[,"Seascapenr"]==1,]
  draw_v1=sv1[sv1[,"draw"]==1,]
  #go iteratively through smooth pair combinations and simulate autocorrelation
  sim_autocor=matrix(ncol=4,nrow=0)
  draw_v1[,8]
  autocor=apply(combi_width,MARGIN = 1,
          FUN=function(x){
          V1=x[1]#access element first column
          V2=x[2]#access element second column
          c_var=paste(V1,"&",V2,sep=" ")
          print(c_var)
          
          sv1=sim_dat2[sim_dat2[,"Seascapenr"]==V1,] #subset seascape x
          sv2=sim_dat2[sim_dat2[,"Seascapenr"]==V2,] #subset seascape y
          
          simulations=c(seq(1:100))
          
          sapply(simulations,FUN=function(y){

            draw_v1=sv1[sv1[,"draw"]==y[1],] # draw n subset seascape x
            draw_v2=sv2[sv2[,"draw"]==y[1],] # draw n subset seascape y
            
            sim=as.numeric(y[1]) 
            
             dat=L_AXHA(draw_v1[,8],draw_v2[,8],lags=c(0,1,2,3,4,5,6,7,8,9,10),L=5) #temporal correlation
             dat=cbind(dat,c(0,1,2,3,4,5,6,7,8,9,10)) # lags
             dat=cbind(dat,c(sim,sim,sim,sim,sim,sim,sim,sim,sim,sim,sim)) # draw
             dat=cbind(dat,c(c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var)) #seascapes

             return(dat)
            })
         }
        )
    
  #prepare as csv
  df  <-  as.data.frame(autocor)
  df=data.frame(t(df))
  df=as.data.frame(lapply(df,unlist))
  colnames(df)=c("corr","lag","lag_dup","sim","Seascapes")
  df$lag_dup=NULL
  write.csv(df,paste(path,sprintf("Seascapes/Q1/%s/GAM/CCF/autocor_lag.csv",species_info$file_name[j]),sep=""))
}
#     combi_sim=function(x){
#     V1=x[,1]
#     V2=x[,2]
#     ci=c(V1,V2)
#     return(ci)
#     
#   }
#   ci
#   for (item in 1:nrow(combi_width)){ 
#     print(item)
#     V1=sprintf("s(Year):Seascapenr%s",combi_width$V1[item])
#     V2=sprintf("s(Year):Seascapenr%s",combi_width$V2[item])
#     
#     sim_dat_V1=subset(sim_dat,term==V1)
#     sim_dat_V2=subset(sim_dat,term==V2)
#     c_var=paste(combi_width$V1[item],"&",combi_width$V2[item],sep=" ")
# 
#     
#     for (sim in 1:100){
#         #print(sim)
#         draw_V1=subset(sim_dat_V1,draw==sim)
#         draw_V2=subset(sim_dat_V2,draw==sim)
#         
#         #Call temporal autocorrelation function
#         dat=lagged_DCCA_CC(x=draw_V1$value,y=draw_V2$value,lags=c(0,1,2,3,4,5,6,7,8,9,10),k=5)
#         dat=cbind(dat,c(0,1,2,3,4,5,6,7,8,9,10))
#         dat=cbind(dat,c(sim,sim,sim,sim,sim,sim,sim,sim,sim,sim,sim))
#         dat=cbind(dat,c(c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var,c_var))
#         sim_autocor=rbind(sim_autocor,dat)
#         #dat$Seascapes=Seascapes
#         #sim_autocor=rbind(sim_autocor,dat)
#       } 
#   }
#   sim_autocor=as.data.frame(sim_autocor)
#   #write autocor file to csv
#   write.csv(sim_autocor,paste(path,sprintf("Seascapes/Q1/%s/GAM/CCF/autocor_lag.csv",species_info$file_name[j]),sep=""))
# }


#######################################
### Create autocor graphs from file ###----
#######################################

#seascape combinations to recall
combi=c(1,2,3,4,5,6,7,8,9,10)
combi_width=as.data.frame(combinations(n=10,r=2,v=combi))

for (j in 1:16){
  data=read.csv(paste(path,sprintf("Seascapes/Q1/%s/GAM/CCF/autocor_lag.csv",species_info$file_name[j]),sep=""))
  names(data)[1:5]=c("ID","cor","lag","sim","Seascapes")
  print(species_info$file_name[j])
  for (i in 1:nrow(combi_width)){
    sc_combi=sprintf("%s & %s",combi_width$V1[i],combi_width$V2[i])
  
    #subset data
    combi_data=subset(data, Seascapes==sc_combi)
    combi_data=select(combi_data,-ID,sim)#drop unneccesary columns
    
    # aggregate autocorrelation per lag (mean per lag) for subset
    combi_data2=do.call(data.frame, aggregate(.~lag+Seascapes,data=combi_data,function(x)c(mean=mean(x),se=std_err(x))))
    combi_data2$upr_conf=combi_data2$cor.mean+(2*combi_data2$cor.se)
    combi_data2$lwr_conf=combi_data2$cor.mean-(2*combi_data2$cor.se)
    
    #create graph
    autocor_plot=ggplot(combi_data2, aes(x = lag, y = cor.mean))+
                 theme_classic()+
                 geom_bar(stat="identity",col="black",fill="gray",alpha=0.5,width=0.5)+#col="black",fill="gray",width=0.5)+
                 geom_errorbar(aes(ymin=lwr_conf,ymax=upr_conf),width=0.2)+

                 geom_hline(yintercept=0.3, linetype="dashed", color = "red",size=1,alpha=0.5)+
                 geom_hline(yintercept=-0.3, linetype="dashed", color = "red",size=1,alpha=0.5)+
                 geom_hline(yintercept=0, color = "black",size=0.5)+


                 ggtitle(sprintf("%s correlation S%s & S%s",species_info$Scientific.name[j],combi_width$V1[i],combi_width$V2[i]))+
                 labs(x = "Lag (Year)", y = 'AXHA correlation coefficient')+
                 scale_x_continuous(breaks = round(seq(min(combi_data2$lag), max(combi_data2$lag), by = 1),1))+
                 #scale_y_continuous(breaks = round(seq(min(combi_data$cor), max(combi_data$cor), by = 0.05),1))+

                  theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
                  theme(axis.text.x = element_text(face="bold",size=9))+
                  theme(axis.text.y = element_text(face="bold",size=9))+
                  theme(axis.title.x = element_text(face="bold",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
                  theme(axis.title.y = element_text(face="bold",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
                  ylim(-1,1)
    autocor_plot
    ggsave(autocor_plot,file = paste(path,sprintf("Seascapes/Q1/%s/GAM/CCF/CCF_smooth_%s_%s_ANSC.png",species_info$file_name[j],combi_width$V1[i],combi_width$V2[i]),sep=""), height=3.5, width=4.5,dpi = 600)
    
  } 
}

####################################################################
###figure to illustrate fitting of gam and posterior simulations####
####################################################################
data=read.csv(paste(path,sprintf("Data/Filtered_Limanda_limanda.csv",species_info$file_name[1]),sep=""))

#dataprep----
names(data)[names(data) == "id"] <- "Seascapenr"
print(species_info$file_name[1])

# Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
data$ID = cumsum(!duplicated(data[3:7]))

ssc_data=completeFun(data,"Seascapenr")
ssc_data=subset(ssc_data,Quarter==1)
data_agg=ssc_data
data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species, data=ssc_data, sum)

# # #set outlier weight limit (99th percentile)
#   nf_p=quantile(data_agg$Total_wgt,0.99)
#   if (nf_p > 0){
#    #subset of data excluding outlier weights
#     data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
#   }

# convert weights from grams to kg
data_agg$biomass_kg=data_agg$Total_wgt/1000

#separate subset per seascape, anscombe transform, max scale it and finally rbind back into single dataset----
  #1----
  data_agg1=subset(data_agg,Seascapenr==1)
  data_agg1$ansc_biomass= nthroot(data_agg1$biomass_kg+(3/8),2)
  data_agg1$scaled_ansc_biomass=data_agg1$ansc_biomass/max(data_agg1$ansc_biomass)
  data_agg1$scaled_ansc_biomass=scaleTR(data_agg1$ansc_biomass)
  data_agg1$scaled_ansc_biomass[is.na(data_agg1$scaled_ansc_biomass)] <- 0.0002
  
  #2-----
  data_agg2=subset(data_agg,Seascapenr==2)
  data_agg2$ansc_biomass= nthroot(data_agg2$biomass_kg+(3/8),2)
  data_agg2$scaled_ansc_biomass=data_agg2$ansc_biomass/max(data_agg2$ansc_biomass)
  data_agg2$scaled_ansc_biomass=scaleTR(data_agg2$ansc_biomass)
  data_agg2$scaled_ansc_biomass[is.na(data_agg2$scaled_ansc_biomass)] <- 0.0002
  
  #3----
  data_agg3=subset(data_agg,Seascapenr==3)
  data_agg3$ansc_biomass= nthroot(data_agg3$biomass_kg+(3/8),2)
  data_agg3$scaled_ansc_biomass=data_agg3$ansc_biomass/max(data_agg3$ansc_biomass)
  data_agg3$scaled_ansc_biomass=scaleTR(data_agg3$ansc_biomass)
  data_agg3$scaled_ansc_biomass[is.na(data_agg3$scaled_ansc_biomass)] <- 0.0002
  
  #4----
  data_agg4=subset(data_agg,Seascapenr==4)
  data_agg4$ansc_biomass= nthroot(data_agg4$biomass_kg+(3/8),2)
  data_agg4$scaled_ansc_biomass=data_agg4$ansc_biomass/max(data_agg4$ansc_biomass)
  data_agg4$scaled_ansc_biomass=scaleTR(data_agg4$ansc_biomass)
  data_agg4$scaled_ansc_biomass[is.na(data_agg4$scaled_ansc_biomass)] <- 0.0002
  
  #5----
  data_agg5=subset(data_agg,Seascapenr==5)
  data_agg5$ansc_biomass= nthroot(data_agg5$biomass_kg+(3/8),2)
  data_agg5$scaled_ansc_biomass=data_agg5$ansc_biomass/max(data_agg5$ansc_biomass)
  data_agg5$scaled_ansc_biomass=scaleTR(data_agg5$ansc_biomass)
  data_agg5$scaled_ansc_biomass[is.na(data_agg5$scaled_ansc_biomass)] <- 0.0002
  
  
  #6----
  data_agg6=subset(data_agg,Seascapenr==6)
  data_agg6$ansc_biomass= nthroot(data_agg6$biomass_kg+(3/8),2)
  data_agg6$scaled_ansc_biomass=data_agg6$ansc_biomass/max(data_agg6$ansc_biomass)
  data_agg6$scaled_ansc_biomass=scaleTR(data_agg6$ansc_biomass)
  data_agg6$scaled_ansc_biomass[is.na(data_agg6$scaled_ansc_biomass)] <- 0.0002
  
  #7----
  data_agg7=subset(data_agg,Seascapenr==7)
  data_agg7$ansc_biomass= nthroot(data_agg7$biomass_kg+(3/8),2)
  data_agg7$scaled_ansc_biomass=data_agg7$ansc_biomass/max(data_agg7$ansc_biomass)
  data_agg7$scaled_ansc_biomass=scaleTR(data_agg7$ansc_biomass)
  data_agg7$scaled_ansc_biomass[is.na(data_agg7$scaled_ansc_biomass)] <- 0.0002
  
  #8----
  data_agg8=subset(data_agg,Seascapenr==8)
  data_agg8$ansc_biomass= nthroot(data_agg8$biomass_kg+(3/8),2)
  data_agg8$scaled_ansc_biomass=data_agg8$ansc_biomass/max(data_agg8$ansc_biomass)
  data_agg8$scaled_ansc_biomass=scaleTR(data_agg8$ansc_biomass)
  data_agg8$scaled_ansc_biomass[is.na(data_agg8$scaled_ansc_biomass)] <- 0.0002
  
  #9----
  data_agg9=subset(data_agg,Seascapenr==9)
  data_agg9$ansc_biomass= nthroot(data_agg9$biomass_kg+(3/8),2)
  data_agg9$scaled_ansc_biomass=data_agg9$ansc_biomass/max(data_agg9$ansc_biomass)
  data_agg9$scaled_ansc_biomass=scaleTR(data_agg9$ansc_biomass)
  data_agg9$scaled_ansc_biomass[is.na(data_agg9$scaled_ansc_biomass)] <- 0.0002
  
  #10----
  data_agg10=subset(data_agg,Seascapenr==10)
  data_agg10$ansc_biomass= nthroot(data_agg10$biomass_kg+(3/8),2)
  data_agg10$scaled_ansc_biomass=data_agg10$ansc_biomass/max(data_agg10$ansc_biomass)
  data_agg10$scaled_ansc_biomass=scaleTR(data_agg10$ansc_biomass)
  data_agg10$scaled_ansc_biomass[is.na(data_agg10$scaled_ansc_biomass)] <- 0.0002
  
  #rbind----
  data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9,data_agg10)
  
####################################################
#GAM model with separate smoother per seascape ----
###################################################
print("gam")
data_agg_ansc$Seascapenr=as.factor(data_agg_ansc$Seascapenr) #Set as factor
data_agg_ansc = start_event(data_agg_ansc, column="Year",event="Seascapenr") #Sep time series per factor level


mod=bam(scaled_ansc_biomass~Seascapenr+s(Year, by=Seascapenr),data=data_agg_ansc,method="REML",family=betar)

sim_dat=as.data.frame(smooth_samples(mod,n=100,new_data=newdat,n_vals=43))

newd <- expand.grid(Year = seq(1977, 2019, length = 43), Seascapenr = c(1,2,3,4,5,6,7,8,9,10))
pred <- predict(mod, newd, se.fit = TRUE)
se.fit <- pred$se.fit
pred <- transform(cbind(data.frame(pred), newd),
                  uprP = fit + (2 * se.fit),
                  lwrP = fit - (2 * se.fit))
###
Vb <- vcov(mod)
newd <- expand.grid(Year = seq(1977, 2019, length = 43), Seascapenr = c(1,2,3,4,5,6,7,8,9,10))
pred <- predict(mod, newd, se.fit = TRUE)
se.fit <- pred$se.fit

set.seed(42)
N <- 10000

BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

#xp matrix where basisfunctions of the model have been evaluated at 400 time point values per area
Cg <- predict(mod, newd, type = "lpmatrix")
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
sim_dat_S1=subset(sim_dat,Seascapenr==1)
sim_dat_S1$uprS=rep(S1$uprS,100)#+2.8,100)
sim_dat_S1$lwrS=rep(S1$lwrS,100)#+2.8,100)
sim_dat_S2=subset(sim_dat,Seascapenr==2)

#GGPLOTS---
#sim1
s1_plot=ggplot(sim_dat_S1, aes(x = .x1, y = value, group=draw))+
  theme_classic()+
  #geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.005,fill="blue")+
  geom_line(alpha=0.1)+
  ggtitle(sprintf("Posterior simulations seascape %s",sim_dat_S1$Seascapenr[1]))+
  labs(x = "Year", y = 'Fish biomass trend')+
  #scale_x_continuous(breaks = round(seq(min(combi_data2$lag), max(combi_data2$lag), by = 1),1))+
  theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
  theme(axis.text.x = element_text(face="bold",size=12))+
  theme(axis.text.y = element_text(face="bold",size=12))+
  theme(axis.title.x = element_text(face="bold",size=12,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  theme(axis.title.y = element_text(face="bold",size=12,margin = margin(t = 0, r = 10, b = 0, l = 0)))
s1_plot
ggsave(s1_plot,file = paste(path,sprintf("images/simulated_trends1.png"),sep=""), height=3.5, width=4.5,dpi = 600)

#sim2
s2_plot=ggplot(sim_dat_S2, aes(x = .x1, y = value, group=draw))+
  theme_classic()+
  geom_line(alpha=0.1)+
  ggtitle(sprintf("Posterior simulations seascape %s",sim_dat_S2$Seascapenr[1]))+
  labs(x = "Year", y = 'Fish biomass trend')+
  #scale_x_continuous(breaks = round(seq(min(combi_data2$lag), max(combi_data2$lag), by = 1),1))+
  theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
  theme(axis.text.x = element_text(face="bold",size=12))+
  theme(axis.text.y = element_text(face="bold",size=12))+
  theme(axis.title.x = element_text(face="bold",size=12,margin = margin(t = 10, r = 0, b = 0, l = 5)))+
  theme(axis.title.y = element_text(face="bold",size=12,margin = margin(t = 0, r = 10, b = 0, l = 0)))
s2_plot
ggsave(s2_plot,file = paste(path,sprintf("images/simulated_trends2.png"),sep=""), height=3.5, width=4.5,dpi = 600)

?smooth_samples

getAnywhere(smooth_samples())
smooth_samples
showMethods(smooth_samples)
getMethod(smooth_samples)
