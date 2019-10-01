#load packages----
library(plyr)
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
species_info=read.table(paste(path,"spec_info.txt",sep=""),sep="\t",header=TRUE)
species_info$file_name =gsub(" ","_",species_info$Scientific.name)

#open seascapes map and convert to wgs84 long,lat----
ssc_map=readOGR("C:/Users/mrademaker/Documents/Research projects/STCNWS/Seascapes/.",layer="seascapes")
ssc_map_wgs84 <- spTransform(ssc_map, CRS("+proj=longlat +datum=WGS84")) 

#filter and adjust data_sets----
for (i in 14:16){
  print(i)
  file_name=paste(path,species_info$file_name[i],sep="")
  
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
  
  #save filtered file
  file_name=sprintf("C:/Users/mrademaker/Documents/Research projects/STCNWS/DATRAS/cpue_length_hour/Filtered_%s.csv",species_info$file_name[i])
  write.csv(data_ssc,file_name)
  #
  
  
  
}



####For loop - Plot trend in biomass (kg) ---- 
for (i in 14:16){
  print(paste(i,species_info$file_name[i],sep=","))
  # Read in data
  data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[i]),sep=""))
  
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  #Aggregate weight over each haul (sum)
  data_agg=aggregate(Total_wgt ~ ID+Year, data=data, sum)
  
  #set outlier weight limit (95th percentile)
  nf_p=quantile(data_agg$Total_wgt,0.95) 
  
  #subset of data excluding outlier weights
  data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000 
  
  #mean per year to get better overview
  data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
  
  #Plot
  p=ggplot(data_agg2,aes(Year,biomass_kg))+
    geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
    geom_point(size=3, alpha=0.5)+
    theme_classic()+
    ggtitle(species_info$Scientific.name[i])+
    labs(x = "Year", y = 'Mean CPUE (kg)')+
    theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
    theme(axis.text.x = element_text(face="bold",size=18))+
    theme(axis.text.y = element_text(face="bold",size=18))+
    theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
    theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    theme(axis.line = element_line(size=1))
  p
  
  #Add species image
  img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p2=p+annotation_custom(g,xmin = 1977, xmax = 1986, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
  ggsave(p2, file = paste(path,sprintf("graphs/Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
}


####Separate plots for all data (mix Q1,Q3), clean Q1: 1977-2018 & clean Q1-Q3: 1991-2018 + set up indexes
for (i in 1:12){
  print(paste(i,species_info$file_name[i],sep=","))
  # Read in data
  data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[i]),sep=""))
  
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
 
  #All data----
  #Aggregate weight over each haul (sum)
  data_agg=aggregate(Total_wgt ~ ID+Year, data=data, sum)
  
  #set outlier weight limit (95th percentile)
  nf_p=quantile(data_agg$Total_wgt,0.95) 
  
  #subset of data excluding outlier weights
  data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000 
  
  #mean per year to get better overview
  data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
  
  
  #Plot actual biomass
  p=ggplot(data_agg2,aes(Year,biomass_kg))+
    geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
    geom_point(size=3, alpha=0.5)+
    theme_classic()+
    ggtitle(species_info$Scientific.name[i])+
    labs(x = "Year", y = 'Mean CPUE (kg)')+
    theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
    theme(axis.text.x = element_text(face="bold",size=18))+
    theme(axis.text.y = element_text(face="bold",size=18))+
    theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
    theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    theme(axis.line = element_line(size=1))
  p
  
  #Add species image
  img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$biomass_kg)-0.15*max(data_agg2$biomass_kg),ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
  ggsave(p2, file = paste(path,sprintf("graphs/All/Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
  
  #index
  #add 1 to each biomass value to exclude potential problems with zero starting values
  data_agg2$biomass_kg=data_agg2$biomass_kg+1
  base = data_agg2$biomass_kg[1]
  print(data_agg2$biomass_kg[1])
  print(base)
  
  data_agg2$indexed_biomass = (data_agg2$biomass_kg/base)*100
  
  #Plot indexed biomass
  p=ggplot(data_agg2,aes(Year,indexed_biomass))+
    geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
    geom_point(size=3, alpha=0.5)+
    theme_classic()+
    ggtitle(species_info$Scientific.name[i])+
    labs(x = "Year", y = 'Indexed biomass')+
    theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
    theme(axis.text.x = element_text(face="bold",size=18))+
    theme(axis.text.y = element_text(face="bold",size=18))+
    theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
    theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    theme(axis.line = element_line(size=1))
  p
  
  #Add species image
  img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  #different image size for callionymus lyra or will not display proportionally
  if(species_info$file_name[i]!= "Callionymus_lyra"){
  p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.15*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
  ggsave(p2, file = paste(path,sprintf("graphs/All/Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
  } else{
    p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.25*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
    ggsave(p2, file = paste(path,sprintf("graphs/All/Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
    
  }
  
  #Q1 (1977-2018) ---- 
  Q1_data = subset(data,Quarter==1)
  write.csv(Q1_data,paste(path,sprintf("Q1_filtered_%s.csv",species_info$file_name[i]),sep=""))
  #Aggregate weight over each haul (sum)
  data_agg=aggregate(Total_wgt ~ ID+Year, data=Q1_data, sum)
  
  #set outlier weight limit (95th percentile)
  nf_p=quantile(data_agg$Total_wgt,0.95) 
  
  #subset of data excluding outlier weights
  data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000 
  
  #mean per year to get better overview
  data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
  
  #Plot
  p=ggplot(data_agg2,aes(Year,biomass_kg))+
    geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
    geom_point(size=3, alpha=0.5)+
    theme_classic()+
    ggtitle(species_info$Scientific.name[i])+
    labs(x = "Year", y = 'Mean CPUE (kg)')+
    theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
    theme(axis.text.x = element_text(face="bold",size=18))+
    theme(axis.text.y = element_text(face="bold",size=18))+
    theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
    theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    theme(axis.line = element_line(size=1))
  p
  
  #Add species image
  img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p2=p+annotation_custom(g,xmin = 1977, xmax = 1986, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
  ggsave(p2, file = paste(path,sprintf("graphs/Q1/Q1_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
  
  #index
  #add 1 to each biomass value to exclude potential problems with zero starting values
  data_agg2$biomass_kg=data_agg2$biomass_kg+1
  base = data_agg2$biomass_kg[1]
  print(data_agg2$biomass_kg[1])
  print(base)
  
  data_agg2$indexed_biomass = (data_agg2$biomass_kg/base)*100
  
  #Plot indexed biomass
  p=ggplot(data_agg2,aes(Year,indexed_biomass))+
    geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
    geom_point(size=3, alpha=0.5)+
    theme_classic()+
    ggtitle(species_info$Scientific.name[i])+
    labs(x = "Year", y = 'Indexed biomass')+
    theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
    theme(axis.text.x = element_text(face="bold",size=18))+
    theme(axis.text.y = element_text(face="bold",size=18))+
    theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
    theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    theme(axis.line = element_line(size=1))
  p
  
  #Add species image
  img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  #different image size for callionymus lyra or will not display proportionally
  if(species_info$file_name[i]!= "Callionymus_lyra"){
    p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.15*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
    ggsave(p2, file = paste(path,sprintf("graphs/Q1/Q1_Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
  } else{
    p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.25*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
    ggsave(p2, file = paste(path,sprintf("graphs/Q1/Q1_Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
    
  }
  
  
  #Q1-Q3 (1991-2018)----
  Q1_Q3_data = subset(data,Year>=1991)
  write.csv(Q1_Q3_data,paste(path,sprintf("Q1_Q3_filtered_%s.csv",species_info$file_name[i]),sep=""))
  Q1_Q3_data$QY = cumsum(!duplicated(Q1_Q3_data[3:4]))
  
  #Aggregate weight over each haul (sum)
  data_agg=aggregate(Total_wgt ~ ID+Year+QY, data=Q1_Q3_data, sum)
  
  #set outlier weight limit (95th percentile)
  nf_p=quantile(data_agg$Total_wgt,0.95) 
  
  #subset of data excluding outlier weights
  data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000 
  
  #mean per year to get better overview
  data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
  
  #Plot
  p=ggplot(data_agg2,aes(Year,biomass_kg))+
    geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
    geom_point(size=3, alpha=0.5)+
    theme_classic()+
    ggtitle(species_info$Scientific.name[i])+
    labs(x = "Year", y = 'Mean CPUE (kg)')+
    theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
    theme(axis.text.x = element_text(face="bold",size=18))+
    theme(axis.text.y = element_text(face="bold",size=18))+
    theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
    theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    theme(axis.line = element_line(size=1))
  p
  
  #Add species image
  img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p2=p+annotation_custom(g,xmin = 1977, xmax = 1986, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
  ggsave(p2, file = paste(path,sprintf("graphs/Q1_Q3/Q1_Q3_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )

  #index
  #add 1 to each biomass value to exclude potential problems with zero starting values
  data_agg2$biomass_kg=data_agg2$biomass_kg+1
  base = data_agg2$biomass_kg[1]
  print(data_agg2$biomass_kg[1])
  print(base)
  
  data_agg2$indexed_biomass = (data_agg2$biomass_kg/base)*100
  
  #Plot indexed biomass
  p=ggplot(data_agg2,aes(Year,indexed_biomass))+
    geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
    geom_point(size=3, alpha=0.5)+
    theme_classic()+
    ggtitle(species_info$Scientific.name[i])+
    labs(x = "Year", y = 'Indexed biomass')+
    theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
    theme(axis.text.x = element_text(face="bold",size=18))+
    theme(axis.text.y = element_text(face="bold",size=18))+
    theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
    theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    theme(axis.line = element_line(size=1))
  p
  
  #Add species image
  img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  #different image size for callionymus lyra or will not display proportionally
  if(species_info$file_name[i]!= "Callionymus_lyra"){
    p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.15*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
    ggsave(p2, file = paste(path,sprintf("graphs/Q1_Q3/Q1_Q3_Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
  } else{
    p2=p+annotation_custom(g,xmin = 1977, xmax = 1986,ymin= max(data_agg2$indexed_biomass)-0.25*(max(data_agg2$indexed_biomass)-min(data_agg2$indexed_biomass)),ymax = max(data_agg2$indexed_biomass))
    ggsave(p2, file = paste(path,sprintf("graphs/Q1_Q3/Q1_Q3_Indexed_Trend_%s.png",species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
  }
  
}




####Separate plots for all 9 seascapes all data ----
for (i in 16:16){
  print(paste(i,species_info$file_name[i],sep=","))
  # Read in data
  data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[i]),sep=""))

  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  for(j in 1:9){
      print(paste("Seascape",j,sep=" "))
      ssc_data=subset(data,Seascapenr==j)
    
      data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr, data=ssc_data, sum)
      
      #set outlier weight limit (95th percentile)
      nf_p=quantile(data_agg$Total_wgt,0.95) 
      
      #subset of data excluding outlier weights
      #data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
      
      if (nf_p > 0){
        #subset of data excluding outlier weights
        data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
      }
      
      
      # convert weights from grams to kg
      data_agg$biomass_kg=data_agg$Total_wgt/1000 
      
      #mean per year to get better overview
      data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
      
      #Plot
      p=ggplot(data_agg2,aes(Year,biomass_kg))+
        geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
        geom_point(size=3, alpha=0.5)+
        theme_classic()+
        ggtitle(sprintf("%s Seascape %s",species_info$Scientific.name[i],j))+
        labs(x = "Year", y = 'Mean CPUE (kg)')+
        theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
        theme(axis.text.x = element_text(face="bold",size=18))+
        theme(axis.text.y = element_text(face="bold",size=18))+
        theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
        theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
        theme(axis.line = element_line(size=1))
      
      #Add species image
      img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
      g <- rasterGrob(img, interpolate=TRUE)
      
      p2=p+annotation_custom(g,xmin = min(data_agg2$Year), xmax = min(data_agg2$Year)+9, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
      ggsave(p2, file = paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],j,species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
    
    }
  
    #stack in single image
    plot1 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],1,species_info$file_name[i]),sep=""))
    plot2 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],2,species_info$file_name[i]),sep=""))
    plot3 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],3,species_info$file_name[i]),sep=""))
    plot4 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],4,species_info$file_name[i]),sep=""))
    plot5 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],5,species_info$file_name[i]),sep=""))
    plot6 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],6,species_info$file_name[i]),sep=""))
    plot7 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],7,species_info$file_name[i]),sep=""))
    plot8 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],8,species_info$file_name[i]),sep=""))
    plot9 = readPNG(paste(path,sprintf("graphs/Seascapes/All/%s/%s_Trend_%s.png",species_info$file_name[i],9,species_info$file_name[i]),sep=""))
    
    cplot=grid.arrange(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot5),
                 rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),rasterGrob(plot9),ncol=3)
    ggsave(cplot, file = paste(path,sprintf("graphs/Seascapes/All/%s/compiled_Trend_%s.png",species_info$file_name[i],species_info$file_name[i]),sep=""), height=21, width=27,dpi = 600) # adjust dpi )

}




#############GAM INDEX VERSION (single averaged datapoint/year)############## -----
#############################################################################
for (j in 16:16){
  data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[j]),sep=""))
  print(species_info$file_name[j])
  
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  ssc_data=completeFun(data,"Seascapenr")
  ssc_data=subset(ssc_data,Quarter==1)
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species, data=ssc_data, sum)
  
  #set outlier weight limit (95th percentile)
  nf_p=quantile(data_agg$Total_wgt,0.95) 
  if (nf_p > 0){
    #subset of data excluding outlier weights
    data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
  }
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000 
  
  #mean per year and seascape to get better overview
  data_agg_year=aggregate(biomass_kg~Year+Seascapenr+Species,data=data_agg,mean)
  
  #separate subset per seascape, index it and finally rbind back into single dataset----
  #1----
  data_agg1=subset(data_agg_year,Seascapenr==1)
  data_agg1$biomass_kg=data_agg1$biomass_kg+100 #scale up for index
  base = data_agg1$biomass_kg[1]
  print(data_agg1$biomass_kg[1])
  print(base)
  data_agg1$indexed_biomass = (data_agg1$biomass_kg/base)*100
  data_agg1$biomass_kg=data_agg1$biomass_kg-100 #return back for actual value
  
  #2-----
  data_agg2=subset(data_agg_year,Seascapenr==2)
  data_agg2$biomass_kg=data_agg2$biomass_kg+100
  base = data_agg2$biomass_kg[1]
  print(data_agg2$biomass_kg[1])
  print(base)
  data_agg2$indexed_biomass = (data_agg2$biomass_kg/base)*100
  data_agg2$biomass_kg=data_agg2$biomass_kg-100 #return back for actual value
  
  #3----
  data_agg3=subset(data_agg_year,Seascapenr==3)
  data_agg3$biomass_kg=data_agg3$biomass_kg+100
  base = data_agg3$biomass_kg[1]
  print(data_agg3$biomass_kg[1])
  print(base)
  data_agg3$indexed_biomass = (data_agg3$biomass_kg/base)*100
  data_agg3$biomass_kg=data_agg3$biomass_kg-100 #return back for actual value
  
  #4----
  data_agg4=subset(data_agg_year,Seascapenr==4)
  data_agg4$biomass_kg=data_agg4$biomass_kg+100
  base = data_agg4$biomass_kg[1]
  print(data_agg4$biomass_kg[1])
  print(base)
  data_agg4$indexed_biomass = (data_agg4$biomass_kg/base)*100
  data_agg4$biomass_kg=data_agg4$biomass_kg-100 #return back for actual value
  
  #5----
  data_agg5=subset(data_agg_year,Seascapenr==5)
  data_agg5$biomass_kg=data_agg5$biomass_kg+100
  base = data_agg5$biomass_kg[1]
  print(data_agg5$biomass_kg[1])
  print(base)
  data_agg5$indexed_biomass = (data_agg5$biomass_kg/base)*100
  data_agg5$biomass_kg=data_agg5$biomass_kg-100 #return back for actual value
  
  #6----
  data_agg6=subset(data_agg_year,Seascapenr==6)
  data_agg6$biomass_kg=data_agg6$biomass_kg+100
  base = data_agg6$biomass_kg[1]
  print(data_agg6$biomass_kg[1])
  print(base)
  data_agg6$indexed_biomass = (data_agg6$biomass_kg/base)*100
  data_agg6$biomass_kg=data_agg6$biomass_kg-100 #return back for actual value
  
  #7----
  data_agg7=subset(data_agg_year,Seascapenr==7)
  data_agg7$biomass_kg=data_agg7$biomass_kg+100
  base = data_agg7$biomass_kg[1]
  print(data_agg7$biomass_kg[1])
  print(base)
  data_agg7$indexed_biomass = (data_agg7$biomass_kg/base)*100
  data_agg7$biomass_kg=data_agg7$biomass_kg-100 #return back for actual value
  
  #8----
  data_agg8=subset(data_agg_year,Seascapenr==8)
  data_agg8$biomass_kg=data_agg8$biomass_kg+100
  base = data_agg8$biomass_kg[1]
  print(data_agg8$biomass_kg[1])
  data_agg8$indexed_biomass = (data_agg8$biomass_kg/base)*100
  data_agg8$biomass_kg=data_agg8$biomass_kg-100 #return back for actual value
  
  #9----
  data_agg9=subset(data_agg_year,Seascapenr==9)
  data_agg9$biomass_kg=data_agg9$biomass_kg+100
  base = data_agg9$biomass_kg[1]
  print(data_agg9$biomass_kg[1])
  print(base)
  data_agg9$indexed_biomass = (data_agg9$biomass_kg/base)*100
  data_agg9$biomass_kg=data_agg9$biomass_kg-100 #return back for actual value
  
  #rbind----
  data_agg_indexed=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9)
  
  
  ####################################################
  #GAM model with separate smoother per seascape
  ###################################################
  
  data_agg_indexed$Seascapenr=as.factor(data_agg_indexed$Seascapenr) #Set as factor
  data_agg_indexed = start_event(data_agg_indexed, column="Year",event="Seascapenr") #Sep time series per factor level
  
  #m=bam(biomass_kg~Seascapenr+s(Year, by=Seascapenr),data=data_agg_indexed,method="REML")
  #(valRho <- acf(resid(m), plot=FALSE)$acf[2])
  #AIC(m)
  
  mi=bam(indexed_biomass~Seascapenr+s(Year, by=Seascapenr),data=data_agg_indexed,method="REML")
  (valRho <- acf(resid(mi), plot=FALSE)$acf[2])
  AIC(mi)
  
  
  mi1=bam(indexed_biomass~Seascapenr+s(Year,by=Seascapenr),data=data_agg_indexed,method="REML",AR.start=data_agg_indexed$start.event,rho=valRho)
  AIC(mi1)
  
  par(mfrow=c(2,2))
  gam.check(mi1)
  #savePlot(acf,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ACF.png",species_info$file_name[j]),sep=""))
  dev.copy(png,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_diagnostics.png",species_info$file_name[j]),sep=""))
  dev.off()

  
  
  ##############################################################  #commented out/obsolete
  #Manual Trial difference smooths Seascape 1 - Seascape 2----
    # pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
    # #Control - Nurse plants
    # xp <- predict(mi1, pdat, type = 'lpmatrix')
    # col1 <- grepl("Seascapenr4", colnames(xp))
    # col2 <- grepl("Seascapenr7", colnames(xp))
    # row1 <- pdat[["Seascapenr"]] == 4
    # row2 <- pdat[["Seascapenr"]] == 7
    # 
    # 
    # ## difference rows of pred for data from comparison
    # X <- xp[row1, ] - xp[row2, ]
    # ## zero out cols of X related to splines for other lochs
    # X[, ! (col1 | col2)] <- 0
    # ## zero out the parametric cols
    # X[, !grepl('^s\\(', colnames(xp))] <- 0
    # 
    # #predicted values for differences
    # Vb=vcov(mi1)
    # dif <- X %*% coef(mi1)
    # se <- sqrt(rowSums((X %*% Vb) * X))
    # pred=cbind(dif,se)
    # 
    # ## CI of the estimates
    # set.seed(42)
    # N <- 10000
    # 
    # #predi <- predict(m1, pdat, se.fit = TRUE)
    # #se.fit <- predi$se.fit
    # 
    # BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)    
    # Cg <- predict(mi1, pdat, type = "lpmatrix")
    # simDev <- Cg %*% t(BUdiff)
    # 
    # absDev <- abs(sweep(simDev, 1, se, FUN = "/"))
    # masd <- apply(absDev, 2L, max)
    # crit <- quantile(masd, prob = 0.95, type = 8)  
    # 
    # newd <- expand.grid(Year = seq(1977, 2019, length = 126))
    # 
    # S12pred <- transform(cbind(data.frame(pred), newd),
    #                      uprP = V1 + (2 * se),
    #                      lwrP = V1 - (2 * se),
    #                      uprS = V1 + (crit * se),
    #                      lwrS = V1 - (crit * se))  
    # S12pred$sigdif <-ifelse(S12pred$lwrS<0 & S12pred$uprS<0,"lower",ifelse(S12pred$lwrS>0 & S12pred$uprS>0,"higher","no"))
    # S12higher=S12pred
    # S12higher["V1"][S12pred["sigdif"]!="higher"]=NA
    # S12lower=S12pred
    # S12lower["V1"][S12pred["sigdif"]!="lower"]=NA
    # 
    # #S12higher[S12higher$sigdif!="higher",]=NA
    # 
    # #S12higher=S12pred[which(S12pred$sigdif=="higher"),]
    # #S12pred$higher=S12pred$sigdif
    # #S12pred$higher[S12pred$higher!="higher"]=NA
    # #S12pred$lower=S12pred$sigdif
    # #S12pred$lower[S12pred$lower!="lower"]=NA
    # 
    # #mean(S12higher$V1)
    # #sd(S12higher$V1)
    # #S12lower=as.data.frame(S12pred$lower)
    # #S12higher=as.data.frame(S12pred$higher)
    # #mean(HCNlower$V1)
    # #sd(HCNlower$V1)
    # S12pred$sigdif=factor(levels="lower","higher","no")
    # 
    # DifS12plot=
    #   ggplot(S12pred, aes(x = Year, y = V1)) +
    #   theme_classic()+
    #   #geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
    #   geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.1) +
    #   geom_line(data=S12pred,aes(x=Year,y=V1),col="black",size=1) +
    #   geom_line(data=S12higher,aes(x=Year,y=V1),col="springgreen1",size=2)+
    #   geom_line(data=S12lower,aes(x=Year,y=V1),col="red",size=2)+
    #   ggtitle("Seascape 1 - 2")+
    #   labs(x = "Year", y = 'Difference Fish biomass index')+
    #   theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
    #   theme(axis.text.x = element_text(face="bold",size=18))+
    #   theme(axis.text.y = element_text(face="bold",size=18))+
    #   theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
    #   theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    #   theme(axis.line = element_line(size=1))
    # 
    # DifS12plot
  
  
  ##############################################################
  
   Vb <- vcov(mi1)
   newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
   pred <- predict(mi1, newd, se.fit = TRUE)
   se.fit <- pred$se.fit
  # 
   set.seed(42)
   N <- 10000
  # 
   BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  # 
  # xp matrix where basisfunctions of the model have been evaluated at 126 time point values per area
   Cg <- predict(mi1, newd, type = "lpmatrix")
   simDev <- Cg %*% t(BUdiff)
  # # 
   absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
   masd <- apply(absDev, 2L, max)
   crit <- quantile(masd, prob = 0.95, type = 8)
   pred <- transform(cbind(data.frame(pred), newd),
                     uprP = fit + (2 * se.fit),
                     lwrP = fit - (2 * se.fit),
                     uprS = fit + (crit * se.fit),
                     lwrS = fit - (crit * se.fit))
  # # 
  # # 
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
  
  # # # 
   s1plot=ggplot(data=S1,aes(x=Year))+
     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
     geom_line(aes(y=fit),col="black",size=1)+
     ggtitle("S1")+
     labs(x="Year",y="Indexed fish biomass trend")
   s2plot=ggplot(data=S2,aes(x=Year))+
     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
     geom_line(aes(y=fit),col="black",size=1)+
     ggtitle("S2")+
     labs(x="Year",y="Indexed fish biomass trend")
   s3plot=ggplot(data=S3,aes(x=Year))+
     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
     geom_line(aes(y=fit),col="black",size=1)+
     ggtitle("S3")+
     labs(x="Year",y="Indexed fish biomass trend")
   s4plot=ggplot(data=S4,aes(x=Year))+
     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
     geom_line(aes(y=fit),col="black",size=1)+
     ggtitle("S4")+
     labs(x="Year",y="Indexed fish biomass trend")
   s5plot=ggplot(data=S5,aes(x=Year))+
     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
     geom_line(aes(y=fit),col="black",size=1)+
     ggtitle("S5")+
     labs(x="Year",y="Indexed fish biomass trend")
   s6plot=ggplot(data=S6,aes(x=Year))+
     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
     geom_line(aes(y=fit),col="black",size=1)+
     ggtitle("S6")+
     labs(x="Year",y="Indexed fish biomass trend")
   s7plot=ggplot(data=S7,aes(x=Year))+
     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
     geom_line(aes(y=fit),col="black",size=1)+
     ggtitle("S7")+
     labs(x="Year",y="Indexed fish biomass trend")
   s8plot=ggplot(data=S8,aes(x=Year))+
     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
     geom_line(aes(y=fit),col="black",size=1)+
     ggtitle("S8")+
     labs(x="Year",y="Indexed fish biomass trend")
   s9plot=ggplot(data=S9,aes(x=Year))+
     geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
     geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
     geom_line(aes(y=fit),col="black",size=1)+
     ggtitle("S9")+
     labs(x="Year",y="Indexed fish biomass trend")
 
  cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,ncol=3,nrow=3)
  ggsave(cplot,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_smooths.png",species_info$file_name[j]),sep=""),height=7, width=9,dpi = 600)# 
  
  #save residual autocorrelation plot
  par(mfrow=c(1,1))
  acf_resid(mi1,split_pred="Seascapenr",main="ACF resid(m1)")
  dev.copy(png,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_autocorrelation.png",species_info$file_name[j]),sep=""))
  dev.off()
  
  ######################################
  #Compute difference smooths ----
  ######################################
  combi=c(1,2,3,4,5,6,7,8,9)
  combi_width=as.data.frame(combinations(n=9,r=2,v=combi))
  pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
  xp <- predict(mi1, pdat, type = 'lpmatrix')
  
  sprintf("Seascapenr%s",combi_width$V1[2])
  sprintf("Seascapenr%s",combi_width$V2[2])
  
  #for (i in 1:nrow(combi_width)){
  for (i in 1:nrow(combi_width)){
    print(i)
    pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
    xp <- predict(mi1, pdat, type = 'lpmatrix')
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
    Vb=vcov(mi1)
    dif <- X %*% coef(mi1)
    se <- sqrt(rowSums((X %*% Vb) * X))
    pred=cbind(dif,se)
    
    ## CI of the estimates
    set.seed(42)
    N <- 10000
    
    #predi <- predict(m1, pdat, se.fit = TRUE)
    #se.fit <- predi$se.fit
    
    BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)    
    Cg <- predict(mi1, pdat, type = "lpmatrix")
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
    
    DifS12plot=
      ggplot(S12pred, aes(x = Year, y = V1)) +
      theme_bw()+
      #geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
      geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
      geom_line(data=S12pred,aes(x=Year,y=V1),col="black",size=0.5) +
      geom_line(data=S12higher,aes(x=Year,y=V1),col="seagreen1",size=1)+
      geom_line(data=S12lower,aes(x=Year,y=V1),col="red",size=1)+
      ggtitle(sprintf("%s Seascape %s - %s",species_info$Scientific.name[j],combi_width$V1[i],combi_width$V2[i]))+
      labs(x = "Year", y = 'Difference Fish biomass index')+
      theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
      theme(axis.text.x = element_text(face="bold",size=9))+
      theme(axis.text.y = element_text(face="bold",size=9))+
      theme(axis.title.x = element_text(face="bold",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
      theme(axis.title.y = element_text(face="bold",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
      theme(axis.line = element_line(size=1))
    ggsave(DifS12plot,file = paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooth_%s_%s.png",species_info$file_name[j],combi_width$V1[i],combi_width$V2[i]),sep=""), height=3.5, width=4.5,dpi = 600)

  }
}  

#save different smooths in single large image            
for (j in 16:16){
  print(j)
  p=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth",species_info$file_name[j]),sep="")
  filenames=list.files(path=p,pattern="*.png",full.names=TRUE)
  
  rl <- lapply(filenames, png::readPNG)
  gl <- lapply(rl, grid::rasterGrob)
  cp=grid.arrange(gl[[1]],gl[[2]],gl[[3]],gl[[4]],gl[[5]],gl[[6]],
                  gl[[7]],gl[[8]],gl[[9]],gl[[10]],gl[[11]],gl[[12]],
                  gl[[13]],gl[[14]],gl[[15]],gl[[16]],gl[[17]],gl[[18]],
                  gl[[19]],gl[[20]],gl[[21]],gl[[22]],gl[[23]],gl[[24]],
                  gl[[25]],gl[[26]],gl[[27]],gl[[28]],gl[[29]],gl[[30]],
                  gl[[31]],gl[[32]],gl[[33]],gl[[34]],gl[[35]],gl[[36]],ncol=6,nrow=6)
  ggsave(cp,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooths.png",species_info$file_name[j]),sep=""),height=21, width=27,dpi = 300)
}




################################################################################
####GAM Anscombe transformed (sqrt(x+3/8)) and scaled by max (BETA errror)   ###-----
################################################################################
list_AIC_beta_2= vector("list",16)


for (j in 1:16){
  data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[j]),sep=""))
  print(species_info$file_name[j])
  
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  ssc_data=completeFun(data,"Seascapenr")
  ssc_data=subset(ssc_data,Quarter==1)
  data_agg=ssc_data
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species, data=ssc_data, sum)
  # 
   # #set outlier weight limit (95th percentile)
   #  nf_p=quantile(data_agg$Total_wgt,0.95)
   #  if (nf_p > 0){
   #   #subset of data excluding outlier weights
   #    data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
   #  }
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000
  
  #separate subset per seascape, anscombe transform, max scale it and finally rbind back into single dataset----
  #1----
  data_agg1=subset(data_agg,Seascapenr==1)
  data_agg1$ansc_biomass= nthroot(data_agg1$biomass_kg+(3/8),2)
  data_agg1$scaled_ansc_biomass=data_agg1$ansc_biomass/max(data_agg1$ansc_biomass)

  #2-----
  data_agg2=subset(data_agg,Seascapenr==2)
  data_agg2$ansc_biomass= nthroot(data_agg2$biomass_kg+(3/8),2)
  data_agg2$scaled_ansc_biomass=data_agg2$ansc_biomass/max(data_agg2$ansc_biomass)
  
  #3----
  data_agg3=subset(data_agg,Seascapenr==3)
  data_agg3$ansc_biomass= nthroot(data_agg3$biomass_kg+(3/8),2)
  data_agg3$scaled_ansc_biomass=data_agg3$ansc_biomass/max(data_agg3$ansc_biomass)
 
  #4----
  data_agg4=subset(data_agg,Seascapenr==4)
  data_agg4$ansc_biomass= nthroot(data_agg4$biomass_kg+(3/8),2)
  data_agg4$scaled_ansc_biomass=data_agg4$ansc_biomass/max(data_agg4$ansc_biomass)

  #5----
  data_agg5=subset(data_agg,Seascapenr==5)
  data_agg5$ansc_biomass= nthroot(data_agg5$biomass_kg+(3/8),2)
  data_agg5$scaled_ansc_biomass=data_agg5$ansc_biomass/max(data_agg5$ansc_biomass)
  
  #6----
  data_agg6=subset(data_agg,Seascapenr==6)
  data_agg6$ansc_biomass= nthroot(data_agg6$biomass_kg+(3/8),2)
  data_agg6$scaled_ansc_biomass=data_agg6$ansc_biomass/max(data_agg6$ansc_biomass)
  
  #7----
  data_agg7=subset(data_agg,Seascapenr==7)
  data_agg7$ansc_biomass= nthroot(data_agg7$biomass_kg+(3/8),2)
  data_agg7$scaled_ansc_biomass=data_agg7$ansc_biomass/max(data_agg7$ansc_biomass)
  
  #8----
  data_agg8=subset(data_agg,Seascapenr==8)
  data_agg8$ansc_biomass= nthroot(data_agg8$biomass_kg+(3/8),2)
  data_agg8$scaled_ansc_biomass=data_agg8$ansc_biomass/max(data_agg8$ansc_biomass)
  
  
  #9----
  data_agg9=subset(data_agg,Seascapenr==9)
  data_agg9$ansc_biomass= nthroot(data_agg9$biomass_kg+(3/8),2)
  data_agg9$scaled_ansc_biomass=data_agg9$ansc_biomass/max(data_agg9$ansc_biomass)
  
  #rbind----
  data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9)
  
  ####################################################
  #GAM model with separate smoother per seascape
  ###################################################
  print("gam")
  data_agg_ansc$Seascapenr=as.factor(data_agg_ansc$Seascapenr) #Set as factor
  data_agg_ansc = start_event(data_agg_ansc, column="Year",event="Seascapenr") #Sep time series per factor level
  
  
  mod=bam(scaled_ansc_biomass~Seascapenr+s(Year, by=Seascapenr),data=data_agg_ansc,method="REML",family=betar)
  print("gam fitted")
  print(AIC(mod))
  par(mfrow=c(2,2))
  gam.check(mod)
  #list_AIC_beta[[j]]=AIC(mod)
  
  dev.copy(png,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_diagnostics_2.png",species_info$file_name[j]),sep=""))
  dev.off()
  
  
  # ################################################################
  # # Plot predicted smooths with all data                         #----
  # ################################################################
  newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
  pred <- predict(mod, newd, se.fit = TRUE,type="response")
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


  img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[j]),sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  ### plot smooths with all data -----  
  s1plot=ggplot(data=S1,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg1)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S1")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  s1plot=s1plot+annotation_custom(g,xmin = min(data_agg1$Year), xmax = min(data_agg1$Year)+9, ymin = max(data_agg1$scaled_ansc_biomass)-0.15*(max(data_agg1$scaled_ansc_biomass)), ymax = max(data_agg1$scaled_ansc_biomass))

  
  s2plot=ggplot(data=S2,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S2")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s2plot=s2plot+annotation_custom(g,xmin = min(data_agg2$Year), xmax = min(data_agg2$Year)+8, ymin = max(data_agg2$scaled_ansc_biomass)-0.15*(max(data_agg2$scaled_ansc_biomass)), ymax = max(data_agg2$scaled_ansc_biomass))
  
  
  s3plot=ggplot(data=S3,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg3)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S3")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s3plot=s3plot+annotation_custom(g,xmin = min(data_agg3$Year), xmax = min(data_agg3$Year)+8, ymin = max(data_agg3$scaled_ansc_biomass)-0.15*(max(data_agg3$scaled_ansc_biomass)), ymax = max(data_agg3$scaled_ansc_biomass))
  
  s4plot=ggplot(data=S4,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg4)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S4")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s4plot=s4plot+annotation_custom(g,xmin = min(data_agg4$Year), xmax = min(data_agg4$Year)+8, ymin = max(data_agg4$scaled_ansc_biomass)-0.15*(max(data_agg4$scaled_ansc_biomass)), ymax = max(data_agg4$scaled_ansc_biomass))
  
  
  s5plot=ggplot(data=S5,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg5)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S5")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s5plot=s5plot+annotation_custom(g,xmin = min(data_agg5$Year), xmax = min(data_agg5$Year)+8, ymin = max(data_agg5$scaled_ansc_biomass)-0.15*(max(data_agg5$scaled_ansc_biomass)), ymax = max(data_agg5$scaled_ansc_biomass))

  
  s6plot=ggplot(data=S6,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg6)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S6")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s6plot=s6plot+annotation_custom(g,xmin = min(data_agg6$Year), xmax = min(data_agg6$Year)+8, ymin = max(data_agg6$scaled_ansc_biomass)-0.15*(max(data_agg6$scaled_ansc_biomass)), ymax = max(data_agg6$scaled_ansc_biomass))
  
  
  s7plot=ggplot(data=S7,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg7)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S7")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s7plot=s7plot+annotation_custom(g,xmin = min(data_agg7$Year), xmax = min(data_agg7$Year)+8, ymin = max(data_agg7$scaled_ansc_biomass)-0.15*(max(data_agg7$scaled_ansc_biomass)), ymax = max(data_agg7$scaled_ansc_biomass))
  
  s8plot=ggplot(data=S8,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg8)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S8")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s8plot=s8plot+annotation_custom(g,xmin = min(data_agg8$Year), xmax = min(data_agg8$Year)+8, ymin = max(data_agg8$scaled_ansc_biomass)-0.15*(max(data_agg8$scaled_ansc_biomass)), ymax = max(data_agg8$scaled_ansc_biomass))
  
  
  s9plot=ggplot(data=S9,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_point(aes(x=Year,y=scaled_ansc_biomass),alpha=0.2,colour="blue",data=data_agg9)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S9")+
    labs(x="Year",y="ANSC fish biomass trend")
  #Add species image
  #s9plot=s9plot+annotation_custom(g,xmin = min(data_agg9$Year), xmax = min(data_agg9$Year)+8, ymin = max(data_agg9$scaled_ansc_biomass)-0.15*(max(data_agg9$scaled_ansc_biomass)), ymax = max(data_agg9$scaled_ansc_biomass))
  
  
  cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,ncol=3,nrow=3)
  ggsave(cplot,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_data_smooths.png",species_info$file_name[j]),sep=""),height=7, width=9,dpi = 600)#


  ################################################################
  # Plot predicted smooths at smooth scale and simultaneous CI   #----
  ################################################################
  newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
  pred <- predict(mod, newd, se.fit = TRUE,type="response")
  se.fit <- pred$se.fit
  pred <- transform(cbind(data.frame(pred), newd),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit))
  ###
  Vb <- vcov(mod)
  newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
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
  S3=pred[which(pred$Seascapenr==3),]
  S4=pred[which(pred$Seascapenr==4),]
  S5=pred[which(pred$Seascapenr==5),]
  S6=pred[which(pred$Seascapenr==6),]
  S7=pred[which(pred$Seascapenr==7),]
  S8=pred[which(pred$Seascapenr==8),]
  S9=pred[which(pred$Seascapenr==9),]

  ### plot smooths with all data -----
  s1plot=ggplot(data=S1,aes(x=Year))+
    theme_bw()+
    geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2)+
    geom_ribbon(aes(ymin=lwrP,ymax=uprP),alpha=0.2)+
    geom_line(aes(y=fit),col="black",size=1)+
    ggtitle("S1")+
    labs(x="Year",y="ANSC fish biomass trend")

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

  cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,ncol=3,nrow=3)
  ggsave(cplot,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_smooths.png",species_info$file_name[j]),sep=""),height=7, width=9,dpi = 600)#

  #save residual autocorrelation plot
  par(mfrow=c(1,2))
  acf(residuals(mod),main="ACF(mod)")
  pacf(residuals(mod),main="partial ACF(mod)")
  dev.copy(png,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_autocorrelation.png",species_info$file_name[j]),sep=""))
  dev.off()

  #####################################################
  # Compute difference smooths using simultaneous CI's #  ----
  #####################################################
  combi=c(1,2,3,4,5,6,7,8,9)
  combi_width=as.data.frame(combinations(n=9,r=2,v=combi))
  pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
  xp <- predict(mod, pdat, type = 'lpmatrix')

  sprintf("Seascapenr%s",combi_width$V1[2])
  sprintf("Seascapenr%s",combi_width$V2[2])

  for (i in 1:nrow(combi_width)){
    print(i)
    pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
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
    ggsave(DifS12plot,file = paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooth_%s_%s_ANSC.png",species_info$file_name[j],combi_width$V1[i],combi_width$V2[i]),sep=""), height=3.5, width=4.5,dpi = 600)

  }
}


#save different smooths in single large image            
for (j in 1:1){
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




################################################################################
####GAM Anscombe transformed (sqrt(x+3/8)) and scaled by max (BETA errror)   ###-----
################################################################################
list_AIC_beta_2= vector("list",16)


for (j in 1:16){
  data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[1]),sep=""))
  print(species_info$file_name[j])
  
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  ssc_data=completeFun(data,"Seascapenr")
  ssc_data=subset(ssc_data,Quarter==1)
  data_agg=ssc_data
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species, data=ssc_data, sum)
  # 
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
  
  #2-----
  data_agg2=subset(data_agg,Seascapenr==2)
  data_agg2$ansc_biomass= nthroot(data_agg2$biomass_kg+(3/8),2)
  data_agg2$scaled_ansc_biomass=data_agg2$ansc_biomass/max(data_agg2$ansc_biomass)
  data_agg2$scaled_ansc_biomass=scaleTR(data_agg2$ansc_biomass)
  
  #3----
  data_agg3=subset(data_agg,Seascapenr==3)
  data_agg3$ansc_biomass= nthroot(data_agg3$biomass_kg+(3/8),2)
  data_agg3$scaled_ansc_biomass=data_agg3$ansc_biomass/max(data_agg3$ansc_biomass)
  data_agg3$scaled_ansc_biomass=scaleTR(data_agg3$ansc_biomass)
  
  #4----
  data_agg4=subset(data_agg,Seascapenr==4)
  data_agg4$ansc_biomass= nthroot(data_agg4$biomass_kg+(3/8),2)
  data_agg4$scaled_ansc_biomass=data_agg4$ansc_biomass/max(data_agg4$ansc_biomass)
  data_agg4$scaled_ansc_biomass=scaleTR(data_agg4$ansc_biomass)
  
  #5----
  data_agg5=subset(data_agg,Seascapenr==5)
  data_agg5$ansc_biomass= nthroot(data_agg5$biomass_kg+(3/8),2)
  data_agg5$scaled_ansc_biomass=data_agg5$ansc_biomass/max(data_agg5$ansc_biomass)
  data_agg5$scaled_ansc_biomass=scaleTR(data_agg5$ansc_biomass)
  
  #6----
  data_agg6=subset(data_agg,Seascapenr==6)
  data_agg6$ansc_biomass= nthroot(data_agg6$biomass_kg+(3/8),2)
  data_agg6$scaled_ansc_biomass=data_agg6$ansc_biomass/max(data_agg6$ansc_biomass)
  data_agg6$scaled_ansc_biomass=scaleTR(data_agg6$ansc_biomass)
  
  #7----
  data_agg7=subset(data_agg,Seascapenr==7)
  data_agg7$ansc_biomass= nthroot(data_agg7$biomass_kg+(3/8),2)
  data_agg7$scaled_ansc_biomass=data_agg7$ansc_biomass/max(data_agg7$ansc_biomass)
  data_agg7$scaled_ansc_biomass=scaleTR(data_agg7$ansc_biomass)
  
  #8----
  data_agg8=subset(data_agg,Seascapenr==8)
  data_agg8$ansc_biomass= nthroot(data_agg8$biomass_kg+(3/8),2)
  data_agg8$scaled_ansc_biomass=data_agg8$ansc_biomass/max(data_agg8$ansc_biomass)
  data_agg8$scaled_ansc_biomass=scaleTR(data_agg8$ansc_biomass)
  
  
  #9----
  data_agg9=subset(data_agg,Seascapenr==9)
  data_agg9$ansc_biomass= nthroot(data_agg9$biomass_kg+(3/8),2)
  data_agg9$scaled_ansc_biomass=data_agg9$ansc_biomass/max(data_agg9$ansc_biomass)
  data_agg9$scaled_ansc_biomass=scaleTR(data_agg9$ansc_biomass)
  
  #rbind----
  data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9)
  
  ####################################################
  #GAM model with separate smoother per seascape
  ###################################################
  print("gam")
  data_agg_ansc$Seascapenr=as.factor(data_agg_ansc$Seascapenr) #Set as factor
  data_agg_ansc = start_event(data_agg_ansc, column="Year",event="Seascapenr") #Sep time series per factor level
  
  
  mod=bam(scaled_ansc_biomass~Seascapenr+s(Year, by=Seascapenr),data=data_agg_ansc,method="REML",family=betar)
  print("gam fitted")
  print(AIC(mod))
  par(mfrow=c(2,2))
  gam.check(mod)
  #list_AIC_beta[[j]]=AIC(mod)
  
  dev.copy(png,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_diagnostics_3.png",species_info$file_name[j]),sep=""))
  dev.off()
  
  
  # # ################################################################
  # # # Plot predicted smooths with all data                         #----
  # # ################################################################
  # newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
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
  # 
  # 
  img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[j]),sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  # 
  # ### plot smooths with all data -----  
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
  # 
  # cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,ncol=3,nrow=3)
  # ggsave(cplot,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_data_smooths_4.png",species_info$file_name[j]),sep=""),height=7, width=9,dpi = 600)#
  # 

  ################################################################
  # Plot predicted smooths at smooth scale and simultaneous CI   #----
  ################################################################
  newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
  pred <- predict(mod, newd, se.fit = TRUE)
  se.fit <- pred$se.fit
  pred <- transform(cbind(data.frame(pred), newd),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit))
  ###
  Vb <- vcov(mod)
  newd <- expand.grid(Year = seq(1977, 2019, length = 126), Seascapenr = c(1,2,3,4,5,6,7,8,9))
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
  S3=pred[which(pred$Seascapenr==3),]
  S4=pred[which(pred$Seascapenr==4),]
  S5=pred[which(pred$Seascapenr==5),]
  S6=pred[which(pred$Seascapenr==6),]
  S7=pred[which(pred$Seascapenr==7),]
  S8=pred[which(pred$Seascapenr==8),]
  S9=pred[which(pred$Seascapenr==9),]

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

  cplot=grid.arrange(s1plot,s2plot,s3plot,s4plot,s5plot,s6plot,s7plot,s8plot,s9plot,ncol=3,nrow=3)
  ggsave(cplot,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_smooths_3.png",species_info$file_name[j]),sep=""),height=7, width=9,dpi = 600)#

  #save residual autocorrelation plot
  par(mfrow=c(1,2))
  acf(residuals(mod),main="ACF(mod)")
  pacf(residuals(mod),main="partial ACF(mod)")
  dev.copy(png,file=paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/GAM_ANSC_autocorrelation_3.png",species_info$file_name[j]),sep=""))
  dev.off()

  # # #####################################################
  # # # Compute difference smooths using simultaneous CI's #  ----
  # # #####################################################
  # combi=c(1,2,3,4,5,6,7,8,9)
  # combi_width=as.data.frame(combinations(n=9,r=2,v=combi))
  # pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
  # xp <- predict(mod, pdat, type = 'lpmatrix')
  # 
  # sprintf("Seascapenr%s",combi_width$V1[2])
  # sprintf("Seascapenr%s",combi_width$V2[2])
  # 
  # for (i in 1:nrow(combi_width)){
  #   print(i)
  #   pdat=expand.grid(Year=seq(1977,2019, length=126),Seascapenr = c(1,2,3,4,5,6,7,8,9))
  #   xp <- predict(mod, pdat, type = 'lpmatrix')
  #   col1 <- grepl(sprintf("Seascapenr%s",combi_width$V1[i]), colnames(xp))
  #   col2 <- grepl(sprintf("Seascapenr%s",combi_width$V2[i]), colnames(xp))
  #   row1 <- pdat[["Seascapenr"]] == combi_width$V1[i]
  #   row2 <- pdat[["Seascapenr"]] == combi_width$V2[i]
  # 
  #   ## difference rows of pred for data from comparison
  #   X <- xp[row1, ] - xp[row2, ]
  #   ## zero out cols of X related to splines for other lochs
  #   X[, ! (col1 | col2)] <- 0
  #   ## zero out the parametric cols
  #   X[, !grepl('^s\\(', colnames(xp))] <- 0
  # 
  #   #predicted values for differences
  #   Vb=vcov(mod)
  #   dif <- X %*% coef(mod)
  #   se <- sqrt(rowSums((X %*% Vb) * X))
  #   pred=cbind(dif,se)
  # 
  #   ## CI of the estimates
  #   set.seed(42)
  #   N <- 10000
  # 
  #   #predi <- predict(m1, pdat, se.fit = TRUE)
  #   #se.fit <- predi$se.fit
  # 
  #   BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  #   Cg <- predict(mod, pdat, type = "lpmatrix")
  #   simDev <- Cg %*% t(BUdiff)
  # 
  #   absDev <- abs(sweep(simDev, 1, se, FUN = "/"))
  #   masd <- apply(absDev, 2L, max)
  #   crit <- quantile(masd, prob = 0.95, type = 8)
  # 
  #   newd <- expand.grid(Year = seq(1977, 2019, length = 126))
  # 
  #   S12pred <- transform(cbind(data.frame(pred), newd),
  #                        uprP = V1 + (2 * se),
  #                        lwrP = V1 - (2 * se),
  #                        uprS = V1 + (crit * se),
  #                        lwrS = V1 - (crit * se))
  #   S12pred$sigdif <-ifelse(S12pred$lwrS<0 & S12pred$uprS<0,"lower",ifelse(S12pred$lwrS>0 & S12pred$uprS>0,"higher","no"))
  #   S12higher=S12pred
  #   S12higher["V1"][S12pred["sigdif"]!="higher"]=NA
  #   S12lower=S12pred
  #   S12lower["V1"][S12pred["sigdif"]!="lower"]=NA
  # 
  #   DifS12plot=
  #     ggplot(S12pred, aes(x = Year, y = V1)) +
  #     theme_bw()+
  #     #geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
  #     geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  #     geom_line(data=S12pred,aes(x=Year,y=V1),col="black",size=0.5) +
  #     geom_line(data=S12higher,aes(x=Year,y=V1),col="seagreen1",size=1)+
  #     geom_line(data=S12lower,aes(x=Year,y=V1),col="red",size=1)+
  #     ggtitle(sprintf("%s Seascape %s - %s",species_info$Scientific.name[j],combi_width$V1[i],combi_width$V2[i]))+
  #     labs(x = "Year", y = 'Difference Fish biomass trend')+
  #     theme(plot.title = element_text(hjust = 0.25,size=12,face="bold"))+
  #     theme(axis.text.x = element_text(face="bold",size=9))+
  #     theme(axis.text.y = element_text(face="bold",size=9))+
  #     theme(axis.title.x = element_text(face="bold",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  #     theme(axis.title.y = element_text(face="bold",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  #     theme(axis.line = element_line(size=1))
  #   ggsave(DifS12plot,file = paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/Diff_smooth/Diff_smooth_%s_%s_ANSC.png",species_info$file_name[j],combi_width$V1[i],combi_width$V2[i]),sep=""), height=3.5, width=4.5,dpi = 600)
  # 
  # }
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




##############################
#### temporal correlation #### ----
##############################
for (j in 2:16){
  data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[j]),sep=""))
  print(species_info$file_name[j])
  
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  ssc_data=completeFun(data,"Seascapenr")
  ssc_data=subset(ssc_data,Quarter==1)
  data_agg=ssc_data
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species, data=ssc_data, sum)
  # 
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
  
  #2-----
  data_agg2=subset(data_agg,Seascapenr==2)
  data_agg2$ansc_biomass= nthroot(data_agg2$biomass_kg+(3/8),2)
  data_agg2$scaled_ansc_biomass=data_agg2$ansc_biomass/max(data_agg2$ansc_biomass)
  data_agg2$scaled_ansc_biomass=scaleTR(data_agg2$ansc_biomass)
  
  #3----
  data_agg3=subset(data_agg,Seascapenr==3)
  data_agg3$ansc_biomass= nthroot(data_agg3$biomass_kg+(3/8),2)
  data_agg3$scaled_ansc_biomass=data_agg3$ansc_biomass/max(data_agg3$ansc_biomass)
  data_agg3$scaled_ansc_biomass=scaleTR(data_agg3$ansc_biomass)
  
  #4----
  data_agg4=subset(data_agg,Seascapenr==4)
  data_agg4$ansc_biomass= nthroot(data_agg4$biomass_kg+(3/8),2)
  data_agg4$scaled_ansc_biomass=data_agg4$ansc_biomass/max(data_agg4$ansc_biomass)
  data_agg4$scaled_ansc_biomass=scaleTR(data_agg4$ansc_biomass)
  
  #5----
  data_agg5=subset(data_agg,Seascapenr==5)
  data_agg5$ansc_biomass= nthroot(data_agg5$biomass_kg+(3/8),2)
  data_agg5$scaled_ansc_biomass=data_agg5$ansc_biomass/max(data_agg5$ansc_biomass)
  data_agg5$scaled_ansc_biomass=scaleTR(data_agg5$ansc_biomass)
  
  #6----
  data_agg6=subset(data_agg,Seascapenr==6)
  data_agg6$ansc_biomass= nthroot(data_agg6$biomass_kg+(3/8),2)
  data_agg6$scaled_ansc_biomass=data_agg6$ansc_biomass/max(data_agg6$ansc_biomass)
  data_agg6$scaled_ansc_biomass=scaleTR(data_agg6$ansc_biomass)
  
  #7----
  data_agg7=subset(data_agg,Seascapenr==7)
  data_agg7$ansc_biomass= nthroot(data_agg7$biomass_kg+(3/8),2)
  data_agg7$scaled_ansc_biomass=data_agg7$ansc_biomass/max(data_agg7$ansc_biomass)
  data_agg7$scaled_ansc_biomass=scaleTR(data_agg7$ansc_biomass)
  
  #8----
  data_agg8=subset(data_agg,Seascapenr==8)
  data_agg8$ansc_biomass= nthroot(data_agg8$biomass_kg+(3/8),2)
  data_agg8$scaled_ansc_biomass=data_agg8$ansc_biomass/max(data_agg8$ansc_biomass)
  data_agg8$scaled_ansc_biomass=scaleTR(data_agg8$ansc_biomass)
  
  
  #9----
  data_agg9=subset(data_agg,Seascapenr==9)
  data_agg9$ansc_biomass= nthroot(data_agg9$biomass_kg+(3/8),2)
  data_agg9$scaled_ansc_biomass=data_agg9$ansc_biomass/max(data_agg9$ansc_biomass)
  data_agg9$scaled_ansc_biomass=scaleTR(data_agg9$ansc_biomass)
  
  #rbind----
  data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9)
  
  ####################################################
  #GAM model with separate smoother per seascape
  ###################################################
  print("gam")
  data_agg_ansc$Seascapenr=as.factor(data_agg_ansc$Seascapenr) #Set as factor
  data_agg_ansc = start_event(data_agg_ansc, column="Year",event="Seascapenr") #Sep time series per factor level
  
  
  mod=bam(scaled_ansc_biomass~Seascapenr+s(Year, by=Seascapenr),data=data_agg_ansc,method="REML",family=betar)
  
  print("gam fitted")
  
  combi=c(1,2,3,4,5,6,7,8,9)
  combi_width=as.data.frame(combinations(n=9,r=2,v=combi))
  
  newdat = as.data.frame(1977:2019)
  colnames(newdat)="Year"
  
  #1000 simulations from model posterior distribution
  sim_dat=as.data.frame(smooth_samples(mod,n=1000,new_data=newdat,n_vals=43))
  sim_dat$Year=floor(sim_dat$.x1)
  plot(x=sim_dat$Year,y=sim_dat$value)
  
  #go iteratively through smooth pair combinations and simulate autocorrelation
  sim_autocor=setNames(data.frame(matrix(ncol=6,nrow=0)),c("Seascapes","lag-year","lag","cor","lwr_conf","upr_conf"))
  for (i in 1:nrow(combi_width)){
    print(i)
    V1= sprintf("s(Year):Seascapenr%s",combi_width$V1[i])
    V2= sprintf("s(Year):Seascapenr%s",combi_width$V2[i])
    Seascapes=sprintf("%s & %s",combi_width$V1[i],combi_width$V2[i])
    
    sim_dat_V1=subset(sim_dat,term==V1)
    sim_dat_V2=subset(sim_dat,term==V2)
  
    par(mfrow=c(2,2))
    #plot(sim_dat_V1$.x1,sim_dat_V1$value)
    #plot(sim_dat_V2$.x1,sim_dat_V2$value)
    
    #test correlation between V1 & V2 at lag-1 using all 1000 simulations
    
    for (lag in 1:10){
      for (y in 1977:2018){
        if (y+lag <= 2019){
          dat_V1=subset(sim_dat_V1,Year==y)
          dat_V2=subset(sim_dat_V2,Year==y+lag)
          lag_year=sprintf("%s - %s", y,y+lag)
          
          cor=cor.test(dat_V1$value,dat_V2$value)
          est_cor=cor$estimate
          lwr_conf=cor$conf.int[2]
          upr_conf=cor$conf.int[1]
          sim_autocor[nrow(sim_autocor)+1,]=list(Seascapes,lag_year,lag,est_cor,lwr_conf,upr_conf)
          }
        }
      }
    
    # ccf(sim_dat_V1$value,sim_dat_V2$value)
    # autocor_list=list()
    # sim_autocor=setNames(data.frame(matrix(ncol=21,nrow=1000)),c(-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10))
    # 
    # #run 1000 cross correlations
    # for (k in 1:10){
    #   print(k)
    #   draw_n_V1=subset(sim_dat_V1,draw==k)
    #   draw_n_V2=subset(sim_dat_V2,draw==k)
    #   
    #   ts_V1=ts(draw_n_V1$value)
    #   ts_V2=ts(draw_n_V2$value)
    #   
    #   par(mfrow=c(2,2))
    #   plot(ts_V1)
    #   plot(ts_V2)
    #   plot(diff(ts_V1))
    #   plot(diff(ts_V2))
    #   
    #   #decompose_V1=stl(ts_V1,"multiplicative")
    #   #decompose_V2=decompose(ts_V2)
    #   #diff_draw_V1=diff(draw_n_V1$value,lag=1,difference=1)
    #   #diff_draw_V2=diff(draw_n_V2$value,lag=1,difference=1)
    #   
    #   acf_plot(draw_n_V1$value,max_lag=40)
    #   acf_plot(draw_n_V2$value,max_lag=40)
    #   acf(dv1$x)
    #   acf(dv2$x)
    #   
    # 
    #   ccfV12=ccf(diff(ts_V1),diff(ts_V2),pl=TRUE,lag.max = 10)
    #   sim_autocor[k,1:21]=ccfV12$acf[1:21] #add acf values to dataframe
    # }
    # 
    # #transpose dataframe
    # sim_autocor=as.data.frame(t(sim_autocor))
    # sim_autocor=tibble::rownames_to_column(sim_autocor,"lag")
    # 
    # autocor_unc=setNames(data.frame(matrix(ncol=5,nrow=21)),c("lag","mean","se","uprP","lwrP"))
    # 
    # autocor_unc$lag=as.numeric(sim_autocor$lag)
    # autocor_unc$mean=rowMeans(sim_autocor[,2:1001])
    # autocor_unc$se = apply(sim_autocor[,2:1001],1,se)
    # autocor_unc$uprP = autocor_unc$mean + (autocor_unc$se*2)
    # autocor_unc$lwrP = autocor_unc$mean - (autocor_unc$se*2)
    # 
    # autocor_plot=ggplot(autocor_unc, aes(x = lag, y = mean))+
    #              theme_classic()+
    #              geom_bar(stat="identity",col="black",fill="gray",width=0.5)+
    #              geom_errorbar(aes(ymin=lwrP,ymax=uprP),width=0.2)+
    #                         
    #              geom_hline(yintercept=0.3, linetype="dashed", color = "red",size=1)+
    #              geom_hline(yintercept=-0.3, linetype="dashed", color = "red",size=1)+
    #              geom_hline(yintercept=0, color = "black",size=1)+
    #             
    #              
    #              ggtitle(sprintf("%s cross-correlation S%s & S%s",species_info$Scientific.name[j],combi_width$V1[i],combi_width$V2[i]))+
    #              labs(x = "Lag (Year)", y = 'Autocorrelation factor')+
    #              scale_x_continuous(breaks = round(seq(min(autocor_unc$lag), max(autocor_unc$lag), by = 1),1))+
    #              scale_y_continuous(breaks = round(seq(min(autocor_unc$mean), max(autocor_unc$mean), by = 0.1),1))+
    #   
    #               theme(plot.title = element_text(hjust = 0.25,size=11,face="bold"))+
    #               theme(axis.text.x = element_text(face="bold",size=9))+
    #               theme(axis.text.y = element_text(face="bold",size=9))+
    #               theme(axis.title.x = element_text(face="bold",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)))+
    #               theme(axis.title.y = element_text(face="bold",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)))
    # 
    # #autocor_plot
    # ggsave(autocor_plot,file = paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/CCF/CCF_smooth_%s_%s_ANSC.png",species_info$file_name[j],combi_width$V1[i],combi_width$V2[i]),sep=""), height=3.5, width=4.5,dpi = 600)
    
  }
    #write autocor file to csv
  write.csv(sim_autocor,paste(path,sprintf("graphs/Seascapes/Q1/%s/GAM/CCF/autocor_lag.csv",species_info$file_name[j]),sep=""))
}
#model=nlin_causality.test(draw_n_V1$value,draw_n_V2$value,lag=2, c(2, 10), c(4, 20), 500, TRUE)
#model$summary()






















####Q1_Q3 Separate plots for all 9 seascapes all data ----
for (i in 1:12){
  print(paste(i,species_info$file_name[i],sep=","))
  # Read in data
  data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[i]),sep=""))
  
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  for(j in 1:9){
    print(paste("Seascape",j,sep=" "))
    ssc_data=subset(data,Seascapenr==j)
    ssc_data=subset(ssc_data, Year >= 1991)
    
    data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr, data=ssc_data, sum)
    
    #set outlier weight limit (95th percentile)
    nf_p=quantile(data_agg$Total_wgt,0.95) 
    
    #subset of data excluding outlier weights
    #data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
    
    if (nf_p > 0){
      #subset of data excluding outlier weights
      data_agg=subset(data_agg,Total_wgt < nf_p[[1]])
    }
    
    
    # convert weights from grams to kg
    data_agg$biomass_kg=data_agg$Total_wgt/1000 
    
    #mean per year to get better overview
    data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
    
    #Plot
    p=ggplot(data_agg2,aes(Year,biomass_kg))+
      geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
      geom_point(size=3, alpha=0.5)+
      theme_classic()+
      ggtitle(sprintf("Q1 & Q3 %s Seascape %s",species_info$Scientific.name[i],j))+
      labs(x = "Year", y = 'Mean CPUE (kg)')+
      theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"))+
      theme(axis.text.x = element_text(face="bold",size=18))+
      theme(axis.text.y = element_text(face="bold",size=18))+
      theme(axis.title.x = element_text(face="bold",size=20,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
      theme(axis.title.y = element_text(face="bold",size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
      theme(axis.line = element_line(size=1))
    
    #Add species image
    img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[i]),sep=""))
    g <- rasterGrob(img, interpolate=TRUE)
    
    p2=p+annotation_custom(g,xmin = 1991, xmax = 1996, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg)+0.1*max(data_agg2$biomass_kg))
    ggsave(p2, file = paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_%s_Trend_%s.png",species_info$file_name[i],j,species_info$file_name[i]),sep=""), height=7, width=9,dpi = 600) # adjust dpi )
    
  }
  
  #stack in single image
  plot1 = readPNG(paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_%s_Trend_%s.png",species_info$file_name[i],1,species_info$file_name[i]),sep=""))
  plot2 = readPNG(paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_%s_Trend_%s.png",species_info$file_name[i],2,species_info$file_name[i]),sep=""))
  plot3 = readPNG(paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_%s_Trend_%s.png",species_info$file_name[i],3,species_info$file_name[i]),sep=""))
  plot4 = readPNG(paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_%s_Trend_%s.png",species_info$file_name[i],4,species_info$file_name[i]),sep=""))
  plot5 = readPNG(paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_%s_Trend_%s.png",species_info$file_name[i],5,species_info$file_name[i]),sep=""))
  plot6 = readPNG(paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_%s_Trend_%s.png",species_info$file_name[i],6,species_info$file_name[i]),sep=""))
  plot7 = readPNG(paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_%s_Trend_%s.png",species_info$file_name[i],7,species_info$file_name[i]),sep=""))
  plot8 = readPNG(paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_%s_Trend_%s.png",species_info$file_name[i],8,species_info$file_name[i]),sep=""))
  plot9 = readPNG(paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_%s_Trend_%s.png",species_info$file_name[i],9,species_info$file_name[i]),sep=""))
  
  cplot=grid.arrange(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot5),
                     rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),rasterGrob(plot9),ncol=3)
  ggsave(cplot, file = paste(path,sprintf("graphs/Seascapes/Q1_Q3/%s/Q1_Q3_compiled_Trend_%s.png",species_info$file_name[i],species_info$file_name[i]),sep=""), height=21, width=27,dpi = 600) # adjust dpi )
  
}
























plot.gam(mi1)
plot_gam_check(mi1)


###MANUAL TRIAL 1 SPECIES LIMANDA LIMANDA----
#work with filtered datasets #1 Dab - Limanda limanda
data=read.csv(paste(path,sprintf("filtered_%s.csv",species_info$file_name[1]),sep=""))

#data%>%mutate(ID = group_indices_(data, .dots=c("Year", "Ship","HaulNo"))) 
data$ID = cumsum(!duplicated(data[3:7])) # Unique combination between vessel, haulNo. and Year, i.e. each unique haul

data_agg=aggregate(Total_wgt ~ ID+Year, data=data, sum) #aggregate weight over each haul (sum)

nf_p=quantile(data_agg$Total_wgt,0.95) # set outlier weight limit (95th percentile)
print(nf_p[[1]])

data_agg=subset(data_agg,Total_wgt < nf_p[[1]]) #subset of data excluding outlier weights
data_agg$biomass_kg=data_agg$Total_wgt/1000 # convert weights from grams to kg

#Plot of biomass of Limanda limanda per haul per year (in kg)
plot(x=data_agg$Year, y=data_agg$biomass_kg)

#mean per year to get better overview
data_agg2=aggregate(biomass_kg~Year,data=data_agg,mean)
plot(x=data_agg2$Year,y=data_agg2$biomass_kg)


#nicer ggplot with smoother and CI's
p=ggplot(data_agg,aes(Year,biomass_kg))+
  geom_point()+
  geom_smooth(method="auto",se=TRUE)
p


p2=ggplot(data_agg2,aes(Year,biomass_kg))+
  geom_smooth(stat = 'smooth', color = 'blue', method = 'gam',size=1.5, formula = y ~ s(x, bs = "cs"))+
  geom_point(size=3, alpha=0.5)+
  theme_classic()+
  ggtitle("Limanda limanda")+
  labs(x = "Year", y = 'CPUE (kg)')+
  theme(plot.title = element_text(hjust = 0.5,size=18,face="bold"))+
  theme(axis.text.x = element_text(face="bold",size=12))+
  theme(axis.text.y = element_text(face="bold",size=12))+
  theme(axis.title.x = element_text(face="bold",size=14))+
  theme(axis.title.y = element_text(face="bold",size=14))+
  theme(axis.line = element_line(size=1))
p2

#Add species image
img=readPNG(source=paste(path,sprintf("graphs/spec_%s.png",species_info$file_name[1]),sep=""))
g <- rasterGrob(img, interpolate=TRUE)

p3=p2+annotation_custom(g,xmin = 1977, xmax = 1985, ymin = max(data_agg2$biomass_kg)-0.15*(max(data_agg2$biomass_kg)), ymax = max(data_agg2$biomass_kg))
p3
ggsave(p3, file = paste(path,sprintf("graphs/%s.png",species_info$file_name[1]),sep=""), dpi = 600) # adjust dpi )






