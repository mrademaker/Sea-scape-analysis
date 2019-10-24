library(devtools)
library(windowscanr)
library(rollRegres)
library(NlinTS)
###############################################
# Time lagged DCCA_cc coefficient calculation #
###############################################

## data_1
x= c(-1.042061,-0.669056,-0.685977,-0.067925,0.808380,1.385235,1.455245,0.540762 ,0.139570,-1.038133,0.080121,-0.102159,-0.068675,0.515445,0.600459,0.655325,0.610604,0.482337,0.079108,-0.118951,-0.050178,0.007500,-0.200622)
## data_2
y= c(-2.368030,-2.607095,-1.277660,0.301499,1.346982,1.885968,1.765950,1.242890,-0.464786,0.186658,-0.036450,-0.396513,-0.157115,-0.012962,0.378752,-0.151658,0.774253,0.646541,0.311877,-0.694177,-0.412918,-0.338630,0.276635)
## window size = 6


dat=nlin_causality.test(x,y,lag=1,LayersUniv=c(2,2),LayersBiv=c(4,4),500,TRUE)
dat$summary()

L_AXHA=function(x,y,L,lags){
  
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
        fqXY=sign(XY)*abs(XY)
        
        fqXX = (xx[i]-xx[i+L])^2
        fqYY = (yy[i]-yy[i+L])^2
        
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
  LAXHA=as.data.frame(AXHA_list)
  LAXHA$lag=(0:(length(AXHA_list)-1))
  names(LAXHA)[names(pqLAXHA)=="AXHA_list"]="corr"
  return(LAXHA)
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

df=q_L_AXHA(x=x,y=y,L=1,q=2,lags=c(0,1,2,3,4,5,6,7,8,9,10))
barplot(df$corr,xlab="lag",ylab="correlation")

df=L_AXHA(x=x,y=y,L=1,lags=c(0,1,2,3,4,5,6,7,8,9,10))
