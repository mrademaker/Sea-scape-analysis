library(devtools)
library(windowscanr)
library(rollRegres)
###############################################
# Time lagged DCCA_cc coefficient calculation #
###############################################

## data_1
x= c(-1.042061,-0.669056,-0.685977,-0.067925,0.808380,1.385235,1.455245,0.540762 ,0.139570,-1.038133,0.080121,-0.102159,-0.068675,0.515445,0.600459,0.655325,0.610604,0.482337,0.079108,-0.118951,-0.050178,0.007500,-0.200622)
## data_2
y= c(-2.368030,-2.607095,-1.277660,0.301499,1.346982,1.885968,1.765950,1.242890,-0.464786,0.186658,-0.036450,-0.396513,-0.157115,-0.012962,0.378752,-0.151658,0.774253,0.646541,0.311877,-0.694177,-0.412918,-0.338630,0.276635)
## window size = 6
k=5

DCCA_CC=function(x,y,k){
  
  ## calculate cumulated deviation time series Xt
  xx<- cumsum(x - mean(x))
  yy<- cumsum(y - mean(y))
  
  
  data=data.frame(time=seq(1:length(x)))
  data$time=as.numeric(seq(1:length(x)))
  
  data$xx=xx
  data$yy=yy
  
  par(mfrow=c(2,2))
  plot(x=data$time,y=data$xx)
  plot(x=data$time,y=data$yy)
  
  store_f2_x=list()
  store_f2_y=list()
  store_f2_xy=list()
  
  lag_list=lags
  dcca_list=list()
    #### equation 2. fit linear model on sliding window
    for (i in 1:length(xx)){
      # so sliding window step size 1 untill i + window size reaches end of time series
      if (i + k < length(xx)){
        print(i)
        dat=data[i:(i+k),] #subset of x data in window
        
        #for x
        x_act=dat$xx                  #data values in window
        mod_x=lm(xx~time,data=dat)    #fit linear model
        x_hat = mod_x$fitted.values     #fitted values x in window
        f2DCCA_x = sum(((x_act-x_hat)^2))/(k+1)    #equation 2 
        store_f2_x = as.numeric(append(store_f2_x,f2DCCA_x))#add to list
        
        #for y
        y_act=dat$yy
        mod_y=lm(yy~time,data=dat)
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
    pDCCA = F_DCCA_xy/(sqrt(F_DCCA_x)*sqrt(F_DCCA_y)) #equation 1.
  return(pDCCA)
}

DCCA_CC(x=x,y=y,k=5)
