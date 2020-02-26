library(DataCombine)
library(ggplot2)

#set path----
path=("C:/Users/mrademaker/Documents/Research projects/STCNWS/DATRAS/cpue_length_hour/")

#file with species info----
species_info=read.table(paste(path,"Data/spec_info.txt",sep=""),sep="\t",header=TRUE)
species_info$file_name =gsub(" ","_",species_info$Scientific.name)

datpd=read.csv(paste(path,sprintf("Seascapes/Q1/%s/GAM/Diff_smooth/num_diff.csv",species_info$file_name[1]),sep=""))
datmi=read.csv(paste(path,sprintf("Seascapes/Q1/%s/GAM/MI/mutinf_lag.csv",species_info$file_name[1]),sep=""))
#datmi=subset(datmi,lag==0)

combis=c("1 - 2"  ,"1 - 3" ,"2 - 3" ,"2 - 4" ,"3 - 5" ,"3 - 6" ,"4 - 5" ,"5 - 10" ,"5 - 6"  ,"5 - 7" ,"5 - 8" ,"5 - 9" ,"6 - 8" ,"7 - 8" ,"7 - 9" ,"8 - 10" ,"8 - 9" ,"9 - 10")
#row510=datpd[35,]
datpd <- datpd[datpd$combi %in% combis, ]
#datpd=datpd[-c(12),]
#datpd=InsertRow(datpd, NewRow = row510, RowNum = 8)


cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)


#BCMI neigbouring seascapes
for (i in c(1,2,3,4,5,7,12,13,14)){
  datmi=read.csv(paste(path,sprintf("Seascapes/Q1/%s/GAM/MI/mutinf_lag.csv",species_info$file_name[i]),sep=""))
  datmi$Seascapes=gsub("&","-",datmi$Seascapes) 
  print(species_info$file_name[i])
  datmi$Seascapes=factor(datmi$Seascapes, levels = c("1 - 2", "1 - 3", "1 - 4", "1 - 5", "1 - 6", "1 - 7", "1 - 8", "1 - 9", "1 - 10",
                                                   "2 - 3", "2 - 4", "2 - 5", "2 - 6", "2 - 7", "2 - 8", "2 - 9", "2 - 10", "3 - 4",
                                                   "3 - 5", "3 - 6", "3 - 7", "3 - 8", "3 - 9", "3 - 10", "4 - 5", "4 - 6", "4 - 7",
                                                   "4 - 8", "4 - 9", "4 - 10", "5 - 6","5 - 7", "5 - 8", "5 - 9", "5 - 10", "6 - 7",
                                                   "6 - 8", "6 - 9", "6 - 10", "7 - 8", "7 - 9", "7 - 10", "8 - 9", "8 - 10", "9 - 10"))
  #row510=datmi[35,]
  datmi <- datmi[datmi$Seascapes %in% combis, ]
  #datmi=datmi[-c(12),]
  #datmi=InsertRow(datmi, NewRow = row510, RowNum = 8)
  
  p=ggplot(data = datmi, aes(x=Seascapes, y=lag, fill = bcmi))+
    geom_tile(color = "white")+
    scale_fill_gradientn(colors=mypalette,limits=c(-0.1,1))+#(colours=pal,limits=c(-0.1,1)) + #rainbow(21)
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size=9))+
    theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1,size=9))+
    scale_y_continuous(breaks = round(seq(min(datmi$lag), max(datmi$lag), by = 1),1))

ggsave(p,file = paste(path,sprintf("Seascapes/Q1/%s/GAM/MI/NBTiled_bcmi_%s.png",species_info$file_name[i],species_info$file_name[i]),sep=""), height=2.5, width=4,dpi = 600)
}          




#Relationship between MI and Dissim?----
datpd=read.csv(paste(path,sprintf("Seascapes/Q1/%s/GAM/Diff_smooth/num_diff.csv",species_info$file_name[1]),sep=""))
datpd <- datpd[datpd$combi %in% combis, ]
datmi=read.csv(paste(path,sprintf("Seascapes/Q1/%s/GAM/MI/mutinf_lag.csv",species_info$file_name[1]),sep=""))
datmi$Seascapes=gsub("&","-",datmi$Seascapes) 
datmi$Seascapes=factor(datmi$Seascapes, levels = c("1 - 2", "1 - 3", "1 - 4", "1 - 5", "1 - 6", "1 - 7", "1 - 8", "1 - 9", "1 - 10",
                                                   "2 - 3", "2 - 4", "2 - 5", "2 - 6", "2 - 7", "2 - 8", "2 - 9", "2 - 10", "3 - 4",
                                                   "3 - 5", "3 - 6", "3 - 7", "3 - 8", "3 - 9", "3 - 10", "4 - 5", "4 - 6", "4 - 7",
                                                   "4 - 8", "4 - 9", "4 - 10", "5 - 6","5 - 7", "5 - 8", "5 - 9", "5 - 10", "6 - 7",
                                                   "6 - 8", "6 - 9", "6 - 10", "7 - 8", "7 - 9", "7 - 10", "8 - 9", "8 - 10", "9 - 10"))
datmi <- datmi[datmi$Seascapes %in% combis, ]
#datmi = subset(datmi,lag==0)


for (i in c(2,3,4,5,7,12,13,14)){
  datpdd=read.csv(paste(path,sprintf("Seascapes/Q1/%s/GAM/Diff_smooth/num_diff.csv",species_info$file_name[i]),sep=""))
  datpdd <- datpdd[datpdd$combi %in% combis, ]
  datpd=rbind(datpd,datpdd)
  
  datmii=read.csv(paste(path,sprintf("Seascapes/Q1/%s/GAM/MI/mutinf_lag.csv",species_info$file_name[i]),sep=""))
  datmii$Seascapes=gsub("&","-",datmii$Seascapes) 
  print(species_info$file_name[i])
  datmii$Seascapes=factor(datmii$Seascapes, levels = c("1 - 2", "1 - 3", "1 - 4", "1 - 5", "1 - 6", "1 - 7", "1 - 8", "1 - 9", "1 - 10",
                                                     "2 - 3", "2 - 4", "2 - 5", "2 - 6", "2 - 7", "2 - 8", "2 - 9", "2 - 10", "3 - 4",
                                                     "3 - 5", "3 - 6", "3 - 7", "3 - 8", "3 - 9", "3 - 10", "4 - 5", "4 - 6", "4 - 7",
                                                     "4 - 8", "4 - 9", "4 - 10", "5 - 6","5 - 7", "5 - 8", "5 - 9", "5 - 10", "6 - 7",
                                                     "6 - 8", "6 - 9", "6 - 10", "7 - 8", "7 - 9", "7 - 10", "8 - 9", "8 - 10", "9 - 10"))
  datmii <- datmii[datmii$Seascapes %in% combis, ]
  datmii = subset(datmii, lag==0)
  
  datmi=rbind(datmi,datmii)
}


data=cbind(datmi,datpd)
data$X=NULL
s=ggplot(data=data,aes(x=bcmi,y=period_diff))+geom_point(size=2,alpha=0.1)+
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,size=9))+
  theme(axis.text.y = element_text(vjust = 1, hjust = 1,size=9))+
  labs(x = "BCMI", y = 'Dissimilarity')
s
ggsave(s,file = paste(path,sprintf("images/bcmi_dissim.png"),sep=""), height=2.5, width=4,dpi = 600)

#VARRANK signal importance
c1=c(1,1,2,2,3,3,4,5,5,5,5,5,6,7,7,8,8,9)
c2=c(2,3,3,4,5,6,5,6,7,8,9,10,8,8,9,9,10,10)
c3=cbind(c1,c2)

#function to prep data
prep_data=function(dato){
  data=dato
  names(data)[names(data) == "id"] <- "Seascapenr"
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  ssc_data=completeFun(data,"Seascapenr")
  ssc_data=subset(ssc_data,Quarter==1)
  data_agg=ssc_data
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species+Ship+DayNight, data=ssc_data, sum)
  
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
  
  ###############################################################################
  data_agg_ansc=rbind(data_agg1,data_agg2,data_agg3,data_agg4,data_agg5,data_agg6,data_agg7,data_agg8,data_agg9,data_agg10)
  Mean_data=aggregate(biomass_kg~Year+Seascapenr+Species,data=data_agg_ansc,mean)
  return(Mean_data)
}


library(varrank)
data(PimaIndiansDiabetes, package = "mlbench")
varrank.PimaIndiansDiabetes <- varrank(data.df = PimaIndiansDiabetes, method = "estevez", variable.important = "diabetes", discretization.method = "sturges", algorithm = "forward", scheme="mid", verbose = FALSE)
plot(varrank.PimaIndiansDiabetes)

data=read.csv(paste(path,sprintf("Data/Filtered_%s.csv",species_info$file_name[1]),sep=""))
data=prep_data(data)

for (i in 2:16){
  dati=read.csv(paste(path,sprintf("Data/Filtered_%s.csv",species_info$file_name[i]),sep=""))
  dat=prep_data(dati)
  data=rbind(data,dat)
}

dat =data %>% spread(Species,biomass_kg)
write.csv(dat,paste(path,"Data/Mean_biomass_year_all_species.csv",sep=""))
