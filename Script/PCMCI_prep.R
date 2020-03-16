
#set path----
path=("C:/Users/mrademaker/Documents/Research projects/STCNWS/DATRAS/cpue_length_hour/")
#file with species info----
species_info=read.table(paste(path,"Data/spec_info.txt",sep=""),sep="\t",header=TRUE)
species_info$file_name =gsub(" ","_",species_info$Scientific.name)

###########################################
#### compute temporal correlation file #### ----
###########################################
#file with species info----
species_info=read.table(paste(path,"Data/spec_info.txt",sep=""),sep="\t",header=TRUE)
species_info$file_name =gsub(" ","_",species_info$Scientific.name)

#selected species for temporal 
for (k in 1:nrow(species_info)){
  data=read.csv(paste(path,sprintf("Data/Filtered_%s.csv",species_info$file_name[k]),sep=""))
  names(data)[names(data) == "id"] <- "Seascapenr"
  print(species_info$file_name[k])
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  data$ID = cumsum(!duplicated(data[3:7]))
  
  ssc_data=completeFun(data,"Seascapenr")
  ssc_data=subset(ssc_data,Quarter==1)
  data_agg=ssc_data
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species+Ship+DayNight, data=ssc_data, sum)
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000
  
  
  Mean_data=aggregate(biomass_kg~Year+Seascapenr,data=data_agg,mean)
  
  #1----
  data_agg1=subset(Mean_data,Seascapenr==1)
  
  #2-----
  data_agg2=subset(Mean_data,Seascapenr==2)
 
  #3----
  data_agg3=subset(Mean_data,Seascapenr==3)
  
  #4----
  data_agg4=subset(Mean_data,Seascapenr==4)
  
  #5----
  data_agg5=subset(Mean_data,Seascapenr==5)
  
  #6----
  data_agg6=subset(Mean_data,Seascapenr==6)
  
  #7----
  data_agg7=subset(Mean_data,Seascapenr==7)
  
  #8----
  data_agg8=subset(Mean_data,Seascapenr==8)
  
  #9----
  data_agg9=subset(Mean_data,Seascapenr==9)
  
  #10----
  data_agg10=subset(Mean_data,Seascapenr==10)
  
  data_meanyear=cbind(data_agg1$biomass_kg,data_agg2$biomass_kg,data_agg3$biomass_kg,data_agg4$biomass_kg,data_agg5$biomass_kg,data_agg6$biomass_kg,data_agg7$biomass_kg,data_agg8$biomass_kg,data_agg9$biomass_kg,data_agg10$biomass_kg)
  write.csv(data_meanyear,paste(path,sprintf("Data/PCMCI_prep/%s_PCMCI_prep.csv",species_info$file_name[k]),sep=""))
}
