completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
#set path----
path=("/Users/markrademaker/Projects/Sea-scape-analysis/cpue_length_hour/")
#file with species info----
#species_info=read.table(paste(path,"Data/spec_info.txt",sep=""),sep="\t",header=TRUE)
#species_info$file_name =gsub(" ","_",species_info$Scientific.name)

###########################################
#### compute temporal correlation file #### ----
###########################################
#file with species info----
species_info=read.table(paste(path,"spec_info.txt",sep=""),sep="\t",header=TRUE)
species_info$file_name =gsub(" ","_",species_info$Scientific.name)

#selected species for temporal 
for (k in 1:nrow(species_info)){
  data=read.csv(paste(path,sprintf("Filtered_%s.csv",species_info$file_name[k]),sep=""))
  #names(data)[names(data) == "id"] <- "Seascapenr"
  print(species_info$file_name[k])
  # Unique combination between vessel, haulNo. and Year, i.e. each unique haul in dataset
  #print(data[3:7])
  data$ID = cumsum(!duplicated(data[3:7]))
 
  
  ssc_data = subset(data, Year > 1977 )#discard 1977
  ssc_data = completeFun(ssc_data, "Seascapenr")
  print(paste(nrow(ssc_data)," rows"))
  print(paste(length(unique(ssc_data$ID))," unique hauls"))
  #ssc_data=completeFun(data,"Seascapenr")
  #ssc_data=subset(ssc_data,Quarter==1)
  #data_agg=ssc_data
  data_agg=aggregate(Total_wgt ~ ID+Year+Seascapenr+Species+Ship+DayNight, data=ssc_data, sum)
  print(paste(nrow(data_agg), "unique hauls check"))
  
  non_zero = subset(data_agg, Total_wgt > 0)
  print(paste(nrow(non_zero),"non zero hauls"))
  
  # convert weights from grams to kg
  data_agg$biomass_kg=data_agg$Total_wgt/1000
  
  Mean_data=aggregate(biomass_kg~Year+Seascapenr,data=data_agg,mean)
  
  #1----
  data_agg1=subset(Mean_data,Seascapenr==1)
  print(nrow(data_agg1))
  #2-----
  data_agg2=subset(Mean_data,Seascapenr==2)
  print(nrow(data_agg2))
  
  #3----
  data_agg3=subset(Mean_data,Seascapenr==3)
  print(nrow(data_agg3))
  
  #4----
  data_agg4=subset(Mean_data,Seascapenr==4)
  print(nrow(data_agg4))
  
  #5----
  data_agg5=subset(Mean_data,Seascapenr==5)
  print(nrow(data_agg5))
  
  #6----
  data_agg6=subset(Mean_data,Seascapenr==6)
  print(nrow(data_agg6))
  
  #7----
  data_agg7=subset(Mean_data,Seascapenr==7)
  print(nrow(data_agg7))
  
  #8----
  data_agg8=subset(Mean_data,Seascapenr==8)
  print(nrow(data_agg8))
  
  #9----
  data_agg9=subset(Mean_data,Seascapenr==9)
  print(nrow(data_agg9))
  
  #10----
  data_agg10=subset(Mean_data,Seascapenr==10)
  print(nrow(data_agg10))
  
  data_meanyear=cbind(data_agg1$biomass_kg,data_agg2$biomass_kg,data_agg3$biomass_kg,data_agg4$biomass_kg,data_agg5$biomass_kg,data_agg6$biomass_kg,data_agg7$biomass_kg,data_agg8$biomass_kg,data_agg9$biomass_kg,data_agg10$biomass_kg)
  write.csv(data_meanyear,paste(path,sprintf("%s_PCMCI_prep.csv",species_info$file_name[k]),sep=""))
}

