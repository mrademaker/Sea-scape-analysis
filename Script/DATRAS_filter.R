#load packages----
library(dplyr)
library(tidyr)
library(rgdal)



#set path----
path=("/Users/markrademaker/Projects/Sea-scape-analysis/cpue_length_hour/")

#file with species info----
species_info=read.table(paste(path,"spec_info.txt",sep=""),sep="\t",header=TRUE)
species_info$file_name =gsub(" ","_",species_info$Scientific.name)

#open seascapes map and convert to wgs84 long,lat----
ssc_map=readOGR("/Users/markrademaker/Projects/Sea-scape-analysis/cpue_length_hour",layer="seascape_outline")
ssc_map_wgs84 <- spTransform(ssc_map, CRS("+proj=longlat +datum=WGS84")) 

#filter and adjust data_sets----
for (i in 1:nrow(species_info)){
  file_name=paste(path,species_info$file_name[i],sep="")
  print(file_name)
  #Read in data
  data=read.csv(paste(file_name,".csv",sep=""))
  print(paste(nrow(data),"number of rows raw datafile"))
  
  #Subset with same fishing gear
  data=subset(data, Gear == 'GOV')
  data=subset(data, Quarter == 1)
  
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
  #file_name=sprintf("Users/markrademaker/Projects/Sea-scape-analysis/cpue_length_hour/Filtered_%s.csv",species_info$file_name[i])
  print(paste(nrow(data_ssc)," number of rows filtered dataset"))
  write.csv(data_ssc,paste(path,sprintf("Filtered_%s.csv",species_info$file_name[i]),sep=""))  #
}

