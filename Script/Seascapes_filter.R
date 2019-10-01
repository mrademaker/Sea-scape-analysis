library(rgdal)
library(sf)
library(ggplot2)


#############################################
###### Determine position in sea scape#######
#############################################
path=("C:/Users/mrademaker/Documents/Research projects/STCNWS/")

##Trial Limanda - Limanda

#Open filtered dataset
data=read.csv(paste(path,"DATRAS/cpue_length_hour/filtered_Limanda_limanda.csv",sep=""))
coordinates(data)= ~ ShootLong + ShootLat

#Determine for each location in which seascape it is located----
ssc_map=readOGR(paste(path,"Seascapes/.",sep=""),layer="seascapes")
ssc_map_wgs84 <- spTransform(ssc_map, CRS("+proj=longlat +datum=WGS84"))


##for plotting----
proj4string(data) <- proj4string(ssc_map_wgs84)
shapefile_df <- fortify(ssc_map_wgs84)

map <- ggplot() +
  geom_path(data = shapefile_df, 
            aes(x = long, y = lat, group = group),
            color = 'black', size = 1)+
   ggtitle("Seascapes plot and NS-IBTS survey data (Limanda limanda)")
 print(map)
        
p2 = map+geom_point(data=as.data.frame(data),aes(x=ShootLong,y=ShootLat),color="black",fill="red",size=1,shape=21)
print(p2)


#Determine for each location in which seascape it is located----
ssc_map=readOGR(paste(path,"Seascapes/.",sep=""),layer="seascapes")
ssc_map_wgs84 <- spTransform(ssc_map, CRS("+proj=longlat +datum=WGS84"))

data=read.csv(paste(path,"DATRAS/cpue_length_hour/filtered_Limanda_limanda.csv",sep=""))
coordinates(data)= ~ ShootLong + ShootLat
proj4string(data) <- proj4string(ssc_map_wgs84)
data_in=as.data.frame(over(data, ssc_map_wgs84))
data_ssc=as.data.frame(cbind(as.data.frame(data),data_in))
write.csv(data_ssc,paste(path,"DATRAS/cpue_length_hour/filtered_Limanda_limanda.csv",sep=""))
#save as filtered


#Save file