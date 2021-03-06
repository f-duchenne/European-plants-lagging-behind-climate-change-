library(data.table)
library(rgbif)
library(dplyr)
library(raster)
library(rgdal)
library(sp)
memory.limit(size = 1e9)
#load the reference grid, from the temperature raster:
setwd(dir="C:/Users/Francois/Documents/land use change/time_series_bioclim")
r=stack("bioclim_variables_year_56_10km_eumedclim.tif")
crsnew=projection(r)
setwd(dir="C:/Users/Francois/Documents/land use change/biogeographic_regions")
r2=raster("BiogeoRegions2016_raster.tif") #layer that allow to consider only the European region
r=projectRaster(r, projectExtent(r,projection(r2))) #put the r raster on the right projection

#bioreg coding translation table:
bioreg=data.frame(zone=c(1,2,3,4,5,6,7,8,9,10,11,12),region=c("Alpine","Anatolian","Arctic","Atlantic",
"BlackSea","Boreal","Continental","Macaronesia","Mediterranean","Outside","Pannonian","Steppic "))

bioreg=as.data.table(bioreg)
library(sf)
library(data.table)
setwd(dir="D:/")
for(i in 1:10){
tab=fread(paste0("extraction_part_",i,".txt"),header=T) #load ith part of the total database 
tab=tab %>% filter(year>=1951 & year<=2014) #filter data
tab=tab %>% filter(decimalLongitude>=(-13)) #filter data
#make a spatial dataframe:
coord=unique(tab[,c("decimalLongitude","decimalLatitude")])
bidon=sf::st_transform(sf::st_as_sf(coord,coords = c('decimalLongitude','decimalLatitude'),crs=crsnew),
crs=projection(r2))
#extract the biogeographic region:
coord$zone=raster::extract(r2, as(bidon, 'Spatial'))
obj=raster::extract(r, as(bidon, 'Spatial'),cellnumbers=T) # set the occurence dataset at the same resolution than the raster of climatic data
coord=cbind(coord,coordinates(r)[obj[,1],])
coord=coord[!is.na(obj[,2]),]
names(coord)[4:5]=c("Longitude2","Latitude2")
coord=merge(coord,bioreg,by="zone")
coord=coord %>% filter(region!="Outside" & !is.na(region)) #remove records from North Africa
a=nrow(tab)

tab=merge(tab,coord,by=c("decimalLongitude","decimalLatitude")) # merge new coordinates with the database
tab3=tab %>% dplyr::group_by(speciesKey,year,Longitude2,Latitude2) %>% dplyr::count() #count the number of records by grid celle by species by year
if(i==1){res=tab3}else{res=rbind(res,tab3)}
print(paste(i,a-nrow(tab),sep=" - "))
}
# Export the total database:
setwd(dir="C:/Users/Francois/Documents/gabi plantes/data")
fwrite(res,paste0("number_data_by_grid_cell_by_year_coord_round_to_10km.txt"),sep="\t",row.names=F)

####################################################