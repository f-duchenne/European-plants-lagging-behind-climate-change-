#ASSOCIER LES DONNES CLIMATIQUES AUX DONNEES D'OCCURENCE
library(data.table)
library(rgbif)
library(dplyr)
library(raster)
memory.limit(size = 1e9)
setwd(dir="C:/Users/Francois/Documents/land use change/Ecosystem_types_of_Europe_Terrestrial/Ecosystem types of Europe - version 3.1 Full map")
aqua=raster(paste0("aqua.tif"))
Grasslands=raster("Grasslands.tif")
Heathland=raster("Heathland.tif")
Woodland=raster("Woodland.tif")
Agricole=raster("Agricole.tif")
Urbain=raster("Urbain.tif")
Roche=raster("Roche.tif")
setwd(dir="C:/Users/Francois/Documents/land use change/biogeographic_regions")
Bio=raster("BiogeoRegions2016_raster.tif")

setwd(dir="C:/Users/Francois/Documents/gabi plantes/data")
res=fread("number_data_by_grid_cell_by_year_coord_round_to_10km.txt",sep="\t",head=T)
noms=c("ann_mean","diurn_range","max_temp","min_temp","prec","prec_max","prec_min")
all=as.data.frame(res %>% dplyr::group_by(year,Longitude2,Latitude2) %>% dplyr::summarise(nplant_tot=sum(n)))
all[,c(noms,"aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche","zone")]=NA


for(j in unique(all$year)){
ban=which(c(1901:2014)==j) # with thefirst and the last layer corresponding to theyears 1901 and 2014
setwd(dir="C:/Users/Francois/Documents/land use change/time_series_bioclim")
r=stack(paste0("bioclim_variables_year_",ban,"_10km_eumedclim.tif"))
r=projectRaster(r, projectExtent(r,projection(aqua)))
names(r)=noms
r[]<- r[]*0.1
res=subset(all,year==j)
step <- 10000
data_size <- length(res$Longitude2)

for(i in seq(from = 1, to = data_size, by = step)){
  
  i1 <- i
  i2 <- i + (step-1)
  if (i2 > data_size) {i2 <- data_size}
  
  print(c(i1,i2))
  pts <- res[i1:i2,c('Longitude2','Latitude2')]
  sf_pts <- sf::st_transform(sf::st_as_sf(pts,coords=c('Longitude2','Latitude2'),crs=projection(aqua)),raster::projection(r)) # be sure that species obs points and raster are  on the same projection
  res[i1:i2,noms]=raster::extract(r, as(sf_pts, 'Spatial'))
  res[i1:i2,"aqua"]=raster::extract(aqua, as(sf_pts, 'Spatial'))
  res[i1:i2,"Grasslands"]=raster::extract(Grasslands, as(sf_pts, 'Spatial'))
  res[i1:i2,"Heathland"]=raster::extract(Heathland, as(sf_pts, 'Spatial'))
  res[i1:i2,"Woodland"]=raster::extract(Woodland, as(sf_pts, 'Spatial'))
  res[i1:i2,"Agricole"]=raster::extract(Agricole, as(sf_pts, 'Spatial'))
  res[i1:i2,"Urbain"]=raster::extract(Urbain, as(sf_pts, 'Spatial'))
  res[i1:i2,"Roche"]=raster::extract(Roche, as(sf_pts, 'Spatial'))
  res[i1:i2,"zone"]=raster::extract(Bio, as(sf_pts, 'Spatial'))
  
}
if(j==1951){allf=res}else{allf=rbind(allf,res)}
}

dim(allf)
#Petit bout pour récupérer ce qui est un peu loin en mer mais qui match que les raster de température:
setwd(dir="C:/Users/Francois/Documents/land use change/biogeographic_regions")
Bio=raster("BiogeoRegions2016_raster.tif")
Bio<- aggregate(Bio,fact=5,fun=modal,expand=T,na.rm=T)
sf_pts <- sf::st_as_sf(allf,coords=c('Longitude2','Latitude2'),crs=projection(Bio))
allf$zone[is.na(allf$zone)]=raster::extract(Bio, as(sf_pts[is.na(allf$zone),], 'Spatial'))
Bio<- aggregate(Bio,fact=2,fun=modal,expand=T,na.rm=T)
allf$zone[is.na(allf$zone)]=raster::extract(Bio, as(sf_pts[is.na(allf$zone),], 'Spatial'))
points(Latitude2~Longitude2,data=subset(allf,is.na(zone)))
Bio<- aggregate(Bio,fact=2,fun=modal,expand=T,na.rm=T)
allf$zone[is.na(allf$zone)]=raster::extract(Bio, as(sf_pts[is.na(allf$zone),], 'Spatial'))

#rajouter le nom de la zone:
bioreg=data.frame(zone=c(1,2,3,4,5,6,7,8,9,10,11,12),region=c("Alpine","Anatolian","Arctic","Atlantic",
"BlackSea","Boreal","Continental","Macaronesia","Mediterranean","Outside","Pannonian","Steppic"))
bioreg=as.data.table(bioreg)
allf2=merge(allf,bioreg,by="zone",all.x=T,all.y=F)
#exporter le tout
setwd(dir="C:/Users/Francois/Documents/gabi plantes/data")
fwrite(allf2,paste0("climatic_indexes_by_grid_cell_by_year.txt"),sep="\t",row.names=F)

#############################################