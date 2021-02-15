#ASSOCIER LES DONNES CLIMATIQUES AUX DONNEES D'OCCURENCE
library(data.table)
library(rgbif)
library(dplyr)
library(raster)
library(Hmisc)

setwd(dir="C:/Users/Francois/Documents/gabi plantes/data")
all2=fread("climatic_indexes_by_grid_cell_by_year.txt",sep="\t",head=T)
res=fread("number_data_by_grid_cell_by_year_coord_round_to_10km.txt",sep="\t",head=T)
res=merge(all2,res,by=c("Longitude2","Latitude2","year"))
res=res[!is.na(res$speciesKey),]
res=res[res$year>1950,]
res$ann_mean=as.numeric(res$ann_mean)
res$max_temp=as.numeric(res$max_temp)
res$min_temp=as.numeric(res$min_temp)
res$prec=as.numeric(res$prec)
res$prec_max=as.numeric(res$prec_max)
res$prec_min=as.numeric(res$prec_min)
res$prec_min=as.numeric(res$prec_min)

wtd.mean2=function(x,y){
sum(x*y,na.rm=T)/sum(y,na.rm=T)}
res$nb=res$n

bidon=subset(res,!is.na(ann_mean)) %>% group_by(speciesKey,year) %>% 
summarise(nbgrids=length(nb),lat=mean(Latitude2),lon=mean(Longitude2),
mean_bydata=wtd.mean2(ann_mean,nb),ann_mean=wtd.mean2(ann_mean,nb/nplant_tot),
max_temp=wtd.mean2(max_temp,nb/nplant_tot),min_temp=wtd.mean2(min_temp,n/nplant_tot),
prec=wtd.mean2(prec,nb/nplant_tot),prec_max=wtd.mean2(prec_max,nb/nplant_tot),
prec_min=wtd.mean2(prec_min,nb/nplant_tot))

bidon1=subset(bidon,year<=1980 & !is.na(ann_mean)) %>% group_by(speciesKey) %>% summarise(base_ann_mean=mean(ann_mean,na.rm=T),
base_max_temp=mean(max_temp,na.rm=T),base_min_temp=mean(min_temp,na.rm=T),base_prec=mean(prec,na.rm=T),
base_prec_max=mean(prec_max,na.rm=T),base_prec_min=mean(prec_min,na.rm=T))

bidon2=subset(res,!is.na(Roche)) %>% group_by(speciesKey) %>% summarise(Urbain=wtd.mean2(Urbain,nb/nplant_tot),
Agricole=wtd.mean2(Agricole,nb/nplant_tot),aqua=wtd.mean2(aqua,nb/nplant_tot),
Grasslands=wtd.mean2(Grasslands,nb/nplant_tot),Woodland=wtd.mean2(Woodland,nb/nplant_tot),
Heathland=wtd.mean2(Heathland,nb/nplant_tot),Roche=wtd.mean2(Roche,nb/nplant_tot))
summary(bidon2)
bidon2[bidon2$sepciesKey %in% c(2704120,2965316,4121120),]
bidon2$tot=apply(bidon2[,c(2:8)],1,sum)

bidon3=subset(res,region!="") %>% group_by(speciesKey,region) %>% summarise(ratio=sum(nb/nplant_tot))
bidon3=bidon3 %>% group_by(speciesKey) %>% mutate(ratio_tot=sum(ratio))
bidon3$ratio=bidon3$ratio/bidon3$ratio_tot
bidon3=dcast(bidon3,speciesKey~region,value.var="ratio",fill=0)

traits=merge(bidon1,bidon2,by=c("speciesKey"))
traits=merge(traits,bidon3,by=c("speciesKey"))

setwd(dir="C:/Users/Francois/Documents/gabi plantes/data")
tab=fread("0053076-200221144449610.csv")
tab=subset(tab,taxonRank=="SPECIES" & taxonomicStatus=="ACCEPTED" & !is.na(speciesKey)) #CHARGER UNE LISTE D'ESPECE PROPRE AVEC NOMBRE DE DONNEES DANS ZONE D'INTERET (pour final check)
tab=subset(tab,phylum!="" & phylum!="Bryophyta" & phylum!="Chlorophyta" & phylum!="Rhodophyta" &
phylum!="Charophyta" & phylum!="Anthocerotophyta")
tab=tab[!duplicated(tab$speciesKey),] #remove duplicated GBIF ID (species correponding to the same one in the GBIF) to avoid double extraction
tab=subset(tab,numberOfOccurrences>=500)
tab=tab[,c("speciesKey","species","genus","family","order","class","phylum")]
traits=merge(traits,tab,by="speciesKey")

trait_tr8=fread("traits_climatic_debt.txt",sep="\t")
traits=merge(traits,trait_tr8,by="species")
traits=traits[,-which(names(traits) %in% c("li_span","life_span","poll_vect","ell_moist_uk","ell_N",
"poll_vect_B","ell_U_it","ell_N_it","poll_vect_fr","salt","nitrogen","moisture","nitrogen2",
"moisture2","nitrogen3","moisture3","nitrogen4","moisture4","nitrogen5","moisture5","annual","pluri",
"insect","total"))]
fwrite(traits,"traits_final_new.txt",sep="\t",row.names=F)












