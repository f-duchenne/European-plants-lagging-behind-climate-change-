library(data.table)
library(dplyr)

#### ALL THE STUDIED PERIOD:

##EUROPEAN SCALE:
#occupancy trends
setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")
trends=read.table("occupancy_trends_europe.txt",sep="\t",header=T)
trends=subset(trends,param=="as.numeric(year)")
names(trends)[1:2]=c("trend","trend_err")
#Sci trends:
liste=fread("climatic_niche_trends.txt",sep="\t",header=T)
#merge both:
liste=merge(liste,trends,by.x=c("speciesKey"),by.y=c("Esp"),all=T)
#add species traits
trait=read.table("traits_final_new.txt",sep="\t",header=T)
trait[trait==""]=NA
liste=merge(liste,trait,by=c("speciesKey"),all.x=T)

#Select species:
#remove crop species
typeli=fread("liste_invasives_cultivees.txt",sep="\t")
liste=merge(liste,typeli[,c("ID","type")],by.x="speciesKey",by.y="ID",all.x=T,all.y=F)
liste$type[is.na(liste$type)]="wild"
liste=subset(liste,type!="cultivée" & !is.na(Urbain))
#remove crops define at genus level
typeli=fread("potentielles_invasives_cultivees_dapres_leurgenre_GM.txt",sep="\t")
#keep some species
typeli=subset(typeli,suggestion_enlever=="x")
liste=liste[!(liste$speciesKey %in% typeli$speciesKey),]
#remove ornemental species
typeli=fread("liste_filtrage_urbain_agricole_EP.txt",sep="\t")
#enlever de cette liste ce que l'on veut garder:
typeli=subset(typeli,a_enlever=="x")
liste=liste[!(liste$speciesKey %in% typeli$speciesKey),]

fwrite(liste,"SCI_occupancy_traits_Europe.txt",sep="\t")

##REGIONAL SCALE:
#occupancy trends
setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")
trends=read.table("occupancy_trends_region.txt",sep="\t",header=T)
trends=subset(trends,param=="as.numeric(year)")
names(trends)[1:2]=c("trend","trend_err")
#Sci trends:
liste=fread("climatic_niche_trends_by_region.txt",sep="\t",header=T)
#merge both:
liste=merge(liste,trends,by.x=c("speciesKey","region"),by.y=c("Esp","region"),all=T)
#add species traits
trait=read.table("traits_final_new.txt",sep="\t",header=T)
trait[trait==""]=NA
liste=merge(liste,trait,by=c("speciesKey"),all.x=T)

#Select species:
#remove crop species
typeli=fread("liste_invasives_cultivees.txt",sep="\t")
liste=merge(liste,typeli[,c("ID","type")],by.x="speciesKey",by.y="ID",all.x=T,all.y=F)
liste$type[is.na(liste$type)]="wild"
liste=subset(liste,type!="cultivée" & !is.na(Urbain))
#remove crops define at genus level
typeli=fread("potentielles_invasives_cultivees_dapres_leurgenre_GM.txt",sep="\t")
#keep some species
typeli=subset(typeli,suggestion_enlever=="x")
liste=liste[!(liste$speciesKey %in% typeli$speciesKey),]
#remove ornemental species
typeli=fread("liste_filtrage_urbain_agricole_EP.txt",sep="\t")
#enlever de cette liste ce que l'on veut garder:
typeli=subset(typeli,a_enlever=="x")
liste=liste[!(liste$speciesKey %in% typeli$speciesKey),]

fwrite(liste,"SCI_occupancy_traits_region.txt",sep="\t")


#### BEFORE 1990:

##EUROPEAN SCALE:
#occupancy trends
setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")
trends=read.table("occupancy_trends_europe_avant_1990.txt",sep="\t",header=T)
trends=subset(trends,param=="as.numeric(year)")
names(trends)[1:2]=c("trend","trend_err")
#Sci trends:
liste=fread("climatic_niche_trends_before_1990.txt",sep="\t",header=T)
#merge both:
liste=merge(liste,trends,by.x=c("speciesKey"),by.y=c("Esp"),all=T)
#add species traits
trait=read.table("traits_final_new.txt",sep="\t",header=T)
trait[trait==""]=NA
liste=merge(liste,trait,by=c("speciesKey"),all.x=T)

#Select species:
#remove crop species
typeli=fread("liste_invasives_cultivees.txt",sep="\t")
liste=merge(liste,typeli[,c("ID","type")],by.x="speciesKey",by.y="ID",all.x=T,all.y=F)
liste$type[is.na(liste$type)]="wild"
liste=subset(liste,type!="cultivée" & !is.na(Urbain))
#remove crops define at genus level
typeli=fread("potentielles_invasives_cultivees_dapres_leurgenre_GM.txt",sep="\t")
#keep some species
typeli=subset(typeli,suggestion_enlever=="x")
liste=liste[!(liste$speciesKey %in% typeli$speciesKey),]
#remove ornemental species
typeli=fread("liste_filtrage_urbain_agricole_EP.txt",sep="\t")
#enlever de cette liste ce que l'on veut garder:
typeli=subset(typeli,a_enlever=="x")
liste=liste[!(liste$speciesKey %in% typeli$speciesKey),]

fwrite(liste,"SCI_occupancy_traits_Europe_before_1990.txt",sep="\t")

##REGIONAL SCALE:
#occupancy trends
setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")
trends=read.table("occupancy_trends_region_avant_1990.txt",sep="\t",header=T)
trends=subset(trends,param=="as.numeric(year)")
names(trends)[1:2]=c("trend","trend_err")
#Sci trends:
liste=fread("climatic_niche_trends_by_region_before_1990.txt",sep="\t",header=T)
#merge both:
liste=merge(liste,trends,by.x=c("speciesKey","region"),by.y=c("Esp","region"),all=T)
#add species traits
trait=read.table("traits_final_new.txt",sep="\t",header=T)
trait[trait==""]=NA
liste=merge(liste,trait,by=c("speciesKey"),all.x=T)

#Select species:
#remove crop species
typeli=fread("liste_invasives_cultivees.txt",sep="\t")
liste=merge(liste,typeli[,c("ID","type")],by.x="speciesKey",by.y="ID",all.x=T,all.y=F)
liste$type[is.na(liste$type)]="wild"
liste=subset(liste,type!="cultivée" & !is.na(Urbain))
#remove crops define at genus level
typeli=fread("potentielles_invasives_cultivees_dapres_leurgenre_GM.txt",sep="\t")
#keep some species
typeli=subset(typeli,suggestion_enlever=="x")
liste=liste[!(liste$speciesKey %in% typeli$speciesKey),]
#remove ornemental species
typeli=fread("liste_filtrage_urbain_agricole_EP.txt",sep="\t")
#enlever de cette liste ce que l'on veut garder:
typeli=subset(typeli,a_enlever=="x")
liste=liste[!(liste$speciesKey %in% typeli$speciesKey),]


fwrite(liste,"SCI_occupancy_traits_region_before_1990.txt",sep="\t")


