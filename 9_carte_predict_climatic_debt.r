############################################################################
library(data.table)
library(rgbif)
library(raster)
library(Hmisc)
library(car)
library(doBy)
library(ggcorrplot)
library(lme4)
library(ape)
library(picante)
library(doBy)
library(ggtree)
library(phytools)
library(phylotools)
library(caper)
library(dplyr)
library(sf)
library(gridExtra)
library(optimx)

setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")
r2=raster("BiogeoRegions2016_raster.tif")
dat=as.data.frame(r2,xy=T,na.rm=T)
bioreg=data.frame(BiogeoRegions2016_raster=c(1,2,3,4,5,6,7,8,9,10,11,12),region=c("Alpine","Anatolian","Arctic","Atlantic",
"BlackSea","Boreal","Continental","Macaronesia","Mediterranean","Outside","Pannonian","Steppic"))
dat=merge(dat,bioreg,by="BiogeoRegions2016_raster")
sf_pts <- sf::st_as_sf(dat,coords=c('x','y'),crs=projection(r2)) 

noms=c("ann_mean","diurn_range","max_temp","min_temp","prec","prec_max","prec_min")


for(ban in 50:59){
r=stack(paste0("bioclim_variables_year_",ban,"_10km_eumedclim.tif"))
r=projectRaster(r, projectExtent(r,projection(r2)))
r[]<- r[]*0.1
bidon=raster::extract(r, as(sf_pts, 'Spatial'))
if(ban==50){bidonf=bidon}else{bidonf=bidonf+bidon}
}
noms=c("ann_mean","diurn_range","max_temp","min_temp","prec","prec_max","prec_min")
colnames(bidonf)=noms
bidonf=as.data.frame(bidonf)
dat$base_ann_mean=bidonf$ann_mean/10
dat$base_prec=bidonf$prec/10

for(ban in 105:114){
r=stack(paste0("bioclim_variables_year_",ban,"_10km_eumedclim.tif"))
r=projectRaster(r, projectExtent(r,projection(r2)))
r[]<- r[]*0.1
bidon=raster::extract(r, as(sf_pts, 'Spatial'))
if(ban==105){bidonf=bidon}else{bidonf=bidonf+bidon}
}
noms=c("ann_mean","diurn_range","max_temp","min_temp","prec","prec_max","prec_min")
colnames(bidonf)=noms
bidonf=as.data.frame(bidonf)
dat$temp_2005_2014=bidonf$ann_mean/10
dat$prec_2005_2014=bidonf$prec/10

dat$ann_mean=dat$temp_2005_2014-dat$base_ann_mean
dat$prec=dat$prec_2005_2014-dat$base_prec

aqua=raster(paste0("aqua.tif"))
Grasslands=raster("Grasslands.tif")
Heathland=raster("Heathland.tif")
Woodland=raster("Woodland.tif")
Agricole=raster("Agricole.tif")
Urbain=raster("Urbain.tif")
Roche=raster("Roche.tif")
dat[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]=NA
dat$aqua=raster::extract(aqua, as(sf_pts, 'Spatial'))
dat$Grasslands=raster::extract(Grasslands, as(sf_pts, 'Spatial'))
dat$Heathland=raster::extract(Heathland, as(sf_pts, 'Spatial'))
dat$Woodland=raster::extract(Woodland, as(sf_pts, 'Spatial'))
dat$Agricole=raster::extract(Agricole, as(sf_pts, 'Spatial'))
dat$Urbain=raster::extract(Urbain, as(sf_pts, 'Spatial'))
dat$Roche=raster::extract(Roche, as(sf_pts, 'Spatial'))

dat[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]=
dat[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]/
apply(dat[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")],1,sum)


#PHYLOGENIE
a=read.tree("DaPhnE_01.tre")
cod=read.table("code_canons.txt",sep="\t",header=T)
cod=subset(cod,rang=="SPECIES")
names(cod)[1]="noms_phylo"
cod=cod[!duplicated(cod$canon),]

#values
liste=fread("SCI_occupancy_traits_region.txt",header=T)
liste$life_span_cat[liste$life_span_cat==""]=NA
liste=subset(liste,nb_year>=20)
liste=subset(liste,type!="invasive")


liste=merge(liste,cod[,c("noms_phylo","canon")],by.x="species",by.y="canon",all.x=T,all.y=F)
liste$species=gsub(" ","_",liste$species,fixed=T)
liste$noms_phylo=gsub(" ","_",liste$noms_phylo,fixed=T)
length(liste$species[!(liste$species %in% a$tip.label)])
length(liste$espece[!(liste$espece %in% a$tip.label)])
liste$code=liste$species
length(liste$code[(liste$code %in% a$tip.label)])
liste$code[!(liste$code %in% a$tip.label)]=liste$noms_phylo[!(liste$code %in% a$tip.label)]
length(liste$code[(liste$code %in% a$tip.label)])
liste$code[!(liste$code %in% a$tip.label)]=liste$synonyms[!(liste$code %in% a$tip.label)]
length(liste$code[(liste$code %in% a$tip.label)])

liste2=liste[,c("trend","trend_err","code","nsite","nb_year","region","species","family","phylum",
"order","class","max_temp","ann_mean","prec","prec_min","base_ann_mean","base_prec",
"aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche",
"nitrogen_moy","dep_poll","moisture_moy","life_span_cat",
"Alpine","Anatolian","Arctic","Atlantic","BlackSea","Boreal","Continental","Mediterranean","Pannonian",
"Steppic")]
liste2=na.omit(liste2)
lili=c("Atlantic","Continental","Mediterranean","Alpine","Boreal")
dat$fit_lme=NA
dat$dep_poll=NA
dat$life_span_cat="annual"
dat$moisture_moy=NA
dat$nitrogen_moy=NA
dat$fit_pgls=NA
for(i in 1:length(lili)){
bidon=subset(liste2,paste0(region)==lili[i])
bidon[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]=
bidon[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]/
apply(bidon[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")],1,sum)
cor(bidon$prec,bidon$prec_min)

model2=lmer(trend~Urbain+Agricole+Heathland+aqua+Grasslands+Roche+
ann_mean+prec+moisture_moy+nitrogen_moy+dep_poll+
base_ann_mean+base_prec+life_span_cat+(1|family),data=bidon,weights=1/trend_err)


bidon=bidon[!is.na(bidon$code),]
bidon=bidon[!duplicated(bidon$code),]
p3d=comparative.data(a,as.data.frame(bidon),code,vcv=TRUE,na.omit=F)

model=pgls(trend~Urbain+Agricole+Heathland+aqua+Grasslands+Roche+
ann_mean+prec+moisture_moy+nitrogen_moy+dep_poll+
base_ann_mean+base_prec+life_span_cat,data=p3d,lambda='ML')

dat$moisture_moy[dat$region==lili[i]]=mean(bidon$moisture_moy)
dat$nitrogen_moy[dat$region==lili[i]]=mean(bidon$nitrogen_moy)
dat$dep_poll[dat$region==lili[i]]=mean(bidon$dep_poll)

dat$debt_temp[dat$region==lili[i]]=summary(model2)$coeff["ann_mean",1]
dat$debt_prec[dat$region==lili[i]]=summary(model2)$coeff["prec",1]
dat$debt_temp_pgls[dat$region==lili[i]]=summary(model)$coeff["ann_mean",1]
dat$debt_prec_pgls[dat$region==lili[i]]=summary(model)$coeff["prec",1]
dat$fit_lme[dat$region==lili[i]]=predict(model2,newdata=dat[dat$region==lili[i],],re.form=NA)
print(i)
}


dat$debt_lme=dat$ann_mean*dat$debt_temp+dat$prec*dat$debt_prec
dat$debt_pgls=dat$ann_mean*dat$debt_temp_pgls+dat$prec*dat$debt_prec_pgls

fwrite(dat,"climatic_debt_by_regions.txt")

############################################################################
library(data.table)
library(rgbif)
library(raster)
library(Hmisc)
library(car)
library(doBy)
library(ggcorrplot)
library(lme4)
library(ape)
library(picante)
library(doBy)
library(ggtree)
library(phytools)
library(phylotools)
library(caper)
library(dplyr)
library(sf)
library(gridExtra)
library(optimx)

setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")
dat=fread("climatic_debt_by_regions.txt")
dat=subset(dat,!is.na(debt_lme))
dat$debt=apply(dat[,c("debt_pgls","debt_lme")],1,mean)

r2=raster("BiogeoRegions2016_raster.tif")
regions2=as.data.frame(r2,xy=T,na.rm=T)
regions2=subset(regions2,BiogeoRegions2016_raster!=10)
regions2=st_as_sf(regions2,coords=c("x","y"),crs=projection(r2))
regions2=st_transform(regions2,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
regions2=cbind(regions2,st_coordinates(regions2))
regions2=regions2 %>% filter(X>=(-13) & X<=34 & Y>=34 & Y<=75)
regions2=st_transform(regions2,projection(r2))
regions2=cbind(regions2,st_coordinates(regions2))
regions2$X.1=plyr::round_any(regions2$X.1,50000)
regions2$Y.1=plyr::round_any(regions2$Y.1,50000)
conti2=subset(regions2,BiogeoRegions2016_raster==7)
setwd(dir="C:/Users/Francois/Documents/land use change/biogeographic_regions")
regions=st_read("BiogeoRegions2016.shp")
regions=st_transform(regions,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
regions=st_crop(regions,xmin=-13,ymin=34,xmax=34,ymax=75)
regions=st_transform(regions,projection(r2))
regions=subset(regions,code!="Outside")
conti=subset(regions,code=="Continental")




library(rgdal)
setwd(dir="C:/Users/Francois/Documents/papier 1 - données collection/data/shape")
shp <- st_read(".", "Europe_coastline")
setwd(dir="C:/Users/Francois/Documents/land use change/continents")
shp2 <- st_read("continents-of-the-world-merged.shp")
shp2=st_transform(shp2,projection(r2))
shp2=st_crop(shp2,extent(shp))

# fond=as.data.frame(r2,xy=T,na.rm=T)
# fond=subset(fond,BiogeoRegions2016_raster!=10)
# fond=st_as_sf(fond,coords=c("x","y"),crs=projection(r2))
# fond=st_transform(fond,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# fond=cbind(fond,st_coordinates(fond))
# fond=fond %>% filter(X>=(-13) & X<=34 & Y>=34 & Y<=75)
# fond=st_transform(fond,projection(r2))

sf_pts <- sf::st_as_sf(dat,coords=c('x','y'),crs=projection(r2)) 
sf_pts=st_transform(sf_pts,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
sf_pts=cbind(sf_pts,st_coordinates(sf_pts))
sf_pts=sf_pts %>% filter(X>=(-13) & X<=34 & Y>=34 & Y<=75)
setwd(dir="C:/Users/Francois/Documents/land use change/projection_2041-2060")
r1=brick("wc2.1_2.5m_bioc_CNRM-CM6-1_ssp370_2041-2060.tif")
sf_pts$temp_2041_2060=extract(r1[[1]], as(sf_pts, 'Spatial'))
sf_pts$prec_2041_2060=extract(r1[[12]], as(sf_pts, 'Spatial'))
sf_pts$debt_risk_lme=(sf_pts$temp_2041_2060-sf_pts$temp_2005_2014)*sf_pts$debt_temp+
(sf_pts$prec_2041_2060-sf_pts$prec_2005_2014)*sf_pts$debt_prec
sf_pts$debt_risk_pgls=(sf_pts$temp_2041_2060-sf_pts$temp_2005_2014)*sf_pts$debt_temp_pgls+
(sf_pts$prec_2041_2060-sf_pts$prec_2005_2014)*sf_pts$debt_prec_pgls
sf_pts$debt_risk=apply(st_drop_geometry(sf_pts[,c("debt_risk_pgls","debt_risk_lme")]),1,mean)


sf_pts=st_transform(sf_pts,projection(r2))

 library('scales')
mid=abs(min(sf_pts$debt,na.rm=T))/(max(sf_pts$debt,na.rm=T)-min(sf_pts$debt,na.rm=T))
pl1=ggplot()+geom_sf(data=shp2,color="black",fill="grey")+
geom_sf(data=regions,fill="black",col="black")+
geom_sf(data=sf_pts,aes(col=debt),size=0.3)+theme_bw()+
theme(axis.title=element_blank(),axis.text=element_blank(),
axis.ticks=element_blank(),plot.title=element_text(size=14,face="bold"),
panel.grid=element_line(color="darkgrey"),
axis.line=element_blank(),panel.border=element_blank(),panel.background=element_rect(fill="white"))+
scale_colour_gradientn(colors=c("firebrick","pink","white","deepskyblue","dodgerblue4"),
values=c(0,mid-0.1,mid,mid+0.1,1))+
labs(color="")+coord_fixed(ratio = 1)+ggtitle("a")+
coord_map()+
coord_sf(crs = projection(r2))+geom_sf(data=shp2,color="black",fill=NA)+ylim(extent(sf_pts)[3:4])+
xlim(extent(sf_pts)[1:2])


sf_pts$temp_contrib=(c(sf_pts$ann_mean*sf_pts$debt_temp_pgls)+c(sf_pts$ann_mean*sf_pts$debt_temp))/2
sf_pts$prec_contrib=(c(sf_pts$prec*sf_pts$debt_prec_pgls)+c(sf_pts$prec*sf_pts$debt_prec))/2
sf_pts$contrib=sf_pts$temp_contrib/(abs(sf_pts$prec_contrib)+abs(sf_pts$temp_contrib))
mid=abs(min(sf_pts$contrib,na.rm=T))/(max(sf_pts$contrib,na.rm=T)-min(sf_pts$contrib,na.rm=T))
pl2=ggplot()+geom_sf(data=shp2,color="black",fill="grey")+
geom_sf(data=regions,fill="black",col="black")+
geom_sf(data=sf_pts,aes(col=contrib),size=0.3)+theme_bw()+
theme(axis.title=element_blank(),axis.text=element_blank(),
axis.ticks=element_blank(),plot.title=element_text(size=14,face="bold"),
panel.grid=element_line(color="darkgrey"),
axis.line=element_blank(),panel.border=element_blank(),panel.background=element_rect(fill="white"))+
scale_colour_gradientn(colors=c("firebrick","white","dodgerblue4"),values=c(0,mid,1),breaks=c(-0.5,0,0.5),labels = percent(c(-0.5,0,0.5)))+
labs(color="")+coord_fixed(ratio = 1)+ggtitle("b")+
coord_map()+
coord_sf(crs = projection(r2))+geom_sf(data=shp2,color="black",fill=NA)+ylim(extent(sf_pts)[3:4])+
xlim(extent(sf_pts)[1:2])


sf_pts$contrib=sf_pts$prec_contrib/(abs(sf_pts$prec_contrib)+abs(sf_pts$temp_contrib))
mid=abs(min(sf_pts$contrib,na.rm=T))/(max(sf_pts$contrib,na.rm=T)-min(sf_pts$contrib,na.rm=T))
pl3=ggplot()+geom_sf(data=shp2,color="black",fill="grey")+
geom_sf(data=regions,fill="black",col="black")+
geom_sf(data=sf_pts,aes(col=contrib),size=0.3)+theme_bw()+
theme(axis.title=element_blank(),axis.text=element_blank(),
axis.ticks=element_blank(),plot.title=element_text(size=14,face="bold"),
panel.grid=element_line(color="darkgrey"),
axis.line=element_blank(),panel.border=element_blank(),panel.background=element_rect(fill="white"))+
scale_colour_gradientn(colors=c("firebrick","white","dodgerblue4"),values=c(0,mid,1),breaks=c(-0.5,0,0.5),labels = percent(c(-0.5,0,0.5)))+
labs(color="")+coord_fixed(ratio = 1)+ggtitle("c")+
coord_map()+
coord_sf(crs = projection(r2))+geom_sf(data=shp2,color="black",fill=NA)+ylim(extent(sf_pts)[3:4])+
xlim(extent(sf_pts)[1:2])


# grid.arrange(pl1,pl2,pl3,ncol=3)

setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")
png("fig_S9.png",width=1500,height=600,res=120)
grid.arrange(pl1,pl2,pl3,ncol=3)
dev.off();

pdf("fig_S9.pdf",width=15,height=6)
grid.arrange(pl1,pl2,pl3,ncol=3)
dev.off();

#ENCONSIDERANT QUE LES EFFETS SIGNIFICATIFS:
 library('scales')
sf_pts$temp_contrib=(c(sf_pts$ann_mean*sf_pts$debt_temp_pgls)+c(sf_pts$ann_mean*sf_pts$debt_temp))/2
sf_pts$prec_contrib=(c(sf_pts$prec*sf_pts$debt_prec_pgls)+c(sf_pts$prec*sf_pts$debt_prec))/2
sf_pts$temp_contrib[sf_pts$region %in% c("Atlantic","Continental")]=NA
sf_pts$prec_contrib[sf_pts$region %in% c("Boreal","Continental")]=NA
sf_pts$debt=apply(st_drop_geometry(sf_pts[,c("temp_contrib","prec_contrib")]),1,sum,na.rm=T)
sf_pts$debt[sf_pts$region %in% c("Continental")]=NA


na.value.forplot <- 'antiquewhite2'

mid=abs(min(sf_pts$debt,na.rm=T))/(max(sf_pts$debt,na.rm=T)-min(sf_pts$debt,na.rm=T))
pl1=ggplot()+geom_sf(data=shp2,color="black",fill="grey")+
geom_sf(data=regions,fill="black",col="black")+
geom_sf(data=sf_pts,aes(col=debt,fill=na.value.forplot),size=0.3)+theme_bw()+
theme(axis.title=element_blank(),axis.text=element_blank(),
axis.ticks=element_blank(),plot.title=element_text(size=14,face="bold"),
panel.grid=element_line(color="darkgrey"),
axis.line=element_blank(),panel.border=element_blank(),panel.background=element_rect(fill="white"))+
scale_colour_gradientn(colors=c("firebrick","pink","white","slategray2","dodgerblue4"),
values=c(0,mid-0.1,mid,mid+0.1,1),na.value=na.value.forplot)+
labs(color="")+coord_fixed(ratio = 1)+ggtitle("a")+
coord_map()+
coord_sf(crs = projection(r2))+geom_sf(data=shp2,color="black",fill=NA)+ylim(extent(sf_pts)[3:4])+
xlim(extent(sf_pts)[1:2])+scale_fill_discrete(labels="non-significant",
guide=guide_legend("",labels="Unsignificant",override.aes=list(colour="black",
fill='antiquewhite2',shape=22,size=5)))



sf_pts$contrib=sf_pts$temp_contrib/apply(abs(st_drop_geometry(sf_pts[,c("temp_contrib","prec_contrib")])),1,sum,na.rm=T)
mid=abs(min(sf_pts$contrib,na.rm=T))/(max(sf_pts$contrib,na.rm=T)-min(sf_pts$contrib,na.rm=T))
pl2=ggplot()+geom_sf(data=shp2,color="black",fill="grey")+
geom_sf(data=regions,fill="black",col="black")+
geom_sf(data=sf_pts,aes(col=contrib,fill=na.value.forplot),size=0.3)+theme_bw()+
theme(axis.title=element_blank(),axis.text=element_blank(),
axis.ticks=element_blank(),plot.title=element_text(size=14,face="bold"),
panel.grid=element_line(color="darkgrey"),
axis.line=element_blank(),panel.border=element_blank(),panel.background=element_rect(fill="white"))+
scale_colour_gradientn(colors=c("firebrick","white","dodgerblue4"),values=c(0,mid,1),
breaks=c(-0.5,0,0.5),labels = percent(c(-0.5,0,0.5)),,na.value=na.value.forplot)+
labs(color="")+coord_fixed(ratio = 1)+ggtitle("b")+
coord_map()+
coord_sf(crs = projection(r2))+geom_sf(data=shp2,color="black",fill=NA)+ylim(extent(sf_pts)[3:4])+
xlim(extent(sf_pts)[1:2])+scale_fill_discrete(labels="non-significant",
guide=guide_legend("",labels="Unsignificant",override.aes=list(colour="black",
fill='antiquewhite2',shape=22,size=5)))


sf_pts$contrib=sf_pts$prec_contrib/apply(abs(st_drop_geometry(sf_pts[,c("temp_contrib","prec_contrib")])),1,sum,na.rm=T)
mid=abs(min(sf_pts$contrib,na.rm=T))/(max(sf_pts$contrib,na.rm=T)-min(sf_pts$contrib,na.rm=T))
pl3=ggplot()+geom_sf(data=shp2,color="black",fill="grey")+
geom_sf(data=regions,fill="black",col="black")+
geom_sf(data=sf_pts,aes(col=contrib,fill=na.value.forplot),size=0.3)+theme_bw()+
theme(axis.title=element_blank(),axis.text=element_blank(),
axis.ticks=element_blank(),plot.title=element_text(size=14,face="bold"),
panel.grid=element_line(color="darkgrey"),
axis.line=element_blank(),panel.border=element_blank(),panel.background=element_rect(fill="white"))+
scale_colour_gradientn(colors=c("firebrick","white","dodgerblue4"),values=c(0,mid,1),
breaks=c(-0.5,0,0.5),labels = percent(c(-0.5,0,0.5)),,na.value=na.value.forplot)+
labs(color="")+coord_fixed(ratio = 1)+ggtitle("c")+
coord_map()+
coord_sf(crs = projection(r2))+geom_sf(data=shp2,color="black",fill=NA)+ylim(extent(sf_pts)[3:4])+
xlim(extent(sf_pts)[1:2])+scale_fill_discrete(labels="non-significant",
guide=guide_legend("",labels="Unsignificant",override.aes=list(colour="black",
fill='antiquewhite2',shape=22,size=5)))

grid.arrange(pl1,pl2,pl3,ncol=3)

setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")
png("fig3.png",width=1500,height=600,res=120)
grid.arrange(pl1,pl2,pl3,ncol=3)
dev.off();


pdf("figure3.pdf",width=15,height=6)
grid.arrange(pl1,pl2,pl3,ncol=3)
dev.off();

