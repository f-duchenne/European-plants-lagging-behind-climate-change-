#CACLCULER LES SCI TRENDS
library(data.table)
library(rgbif)
library(dplyr)
library(raster)
library(Hmisc)
library(car)
library(doBy)
library(ggcorrplot)
memory.limit(size = 1e9)
setwd(dir="C:/Users/Francois/Documents/gabi plantes/data")
all=fread("climatic_indexes_by_grid_cell_by_year.txt",sep="\t",head=T)
res=fread("number_data_by_grid_cell_by_year_coord_round_to_10km.txt",sep="\t",head=T)
res=merge(all,res,by=c("Longitude2","Latitude2","year"))
res=res[!is.na(res$speciesKey),]
res=res[res$year>1950,]
all=all[all$year>1950,]
res$ann_mean=as.numeric(res$ann_mean)
res$max_temp=as.numeric(res$max_temp)
res$min_temp=as.numeric(res$min_temp)
res$prec=as.numeric(res$prec)
res$prec_max=as.numeric(res$prec_max)
res$prec_min=as.numeric(res$prec_min)

wtd.mean2=function(x,y){
sum(x*y,na.rm=T)/sum(y,na.rm=T)}
names(res)[names(res)=="n"]="nb" #change this name, because it is also a function name in dplyr

#reviewer:
bidon=res %>% dplyr::group_by(year) %>% dplyr::summarise(urba=wtd.mean2(Urbain,nb)*100)
plot(urba~year,data=bidon,ylab="Weighted average of proportion of urban area in visited grid cells (in %)",pch=19)
model=lm(urba~year,data=bidon)
abline(model,lty="dashed")

bidon=res %>% group_by(speciesKey,year) %>% summarise(nbgrids=length(nb),lat=mean(Latitude2),lon=mean(Longitude2),mean_bydata=wtd.mean2(ann_mean,nb),
ann_mean=wtd.mean2(ann_mean,nb/nplant_tot),max_temp=wtd.mean2(max_temp,nb/nplant_tot),min_temp=wtd.mean2(min_temp,nb/nplant_tot),
prec=wtd.mean2(prec,nb/nplant_tot),prec_max=wtd.mean2(prec_max,nb/nplant_tot),prec_min=wtd.mean2(prec_min,nb/nplant_tot))
biche=res %>% group_by(Latitude2,Longitude2,speciesKey) %>% summarise()


liste=data.frame(speciesKey=unique(res$speciesKey))
noms=c("ann_mean","max_temp","min_temp","prec","prec_max","prec_min")
which(liste$speciesKey==3119191)
for(i in 1:nrow(liste)){
bidon2=as.data.frame(bidon[bidon$speciesKey==liste$speciesKey[i],])
bidon2=subset(bidon2,year<=1990)
liste$min_year[i]=min(bidon2$year)
liste$max_year[i]=max(bidon2$year)
liste$nb_year[i]=nrow(bidon2)
for(j in 1:length(noms)){
if(nrow(bidon2)>10){
bidon2$vec=bidon2[,paste(noms[j])]
model=lm(vec~year,weights=sqrt(nbgrids),data=bidon2)
ano=Anova(model)
liste[i,paste0(noms[j])]=model$coef[2]
liste[i,paste0(noms[j],"_sde")]=summary(model)$coef[2,2]
liste[i,paste0(noms[j],"_pval")]=ano[1,4]
}else{
liste[i,paste0(noms[j])]=NA
liste[i,paste0(noms[j],"_sde")]=NA
liste[i,paste0(noms[j],"_pval")]=NA
}
}
#print(boxplot(liste$pref_trend))
}
fwrite(liste,"climatic_niche_trends_before_1990.txt",sep="\t",row.names=F)