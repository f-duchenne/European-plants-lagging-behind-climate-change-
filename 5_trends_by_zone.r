#MAKE OCCUPANCY MATRICES:
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
dim(res)
res=res[!is.na(res$speciesKey),]
dim(res)
res=res %>% group_by(Latitude2,Longitude2) %>% mutate(nyr=n_distinct(year))
res=res[res$nyr>=2,]
# bidon=res %>%  dplyr::select(Latitude2,Longitude2) %>%  distinct
# setwd(dir="C:/Users/Francois/Documents/land use change/elevation")
# elev=raster("elevation1x1_new.tif")
# sf_pts <- sf::st_transform(sf::st_as_sf(bidon, coords = c('Longitude2','Latitude2'), crs=4326), raster::projection(elev)) # match the same projection to species obs points and raster
# bidon$elev=raster::extract(elev, as(sf_pts, 'Spatial'))
# res=merge(as.data.table(res),as.data.table(bidon),by=c('Longitude2','Latitude2'))
res$site=paste(res$Latitude2,res$Longitude2,sep="-")
liste=unique(res$speciesKey)
stepsize=500
for(i in seq(1,length(liste),by=stepsize)){
i2=(i+stepsize)-1
res2=res
res2$speciesKey[!(res2$speciesKey %in% liste[i:i2])]="other"
res2=as.data.table(res2)
objformat=dcast(res2,site+year+region~speciesKey,
fun.aggregate=length,value.var="n")
objformat$tot=apply(objformat[,4:ncol(objformat)],1,sum)
setwd(dir="C:/Users/Francois/Documents/gabi plantes/data")
fwrite(objformat,paste0("objformat_",i,"_",i2,".txt"),sep="\t",row.names=F)
}

#################################################

.libPaths(c("/home/duchenne/R/x86_64-pc-linux-gnu-library/3.4/",.libPaths()))
library(dplyr)
library(data.table)
library(MCMCglmm)
library(doParallel)
library(foreach)
#library(doMPI)
library(igraph)
#cl<-startMPIcluster()
#registerDoMPI(cl)
cl <- makeCluster(96)
registerDoParallel(cl)


liste=fread("species_list_for_trends.txt",sep="\t")
liste=liste$speciesKey


deb=1
end=length(liste)
resultMCMC=foreach(e=deb:end,.combine=rbind)%dopar%{
library(data.table)
library(MCMCglmm)
library(glmmTMB)
library(dplyr)
library(mgcv)
if(e<=500){res=fread(paste0("objformat_1_500.txt"),sep="\t",header=T)}
if(e>500 & e<=1000){res=fread(paste0("objformat_501_1000.txt"),sep="\t",header=T)}
if(e>1000 & e<=1500){res=fread(paste0("objformat_1001_1500.txt"),sep="\t",header=T)}
if(e>1500 & e<=2000){res=fread(paste0("objformat_1501_2000.txt"),sep="\t",header=T)}
if(e>2000 & e<=2500){res=fread(paste0("objformat_2001_2500.txt"),sep="\t",header=T)}
if(e>2500 & e<=3000){res=fread(paste0("objformat_2501_3000.txt"),sep="\t",header=T)}
if(e>3000 & e<=3500){res=fread(paste0("objformat_3001_3500.txt"),sep="\t",header=T)}
if(e>3500 & e<=4000){res=fread(paste0("objformat_3501_4000.txt"),sep="\t",header=T)}
if(e>4000 & e<=4500){res=fread(paste0("objformat_4001_4500.txt"),sep="\t",header=T)}

res$vec=res[,paste0(liste[e]),with=FALSE]
bidon=res %>% dplyr::group_by(site) %>% dplyr::summarise(ndata_site=sum(vec),nyr=n_distinct(year))
bidon=subset(bidon,ndata_site>0)
res=res[res$site %in% bidon$site,]
res$site=as.factor(res$site)
res$year=as.numeric(as.character(res$year))-1951

lili=na.omit(unique(res$region))
resultf=NULL
for(i in 1:length(lili)){
res2=res %>% filter(region==lili[i])
if(nrow(bidon)>=2){
model=glmmTMB(vec~as.numeric(year)+log(tot)+(1|site),data=res2,family=binomial,
control=glmmTMBControl(optCtrl = list(iter.max=100000, eval.max=100000)))
modelnul=glmmTMB::glmmTMB(vec~1,data=res2,family=binomial)
result=as.data.frame(summary(model)$coeff$cond)
result$Esp=liste[e]
result$nsite=length(unique(res2$site))
result$varsite=summary(model)$varcor$cond$site[1]
result$rsquar=1-(logLik(model)/logLik(modelnul))
result$param=rownames(result)
result$region=lili[i]
result$nann=length(unique(res2$year))
}else{
result=as.data.frame(matrix(NA,1,4))
names(result)=c("Estimate","Std. Error","z value","Pr(>|z|)")
result$Esp=liste[e]
result$nsite=nrow(bidon)
result$varsite=NA
result$rsquar=NA
result$param=NA
result$region=lili[i]
result$nann=length(unique(res2$year))
}
resultf=rbind(resultf,result)
}

#names(result2)=paste0("bam_",names(result2))
#result=cbind(result,result2)

#resultMCMC=rbind(result,result)
#print(e)
return(resultf)
}
write.table(resultMCMC,paste("resultsMCMC_lin_zones_",deb,end,".txt",sep="_"),sep="\t",row.names=F)

closeCluster(cl)
#mpi.quit()
