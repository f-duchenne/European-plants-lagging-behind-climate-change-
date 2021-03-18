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
library(optimx)
library(gridExtra)


setwd(dir="C:/Users/Francois/Documents/gabi plantes/data")
liste=fread("SCI_occupancy_traits_Europe.txt",header=T)
liste=subset(liste,type!="invasive")
stat=fread("StatutUICNEspFr.txt",header=T)
liste=merge(liste[,c("Latin name","trend")],stat,by.x="Latin name",by.y="Esp")
boxplot(trend~cat_iucn_europ,data=liste)

setwd(dir="C:/Users/Francois/Documents/gabi plantes/data")
liste=fread("SCI_occupancy_traits_regions.txt",header=T)
liste=subset(liste,type!="invasive")

liste$trend=liste$Estimate
liste$trend_err=liste$Std..Error
library(ggExtra)
pl1=ggplot(data=as.data.frame(liste2),aes(x=min_year,y=trend,group=min_year))+geom_boxplot()+geom_point(col=NA)+
stat_smooth()+xlab("First year with data of the focal species")+ylab("Occupancy trend")+ggtitle("a")
pl1=ggMarginal(pl1, type="histogram",margins ="x")
pl2=ggplot(data=as.data.frame(liste2),aes(x=nb_year,y=trend,group=nb_year))+geom_boxplot()+geom_point(col=NA)+
stat_smooth()+xlab("Number of years with data of the focal species")+ylab("Occupancy trend")+ggtitle("b")
pl2=ggMarginal(pl2, type="histogram",margins ="x")
grid.arrange(pl1,pl2,ncol=2)

#####################################################################
setwd(dir="C:/Users/duche/Documents/3eme papier - climatic_debt/script/climatic_debt/data_and_raster_needed")
liste=fread("SCI_occupancy_traits_region.txt",header=T)
liste$life_span_cat[liste$life_span_cat==""]=NA
liste=subset(liste,type!="invasive")
lili=c("Atlantic","Continental","Mediterranean","Alpine","Boreal")
liste=liste[liste$region %in% lili,]
pl1=ggplot(data=as.data.frame(liste),aes(x=nb_year,fill=region,color=region))+geom_density(alpha=0.5)+
xlab("Number of years with data of the focal species")+ylab("density")+
scale_color_manual(values=c("dodgerblue3","forestgreen","darkslategray3","lightsalmon4","gold3"))+
scale_fill_manual(values=c("dodgerblue3","forestgreen","darkslategray3","lightsalmon4","gold3"))+ggtitle("a")



boxplot(base_ann_mean~region,data=liste)
Rmisc::group.CI(base_ann_mean~region,data=liste)
tab=as.data.frame(expand.grid(lili,lili))
tab=subset(tab,Var1!=Var2)
tab$simi=NA
for(i in 1:nrow(tab)){
bidon=subset(liste,region==tab$Var1[i])
bidon1=subset(liste,region==tab$Var2[i])
tab$simi[i]=length(bidon$speciesKey[bidon$speciesKey %in% bidon1$speciesKey])/length(bidon$speciesKey)
}


##LA MEME CHOSE PAR ZONES:
veco=seq(15,50,5)
for(ii in veco){
liste=fread("SCI_occupancy_traits_region.txt",header=T)
liste$life_span_cat[liste$life_span_cat==""]=NA
liste=subset(liste,type!="invasive")

a=read.tree("DaPhnE_01.tre")
cod=read.table("code_canons.txt",sep="\t",header=T)
cod=subset(cod,rang=="SPECIES")
names(cod)[1]="noms_phylo"
cod=cod[!duplicated(cod$canon),]
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
head(liste[duplicated(liste$code),])

liste2=liste[,c("trend","trend_err","code","nsite","nb_year","region","species","family","phylum",
"order","class","max_temp","ann_mean","prec","prec_min","base_ann_mean","base_prec",
"aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche",
"nitrogen_moy","dep_poll","moisture_moy","life_span_cat",
"Alpine","Anatolian","Arctic","Atlantic","BlackSea","Boreal","Continental","Mediterranean","Pannonian",
"Steppic","type","dep_poll")]
liste2=na.omit(liste2)
# liste2$type=factor(liste2$type,c("wild","invasive"))
liste2=subset(liste2,nb_year>=ii)
summary(as.factor(liste2$region))
lili=c("Atlantic","Continental","Mediterranean","Alpine","Boreal")
datf=NULL
for(i in 1:length(lili)){
bidon=subset(liste2,paste0(region)==lili[i])
bidon[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]=
bidon[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]/
apply(bidon[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")],1,sum)

bidon[,c("max_temp","ann_mean","prec","prec_min",
"base_ann_mean","base_prec",
"nitrogen_moy","dep_poll","moisture_moy")]=as.data.frame(apply(
bidon[,c("max_temp","ann_mean","prec","prec_min",
"base_ann_mean","base_prec",
"nitrogen_moy","dep_poll","moisture_moy")],2,scale))
cor(bidon$prec,bidon$prec_min)

# bidon$ann_mean2=0
# bidon$prec2=0
# bidon$ann_mean2[bidon$ann_mean>0]=bidon$ann_mean[bidon$ann_mean>0]
# bidon$prec2[bidon$prec>0]=bidon$prec[bidon$prec>0]
model2=lmer(trend~Urbain+Agricole+Heathland+aqua+Grasslands+Roche+
ann_mean+prec+moisture_moy+nitrogen_moy+dep_poll+
base_ann_mean+base_prec+life_span_cat+
(1|family),data=bidon,weights=1/trend_err)

dat1b=as.data.frame(summary(model2)$coeff)
dat1b$varia=rownames(dat1b)
dat1b=rbind(dat1b,dat1b[1,],dat1b[1,])
dat1b$varia[rownames(dat1b)=="(Intercept)1"]="Woodland"
dat1b$varia[rownames(dat1b)=="(Intercept)2"]="life_span_catannual"
dat1b$Estimate2=NA
dat1b$Estimate2[dat1b$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")]=
dat1b$Estimate[dat1b$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")]+dat1b$Estimate[dat1b$varia=="(Intercept)"]
mat=as.matrix(vcov(model2))
dat1b[dat1b$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial"),"Std. Error"]=
sqrt(diag(mat)[c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")]+mat["(Intercept)","(Intercept)"]+
2*mat["(Intercept)",
c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")])
dat1b$Estimate2[dat1b$varia %in% c("Woodland","life_span_catannual")]=dat1b$Estimate[dat1b$varia %in% c("Woodland","life_span_catannual")]
dat1b$Estimate2[dat1b$varia %in% c("life_span_catannual","life_span_catintermediate","life_span_catperennial")]=
dat1b$Estimate2[dat1b$varia %in% c("life_span_catannual","life_span_catintermediate","life_span_catperennial")]+
sum(apply(liste2[,c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche")],2,mean)*
dat1b$Estimate[dat1b$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche")])
dat1b$Model="LME"
dat1b[,1:2]=dat1b[,1:2]
dat1b$physig=NA
datb=dat1b
datb$Estimate[!is.na(datb$Estimate2)]=datb$Estimate2[!is.na(datb$Estimate2)]
datb$region=lili[i]
datf=rbind(datf,datb)
print(i)
}

datf$lwr=datf$Estimate-1.96*datf[,2]
datf$upr=datf$Estimate+1.96*datf[,2]
datf$signi=">0.05"
datf$signi[datf$upr<0 | datf$lwr>0]="<0.05"
labs=data.frame(varia=c("base_ann_mean","base_prec","ann_mean","prec","prec_min",
"ellenberg_salt_baseflor","nitrogen_moy",
"moisture_moy","human_dens","dep_poll","typeinvasive","life_span_catintermediate","life_span_catperennial",
"aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche"),
vari=c("Annual T°C","Annual Prec.","Annual T°C","Annual Prec.","Prec. min.",
"Salt EIV","Nitrogen EIV","Moisture EIV","Urbanity","Pollinator dep.","Invasive (ref=native)",
"intermediate","perennial",
"Aquatic and wetlands","Grasslands","Heathlands and tundra","Woodlands",
"Agricultural","Urban","Cliffs and screes"),
cate=c("Historical \n Climatic niche","Historical \n Climatic niche","Climatic debt","Climatic debt","Climatic debt",
"Species traits","Species traits","Species traits","Species traits","Species traits","Species traits",
"Lifespan \n (ref=annual)","Lifespan \n (ref=annual)","Habitat affinity (ref=Woodland)",
"Habitat affinity (ref=Woodland)","Habitat affinity (ref=Woodland)","Habitat affinity (ref=Woodland)",
"Habitat affinity (ref=Woodland)","Habitat affinity (ref=Woodland)","Habitat affinity (ref=Woodland)"))
datf=merge(datf,labs,by="varia")
datf$vari=factor(datf$vari,c("Annual T°C","Annual Prec.","Prec. min.",
"Aquatic and wetlands","Grasslands","Heathlands and tundra","Woodlands",
"Agricultural","Urban","Cliffs and screes",
"Salt EIV","Nitrogen EIV","Moisture EIV",
"Urbanity","Pollinator dep.","Invasive (ref=native)","intermediate","perennial"))
datf$cate=factor(datf$cate,c("Climatic debt","Historical \n Climatic niche","Habitat affinity (ref=Woodland)",
"Species traits","Lifespan \n (ref=annual)"))
datf$threshold=ii
if(ii==15){datff=datf}else{datff=rbind(datff,datf)}
}


dat2=datff[datff$varia %in% c("ann_mean","prec","prec_min"),]
pl2=ggplot(data=dat2,aes(x=threshold,y=Estimate,col=region,group=region,shape=Model))+
facet_wrap(~vari+region,ncol=5)+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge2(width=0.3))+theme_bw()+
geom_hline(yintercept=0,linetype="dashed")+theme(axis.text.x=element_text(angle=0, hjust=0.5,size=10),
strip.background=element_blank(),plot.title=element_text(size=14,face="bold"),legend.position="none")+
labs(color="")+ylab("Standardised effects on occupancy trends")+
ggtitle("b")+xlab("Threshold on the number of years with data")+
scale_color_manual(values=c("dodgerblue3","forestgreen","darkslategray3","lightsalmon4","gold3"))

grid.arrange(pl1,pl2,ncol=1)
png("fig_referee.png",width=1000,height=1500,res=130)
grid.arrange(pl1,pl2,ncol=1)
dev.off();