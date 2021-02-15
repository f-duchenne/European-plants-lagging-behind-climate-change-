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

setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")


#PHYLOGENIE
a=read.tree("DaPhnE_01.tre")
cod=read.table("code_canons.txt",sep="\t",header=T)
cod=subset(cod,rang=="SPECIES")
names(cod)[1]="noms_phylo"
cod=cod[!duplicated(cod$canon),]


##LA MEME CHOSE PAR ZONES:
setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")
liste=fread("SCI_occupancy_traits_region_before_1990.txt",header=T)
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
"Steppic","type","dep_poll")]
liste2=na.omit(liste2)
# liste2$type=factor(liste2$type,c("wild","invasive"))
liste2=subset(liste2,nb_year>=20)
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

bidon=bidon[!is.na(bidon$code),]
bidon=bidon[!duplicated(bidon$code),]
p3d=comparative.data(a,as.data.frame(bidon),code,vcv=TRUE)

model=pgls(trend~Urbain+Agricole+Heathland+aqua+Grasslands+Roche+
ann_mean+prec+moisture_moy+nitrogen_moy+dep_poll+
base_ann_mean+base_prec+life_span_cat,data=p3d,lambda='ML')
datb=as.data.frame(summary(model)$coeff)
datb$varia=rownames(datb)
datb=rbind(datb,datb[1,],datb[1,])
datb$varia[rownames(datb)=="(Intercept)1"]="Woodland"
datb$varia[rownames(datb)=="(Intercept)2"]="life_span_catannual"
datb$Estimate2=NA
datb$Estimate2[datb$varia %in% c("Woodland","life_span_catannual")]=datb$Estimate[datb$varia %in% c("Woodland","life_span_catannual")]
datb$Estimate2[datb$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")]=
datb$Estimate[datb$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")]+datb$Estimate[datb$varia=="(Intercept)"]
datb[datb$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial"),"Std. Error"]=
sqrt(datb[datb$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial"),"Std. Error"]^2+datb[datb$varia=="(Intercept)","Std. Error"]^2)
datb$Estimate2[datb$varia %in% c("life_span_catannual","life_span_catintermediate","life_span_catperennial")]=
datb$Estimate2[datb$varia %in% c("life_span_catannual","life_span_catintermediate","life_span_catperennial")]+
sum(apply(liste2[,c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche")],2,mean)*
datb$Estimate[datb$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche")])
datb$Model="PGLS"
datb$physig=paste(model$param[2],paste(model$param.CI$lambda$ci.val,collapse=" - "),sep=" - ")
datb=rbind(datb[,-4],dat1b)
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
"moisture_moy","human_dens","dep_poll","typeinvasive","life_span_catannual","life_span_catintermediate","life_span_catperennial",
"aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche"),
vari=c("                   Annual T°C","Annual Prec.","                   Annual T°C","Annual Prec.","Prec. min.",
"Salt EIV","Nitrogen EIV","Moisture EIV","Urbanity","Pollinator dep.","Invasive (ref=native)",
"annual","intermediate","perennial",
"Aquatic and wetland","Grassland","Heathland and tundra","Woodland",
"Farmland","Urban areas","Sparsely vegetated land"),
cate=c("Historical \n Climatic niche","Historical \n Climatic niche","Climatic debt","Climatic debt","Climatic debt",
"Species traits","Species traits","Species traits","Species traits","Species traits","Species traits",
"Lifespan","Lifespan","Lifespan","Habitat affinity\n (for score = 1)",
"Habitat affinity\n (for score = 1)","Habitat affinity\n (for score = 1)","Habitat affinity\n (for score = 1)",
"Habitat affinity\n (for score = 1)","Habitat affinity\n (for score = 1)","Habitat affinity\n (for score = 1)"))
datf=merge(datf,labs,by="varia")
datf$vari=factor(datf$vari,c("                   Annual T°C","Annual Prec.","Prec. min.",
"Aquatic and wetland","Grassland","Heathland and tundra","Woodland",
"Farmland","Urban areas","Sparsely vegetated land",
"Salt EIV","Nitrogen EIV","Moisture EIV",
"Urbanity","Pollinator dep.","Invasive (ref=native)","annual","intermediate","perennial"))
datf$cate=factor(datf$cate,c("Climatic debt","Historical \n Climatic niche","Habitat affinity\n (for score = 1)",
"Species traits","Lifespan"))



supl1=ggplot(data=datf[datf$cate %in% c("Climatic debt"),],
aes(x=vari,y=Estimate,col=region,group=region,shape=Model))+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge2(width=0.3))+theme_bw()+
geom_hline(yintercept=0,linetype="dashed")+theme(axis.text.x=element_text(angle=45, hjust=1,size=12),
axis.title.x=element_blank(),strip.background=element_blank(),plot.title=element_text(size=14,face="bold"))+
facet_grid(~cate,scales="free_x",space="free_x")+labs(color="")+ylab("Standardised effects on occupancy trends")+
ggtitle("")+
scale_color_manual(values=c("dodgerblue3","forestgreen","darkslategray3","lightsalmon4","gold3"))

png("figure_S6.png",width=1000, height=900,res=120)
supl1
dev.off();






