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

setwd(dir="C:/Users/duche/Documents/3eme papier - climatic_debt/script/climatic_debt/data_and_raster_needed")
liste=fread("SCI_occupancy_traits_Europe.txt",header=T)
liste$life_span_cat[liste$life_span_cat==""]=NA
liste=subset(liste,nb_year>=20)
dim(liste)
mean(liste$trend)
sd(liste$trend)/sqrt(nrow(liste))
dim(subset(liste,trend<0 & Pr...z..<0.05))

liste$type[liste$type=="wild"]="non-invasive"
ggplot(data=liste,aes(x=trend,fill=type))+geom_density(alpha=0.6)+
geom_vline(xintercept=0,col="red")+scale_fill_manual(values=c("gold4","dodgerblue4"))+theme_bw()+
theme(panel.grid.minor=element_blank(),strip.background = element_blank(),
legend.title = element_blank(),plot.title=element_text(size=14,face="bold"))+
xlab("Occupancy trends (/year)")
liste=subset(liste,type!="invasive")


library(corrplot)
labs=data.frame(varia=c("ann_mean","prec","prec_min","min_temp","max_temp","prec_max",
"base_ann_mean","base_max_temp","base_min_temp","base_prec","base_prec_min","base_prec_max",
"nitrogen_moy","dep_poll","moisture_moy",
"aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche"),
vari=c("Trend of SCI annual T°C","Trend of SCI annual Prec.",
"Trend of SCI driest month Prec.","Trend of SCI min. T°C","Trend of SCI max. T°C",
"Trend of SCI wettest month Prec.","Historical annual T°C","Historical max. T°C","Historical min. T°C",
"Historical Annual Prec.","Historical driest month Prec.","Historical wettest month Prec.",
"Nutrients EIV","Pollinator dep.","Moisture EIV",
"Aquatic and wetland","Grassland","Heathland and tundra","Woodland",
"Farmland","Urban areas","Sparsely vegetated land"))

liste2=liste[,c("trend","trend_err","family","phylum","order","genus","class","max_temp","ann_mean","ann_mean","prec",
"prec_min","min_temp","max_temp","prec_max",
"base_ann_mean","base_max_temp","base_min_temp","base_prec","base_prec_min","base_prec_max",
"aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche",
"nitrogen_moy","dep_poll","moisture_moy","life_span_cat",
"Alpine","Anatolian","Arctic","Atlantic","BlackSea","Boreal","Continental","Mediterranean","Pannonian",
"Steppic")]
liste2[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]=
liste2[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]/
apply(liste2[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")],1,sum)
liste2[,c("ann_mean","prec","prec_min","min_temp","max_temp","prec_max",
"base_ann_mean","base_max_temp","base_min_temp","base_prec","base_prec_min","base_prec_max",
"nitrogen_moy","dep_poll","moisture_moy")]=as.data.frame(scale(
liste2[,c("ann_mean","prec","prec_min","min_temp","max_temp","prec_max",
"base_ann_mean","base_max_temp","base_min_temp","base_prec","base_prec_min","base_prec_max",
"nitrogen_moy","dep_poll","moisture_moy")]))
liste2=na.omit(liste2)

listeb=as.data.frame(liste2)
for(i in which(names(listeb) %in% labs$varia)){
names(listeb)[i]=as.character(labs$vari[labs$varia==names(listeb)[i]])}
png("matrice.png",width=1300,height=1200,res=140)
corrplot(cor(na.omit(listeb[,c("Trend of SCI annual T°C","Trend of SCI annual Prec.",
"Trend of SCI driest month Prec.","Trend of SCI min. T°C","Trend of SCI max. T°C",
"Trend of SCI wettest month Prec.","Historical annual T°C","Historical max. T°C","Historical min. T°C",
"Historical Annual Prec.","Historical driest month Prec.","Historical wettest month Prec.",
"Nutrients EIV","Pollinator dep.","Moisture EIV",
"Aquatic and wetland","Grassland","Heathland and tundra","Woodland",
"Farmland","Urban areas","Sparsely vegetated land")])), 
order = "hclust",method="number",number.cex=0.65)
dev.off();

model2=lmer(trend~(Urbain+Agricole+Heathland+aqua+Grasslands+Roche)+
(base_ann_mean+base_prec)+(ann_mean+prec)+moisture_moy+nitrogen_moy+dep_poll+
life_span_cat+(1|class/family),data=liste2,weights=1/trend_err)
print(warnings())

dat1=as.data.frame(summary(model2)$coeff)
dat1$varia=rownames(dat1)
dat1=rbind(dat1,dat1[1,],dat1[1,])
dat1$varia[rownames(dat1)=="(Intercept)1"]="Woodland"
dat1$varia[rownames(dat1)=="(Intercept)2"]="life_span_catannual"
dat1$Estimate2=NA
dat1$Estimate2[dat1$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")]=
dat1$Estimate[dat1$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")]+dat1$Estimate[dat1$varia=="(Intercept)"]
mat=as.matrix(vcov(model2))
dat1[dat1$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial"),"Std. Error"]=
sqrt(diag(mat)[c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")]+mat["(Intercept)","(Intercept)"]+
2*mat["(Intercept)",
c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")])
dat1$Estimate2[dat1$varia %in% c("Woodland","life_span_catannual")]=dat1$Estimate[dat1$varia %in% c("Woodland","life_span_catannual")]
dat1$Estimate2[dat1$varia %in% c("life_span_catannual","life_span_catintermediate","life_span_catperennial")]=
dat1$Estimate2[dat1$varia %in% c("life_span_catannual","life_span_catintermediate","life_span_catperennial")]+
sum(apply(liste2[,c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche")],2,mean)*
dat1$Estimate[dat1$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche")])
dat1$Model="LME"
dat1[,1:2]=dat1[,1:2]

MuMIn::r.squaredGLMM(model2)

#PHYLOGENIE
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
liste2=liste[!is.na(liste$code),]
liste2=liste[!duplicated(liste2$code),]

p3d=comparative.data(a,as.data.frame(liste2[!duplicated(liste2$code),]),code,vcv=TRUE,na.omit=F)
#svl<-setNames(p3d$data$trend,rownames(p3d$data))
library(phylosignal)
library(phylobase)
#p4d=phylo4d(as.phylo(p3d$phy),tip.data=svl)
#obj=phyloSignal(p4d, methods = c("Lambda"), reps=999, W = NULL)



liste2=liste[,c("trend","trend_err","family","phylum","order","max_temp","ann_mean","prec","prec_min",
"base_ann_mean","base_prec",
"aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche",
"nitrogen_moy","dep_poll","moisture_moy","life_span_cat",
"Alpine","Anatolian","Arctic","Atlantic","BlackSea","Boreal","Continental","Mediterranean","Pannonian",
"Steppic","code")]
liste2[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]=
liste2[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")]/
apply(liste2[,c("aqua","Grasslands","Heathland","Woodland","Agricole","Urbain","Roche")],1,sum)
liste2[,c("max_temp","ann_mean","prec","prec_min",
"base_ann_mean","base_prec",
"nitrogen_moy","dep_poll","moisture_moy")]=as.data.frame(scale(
liste2[,c("max_temp","ann_mean","prec","prec_min",
"base_ann_mean","base_prec",
"nitrogen_moy","dep_poll","moisture_moy")]))
liste2=na.omit(liste2)
liste2=liste2[!is.na(liste2$code),]
liste2=liste2[!duplicated(liste2$code),]
p3d=comparative.data(a,as.data.frame(liste2),code,vcv=TRUE)

modelt=pgls(trend~(Urbain+Agricole+Heathland+aqua+Grasslands+Roche)+
(base_ann_mean+base_prec)+(ann_mean+prec)+moisture_moy+nitrogen_moy+dep_poll+
life_span_cat,data=p3d,lambda='ML')


dat=as.data.frame(summary(modelt)$coeff)
dat$varia=rownames(dat)
dat=rbind(dat,dat[1,],dat[1,])
dat$varia[rownames(dat)=="(Intercept)1"]="Woodland"
dat$varia[rownames(dat)=="(Intercept)2"]="life_span_catannual"
dat$Estimate2=NA
dat$Estimate2[dat$varia %in% c("Woodland","life_span_catannual")]=dat$Estimate[dat$varia %in% c("Woodland","life_span_catannual")]
dat$Estimate2[dat$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")]=
dat$Estimate[dat$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial")]+dat$Estimate[dat$varia=="(Intercept)"]
dat[dat$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial"),"Std. Error"]=
sqrt(dat[dat$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche","life_span_catintermediate",
"life_span_catperennial"),"Std. Error"]^2+dat[dat$varia=="(Intercept)","Std. Error"]^2)
dat$Estimate2[dat$varia %in% c("life_span_catannual","life_span_catintermediate","life_span_catperennial")]=
dat$Estimate2[dat$varia %in% c("life_span_catannual","life_span_catintermediate","life_span_catperennial")]+
sum(apply(liste2[,c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche")],2,mean)*
dat$Estimate[dat$varia %in% c("Urbain","Agricole","Heathland","aqua","Grasslands","Roche")])
dat$Model="PGLS"
dat=rbind(dat[,-4],dat1)
dat$Estimate[!is.na(dat$Estimate2)]=dat$Estimate2[!is.na(dat$Estimate2)]
dat$lwr=dat$Estimate-1.96*dat[,2]
dat$upr=dat$Estimate+1.96*dat[,2]
dat$signi=">0.05"
dat$signi[dat$upr<0 | dat$lwr>0]="<0.05"
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
dat=merge(dat,labs,by="varia")
dat$vari=factor(dat$vari,c("                   Annual T°C","Annual Prec.","Prec. min.",
"Aquatic and wetland","Grassland","Heathland and tundra","Woodland",
"Farmland","Urban areas","Sparsely vegetated land",
"Salt EIV","Nitrogen EIV","Moisture EIV",
"Urbanity","Pollinator dep.","Invasive (ref=native)","annual","intermediate","perennial"))
dat$cate=factor(dat$cate,c("Climatic debt","Historical \n Climatic niche","Habitat affinity\n (for score = 1)",
"Species traits","Lifespan"))

pl1=ggplot(data=dat[dat$cate %in% c("Climatic debt","Historical \n Climatic niche","Species traits"),],
aes(x=vari,y=Estimate,col=signi,group=cate,shape=Model))+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge2(width=0.3))+theme_bw()+
geom_hline(yintercept=0,linetype="dashed")+theme(axis.text.x=element_text(angle=45, hjust=1,size=12),
axis.title.x=element_blank(),strip.background=element_blank(),plot.title=element_text(size=14,face="bold"),
legend.position="none")+
facet_grid(~cate,scales="free_x",space="free_x")+scale_color_manual(values=c("red","black"))+
labs(color="Estimate's p-value")+ylab("Standardised effects on occupancy trends")+
guides(colour=FALSE)+ggtitle("b")
pl1b=ggplot(data=dat[dat$cate %in% c("Habitat affinity\n (for score = 1)","Lifespan"),],
aes(x=vari,y=Estimate2,col=signi,group=cate,shape=Model))+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge2(width=0.3))+theme_bw()+
geom_hline(yintercept=0,linetype="dashed")+theme(axis.text.x=element_text(angle=45, hjust=1,size=12),
axis.title.x=element_blank(),strip.background=element_blank(),plot.title=element_text(size=14,face="bold"))+
facet_grid(~cate,scales="free_x",space="free_x")+scale_color_manual(values=c("red","black"))+
labs(color="Estimate's p-value")+ylab("Predicted occupancy trends")+
guides(colour=FALSE)+ggtitle("")

grid.arrange(pl1,pl1b,ncol=2)

liste$pval=">0.05"
liste$pval[liste$Pr...z..<0.05]="<0.05"
liste$pval=factor(liste$pval,c(">0.05","<0.05"))
p=ggplot(data=liste,aes(x=trend,fill=pval))+geom_histogram(col="black")+
geom_vline(xintercept=0,col="red")+geom_vline(xintercept=mean(liste$trend),col="gold3",
linetype="dashed")+scale_fill_manual(values=c("white","black"))+theme_bw()+
theme(panel.grid.minor=element_blank(),strip.background = element_blank(),
legend.title = element_blank(),legend.position="none",plot.title=element_text(size=14,face="bold"))+
xlab("Occupancy trends (logit(occupancy)/year)")+ylab("Number of species")+ggtitle("a")


##LA MEME CHOSE PAR ZONES:
setwd(dir="C:/Users/Francois/Documents/gabi plantes/script/climatic_debt/data_and_raster_needed")
liste=fread("SCI_occupancy_traits_region.txt",header=T)
liste$life_span_cat[liste$life_span_cat==""]=NA
liste=subset(liste,nb_year>=20)
liste=subset(liste,type!="invasive")

ggplot(data=liste[liste$region %in% c("Atlantic","Continental","Mediterranean","Alpine","Boreal")],
aes(x=trend,fill=region,col=region))+geom_density(alpha=0.4)+
geom_vline(xintercept=0,col="red")+theme_bw()+
theme(panel.grid.minor=element_blank(),strip.background = element_blank(),
legend.title = element_blank(),plot.title=element_text(size=14,face="bold"))+
xlab("Occupancy trends (/year)")


ggplot(data=liste[liste$region %in% c("Atlantic","Continental","Mediterranean","Alpine","Boreal")],
aes(x=ann_mean,fill=region,col=region))+geom_density(alpha=0.4)+
geom_vline(xintercept=0,col="red")+theme_bw()+
theme(panel.grid.minor=element_blank(),strip.background = element_blank(),
legend.title = element_blank(),plot.title=element_text(size=14,face="bold"))+
xlab("Température SCI")

ggplot(data=liste[liste$region %in% c("Atlantic","Continental","Mediterranean","Alpine","Boreal")],
aes(y=trend,x=ann_mean,col=region,col=region))+geom_point(alpha=0.4)+
geom_hline(yintercept=0,col="red")+theme_bw()+
theme(panel.grid.minor=element_blank(),strip.background = element_blank(),
legend.title = element_blank(),plot.title=element_text(size=14,face="bold"))+
xlab("Température SCI")+ylab("Occupancy trends")+facet_wrap(~region,scales="free")+stat_smooth(method="lm")


ggplot(data=liste[liste$region %in% c("Atlantic","Continental","Mediterranean","Alpine","Boreal")],
aes(x=prec,fill=region,col=region))+geom_density(alpha=0.4)+
geom_vline(xintercept=0,col="red")+theme_bw()+
theme(panel.grid.minor=element_blank(),strip.background = element_blank(),
legend.title = element_blank(),plot.title=element_text(size=14,face="bold"))+
xlab("Precipitation SCI trend")


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



pl1r=ggplot(data=datf[datf$cate %in% c("Climatic debt","Historical \n Climatic niche","Species traits"),],
aes(x=vari,y=Estimate,col=region,group=region,shape=Model))+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge2(width=0.3),alpha=0.8)+theme_bw()+
geom_hline(yintercept=0,linetype="dashed")+theme(axis.text.x=element_text(angle=45, hjust=1,size=12),
axis.title.x=element_blank(),strip.background=element_blank(),plot.title=element_text(size=14,face="bold"),
legend.position="none")+
facet_grid(~cate,scales="free_x",space="free_x")+
scale_color_manual(values=c("dodgerblue3","forestgreen","darkslategray3","lightsalmon4","gold3"))+
labs(color="Estimate's p-value")+ylab("Standardised effects on occupancy trends")+
guides(colour=FALSE)+ggtitle("")
pl1br=ggplot(data=datf[datf$cate %in% c("Habitat affinity\n (for score = 1)","Lifespan"),],
aes(x=vari,y=Estimate2,col=region,group=region,shape=Model))+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge2(width=0.3),alpha=0.8)+theme_bw()+
geom_hline(yintercept=0,linetype="dashed")+theme(axis.text.x=element_text(angle=45, hjust=1,size=12),
axis.title.x=element_blank(),strip.background=element_blank(),plot.title=element_text(size=14,face="bold"))+
facet_grid(~cate,scales="free_x",space="free_x")+
scale_color_manual(values=c("dodgerblue3","forestgreen","darkslategray3","lightsalmon4","gold3"))+
labs(color="")+ylab("Predicted occupancy trends")+
ggtitle("")

grid.arrange(pl1r,pl1br,ncol=2)

dat2=datf[datf$varia %in% c("ann_mean","prec","prec_min"),]
dat2$vari=as.character(dat2$vari)
dat2$vari[which(dat2$vari=="                   Annual T°C")]="Annual T°C"
dat2$vari=factor(dat2$vari,c("Annual T°C","Annual Prec."))
pl2=ggplot(data=dat2,aes(x=vari,y=Estimate,col=region,group=region,shape=Model))+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge2(width=0.3))+theme_bw()+
geom_hline(yintercept=0,linetype="dashed")+theme(axis.text.x=element_text(angle=0, hjust=0.5,size=10),
axis.title.x=element_blank(),strip.background=element_blank(),plot.title=element_text(size=14,face="bold"))+
labs(color="")+ylab("Standardised effects on occupancy trends")+
ggtitle("c")+
scale_color_manual(values=c("dodgerblue3","forestgreen","darkslategray3","lightsalmon4","gold3"))


library(gridExtra)
pdf("figure2.pdf",width=10, height=9)
grid.arrange(p,pl2,pl1,pl1b,layout_matrix=rbind(c(1,2),c(3,4)),heights=c(1,1.2))
dev.off();


supl1=ggplot(data=datf[datf$cate %in% c("Climatic debt","Historical \n Climatic niche","Species traits"),],
aes(x=vari,y=Estimate,col=region,group=region,shape=Model))+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge2(width=0.3))+theme_bw()+
geom_hline(yintercept=0,linetype="dashed")+theme(axis.text.x=element_text(angle=45, hjust=1,size=12),
axis.title.x=element_blank(),strip.background=element_blank(),plot.title=element_text(size=14,face="bold"))+
facet_grid(~cate,scales="free_x",space="free_x")+labs(color="")+ylab("Standardised effects on occupancy trends")+
ggtitle("")+
scale_color_manual(values=c("dodgerblue3","forestgreen","darkslategray3","lightsalmon4","gold3"))

supl2=ggplot(data=datf[datf$cate %in% c("Lifespan","Habitat affinity\n (for score = 1)"),],
aes(x=vari,y=Estimate,col=region,group=region,shape=Model))+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge2(width=0.3))+theme_bw()+
geom_hline(yintercept=0,linetype="dashed")+theme(axis.text.x=element_text(angle=45, hjust=1,size=12),
axis.title.x=element_blank(),strip.background=element_blank(),plot.title=element_text(size=14,face="bold"))+
facet_grid(~cate,scales="free_x",space="free_x")+labs(color="")+ylab("Standardised effects on occupancy trends")+
ggtitle("")+
scale_color_manual(values=c("dodgerblue3","forestgreen","darkslategray3","lightsalmon4","gold3"))

png("figure_S8.png",width=1500, height=600,res=120)
grid.arrange(supl1,supl2,ncol=2)
dev.off();






