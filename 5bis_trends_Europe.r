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

if(nrow(bidon)>=2){
res$site=as.factor(res$site)
res$year=as.numeric(as.character(res$year))-1951
model=glmmTMB(vec~as.numeric(year)+log(tot)+(1|site),data=res,family=binomial,
control=glmmTMBControl(optCtrl = list(iter.max=100000, eval.max=100000)))
modelnul=glmmTMB::glmmTMB(vec~1,data=res,family=binomial)
result=as.data.frame(summary(model)$coeff$cond)
result$Esp=liste[e]
result$nsite=nrow(bidon)
result$varsite=summary(model)$varcor$cond$site[1]
result$rsquar=1-(logLik(model)/logLik(modelnul))
result$param=rownames(result)
}else{
result=as.data.frame(matrix(NA,1,4))
names(result)=c("Estimate","Std. Error","z value","Pr(>|z|)")
result$Esp=liste[e]
result$nsite=nrow(bidon)
result$varsite=NA
result$rsquar=NA
result$param=NA
}

#names(result2)=paste0("bam_",names(result2))
#result=cbind(result,result2)

#resultMCMC=rbind(result,result)
#print(e)
return(result)
}
write.table(resultMCMC,paste("resultsMCMC_lin",deb,end,".txt",sep="_"),sep="\t",row.names=F)

closeCluster(cl)
#mpi.quit()
