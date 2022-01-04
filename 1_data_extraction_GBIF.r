#Script used with R.3.4.3. With more recent R versions, the rgbif loop for extraction is a bit different,
#mainly to precise study area, yeas window...
library(data.table)
library(rgbif)
memory.limit(size = 1e9)

setwd(dir="C:/Users/Francois/Documents/gabi plantes/data")
tab=fread("0053076-200221144449610.csv") #liste d'espèce du GBIF préalablement téléchargée
tab=subset(tab,taxonRank=="SPECIES" & taxonomicStatus=="ACCEPTED" & !is.na(speciesKey)) 
tab=subset(tab,phylum!="" & phylum!="Bryophyta" & phylum!="Chlorophyta" & phylum!="Rhodophyta" &
phylum!="Charophyta" & phylum!="Anthocerotophyta")
tab=tab[!duplicated(tab$speciesKey),] #remove duplicated GBIF ID (species correponding to the same one in the GBIF) to avoid double extraction
tab=subset(tab,numberOfOccurrences>=500)
tab=tab[,c("speciesKey","species","genus","family","order","class","phylum","numberOfOccurrences")]
tab$csum=cumsum(tab$numberOfOccurrences)
tab$part=cut(tab$csum,10) #cut the extraction in 10 parts

setwd(dir="D:/")

cle_li=list() #define empty list to put extraction keys
my_wkt=gbif_bbox2wkt(minx = -15, miny = 34, maxx = 34, maxy = 75) #define study area
geom_param <- paste("geometry", "within", my_wkt)

vec=unique(tab$part)
#EXTRCATION LOOP
for (i in 2:length(vec)) {
cle=occ_download(paste0("taxonKey = ",gsub(" ","",toString(tab$speciesKey[tab$part==vec[i]]))),
    'hasCoordinate = TRUE', 'hasGeospatialIssue = FALSE',"year >= 1950","year <= 2014",
	geom_param,
    user = "fduchenne", pwd="XXXXX",email = "XXXXXX")
cle_li[[i]]=cle[1]
while(occ_download_meta(cle_li[[i]][1])$status!="SUCCEEDED"){
Sys.sleep(60) 
} #wait until it is ready
#EXTRACT:
occ_download_get(cle_li[[i]][1],overwrite=T)
#need to reduce table dimension keep only (very) useful columns
#fwrite(b,paste0("extraction_part_",i,".txt"),sep="\t",row.names=F)
}

#DEZIPPER LE TOUT ET PRENDRE LES COLONNES D'INTERET:
setwd(dir="D:/")
vec=list.files()
cle_li=vec[grep("9610.zip",vec)]#define list of extraction keys
#EXTRCATION LOOP
for (i in 1:length(cle_li)) {
b=occ_download_import(as.download(path = cle_li[i], key = NULL),quote='',
select=c("gbifID","day","month","year","decimalLatitude","decimalLongitude","speciesKey")) 
#need to reduce table dimension keep only (very) useful columns
fwrite(b,paste0("extraction_part_",i,".txt"),sep="\t",row.names=F) #write the reduced extraction
}
