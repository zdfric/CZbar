library(rgdal)
library(raster)
library(ape)
library(rhierbaps)
library(geosphere)
library(vegan)
library(pegas)
library(bold)
library(ConR)

#########
## read a map in SHP, containing information you want to keep for BAPS and haplotype networks
## for instance countrymap from Arcgis: https://hub.arcgis.com/datasets/esri::world-countries-generalized/explore 

# an example map with biogeography division
mapa<-readOGR("Gisy/Biogeo_simple.shp")

# European country map
mapa.x<-readOGR("E:/ESRI/ESRIDATA/EUROPE/country.shp")

plot(mapa,col="grey")
mapa
# if maps contain more data, select the one with your info (country names)
mapa<-mapa[,2]

########
## read saved data from BOLD
bold.data<-read.delim("otakarci_bold.txt"); grp<-"Papilionidae"

dim(bold.data)
names(bold.data)
# to see which markers are available
unique(bold.data$markercode)

# select only data for COU-5P
bold.data<-bold.data[bold.data$markercode=="COI-5P",]

# check the lenght of sequences
summary(nchar(bold.data$nucleotides))
hist(nchar(bold.data$nucleotides))

# use onlu sequences with >300 bp
bold.data<-bold.data[nchar(bold.data$nucleotides)>340,]
dim(bold.data)

## check if the sequences are aligned, otherwise export them, align and import again
bold.data.alignment<-as.alignment(as.data.frame(rbind(bold.data$nucleotides)))
bold.data.alignment$nam<-bold.data$processid
length(bold.data.alignment$seq)
length(bold.data.alignment$nam)

bold.data.fasta<-rbind(paste(">",bold.data.alignment$nam,sep=""),bold.data.alignment$seq)
write(bold.data.fasta,paste(grp,"bold_data_fasta.fas",sep="")) # writing and exporting the sequences

read.FASTA(paste(grp,"bold_data_fasta.fas",sep="")) # reading the sequences

dim(bold.data.alignment)
bold.data.alignment<-as.DNAbin(bold.data.alignment)
image(bold.data.fasta)


## limit the sequences to 658bb
bold.data$nucleotides<-substr(bold.data$aligned_nucleotides,1,658)
length(bold.data$nucleotides)

# check coordinates
summary(bold.data$lat)
summary(bold.data$lon)
dim(bold.data)
dim(bold.data[!is.na(bold.data$lat),])

# select only data with coordinates, delete duplicate records
bold.data<-bold.data[!is.na(bold.data$lat),]
dim(bold.data[!duplicated(bold.data$sampleid),])
bold.data<-bold.data[!duplicated(bold.data$sampleid),]
dim(bold.data)

write.table(bold.data,paste(grp,"_data_barcoding_filteredin.txt",sep=""),sep="\t")

#### subsetting to species
species<-sort(unique(bold.data$druh))
species

# empty dataframe for nucleotide diversity
nucldiv<-data.frame(matrix(NA,length(species),2))
names(nucldiv)<-c("species","nucl.diversity")

##### generating datasets and basic results species per species
for (y in 1:length(species))
{
spec<-species[y]
spec.data<-bold.data[bold.data$druh==spec,]
nucldiv$species[y]<-species[y]

pdf(paste(spec,".pdf",sep=""))

dim(spec.data)
plot(mapa,col="grey")
points(spec.data$lon,spec.data$lat,pch=16)
spec.points<-SpatialPointsDataFrame(cbind(spec.data$lon,spec.data$lat),spec.data,proj4string = CRS("+init=epsg:4326"))
points(spec.points,pch=16,col="red")
title(paste(spec," dataset distribution"))
  
spec.data$regions<-unlist(over(spec.points,mapa))
spec.data$regions[is.na(spec.data$regions)] <- "other"
# if regions are marked by strings, not numbers
regs<-unique(spec.data$regions)
regs.index<-seq.int(unique(spec.data$regions))
spec.data$regions.index<-NA
for (i in 1:length(regs)) {
  spec.data[spec.data$regions==regs[i],]$regions.index<-regs.index[i]
}

# to plot the regions on the map
plot(mapa,col="grey")
points(spec.data$lon,spec.data$lat,pch=16,col= spec.data$regions.index)
title(paste(spec," split to regions"))


sekvence<-matrix(NA,ncol = 2,nrow = nrow(spec.data))
  
sekvence[,1]<-paste(sub(" ","_",spec.data$druh),sub(" ","_",spec.data$processid),sep="_")
sekvence[,2]<-spec.data$nucleotides

traits<-sekvence
traits[,2]<-spec.data$regions

bapsy<-sekvence

sekvence<-matrix(paste(sekvence[,1],sekvence[,2],sep="     "),ncol = 1)
dim(sekvence)

traits2<-table(traits[,1],traits[,2])
match(traits[,1],rownames(traits2))
traits2<-traits2[match(traits[,1],rownames(traits2)),] # to make the same order as for sequences

rownames(traits2)
colnames(traits2)
write.table(traits2,"traits.txt",row.names = F,sep=",")
traits2<-read.delim("traits.txt")
traits<-cbind(traits,traits2)
colnames(traits)<-c("1","2","3")
traits<-matrix(paste(traits$`1`,traits$`3`),ncol=1)


# productions of nexus files for PopArt
spec.nexus<-rbind("#NEXUS","BEGIN DATA;",
                  paste("DIMENSIONS NTAX=",nrow(spec.data), " NCHAR=",nchar(spec.data$nucleotides[1]),";",sep=""),
                  "FORMAT INTERLEAVE DATATYPE=DNA MISSING=? GAP=-;",
                  "MATRIX","[coi]",
                  sekvence,";","end;",
                  
                  "Begin Traits;",
                  paste("Dimensions NTraits=",length(regs),";",sep=""),
                  "Format labels=yes missing=? separator=Comma;",
                  paste("TraitLabels ",toString(sort(regs)),";",sep=""),
                  "Matrix",
                  traits,
                  ";","end;"
                  )
write(spec.nexus,paste(spec,".nex",sep=""))

#### calculation of BAPS
baps.seq<-read.nexus.data(file=paste(spec,".nex",sep=""))
baps.seq<-as.DNAbin(baps.seq)
image(baps.seq)
  
# Nucleotide diversity
nucldiv$nucl.diversity[y]<-nuc.div(baps.seq)

# continuing with BAPS
baps.matice<-load_fasta(baps.seq)

try({baps.results<-hierBAPS(baps.matice,max.depth = 5, n.pops = 20,quiet = F)
baps.results
write.table(baps.results$partition.df,paste(spec,"_baps_results.txt",sep=""))

spec.data$baps<-baps.results$partition.df[,2] # pripadne zvysit cislo, kdyby bylo treba vyssi level
plot(mapa,col="grey") # global view
points(spec.data$lon,spec.data$lat,col=spec.data$baps,pch=16)
title(paste(spec," BAPS clusters identity"))

#svg(paste(spec,"_baps_map.svg"))
plot(mapa.svet,xlim=c(-31.26575,69.07032), ylim=c(32.39748,81.85737)) #European view, otherwise change the values
#plot(mapa,col="grey",add=T)
points(spec.data$lon,spec.data$lat,col=spec.data$baps,pch=16,cex=1.2)
title(paste(spec," BAPS clusters identity"))
#dev.off()
write.csv(cbind(spec.data$processid,spec.data$lat,spec.data$lon,spec.data$baps),paste("./baps_coords/",spec,"_baps_coordinates.csv",sep=""))

### nexus for baps inclusion into PopArt
bapsy[,2]<-spec.data$baps
bapsy2<-table(bapsy[,1],bapsy[,2])
match(bapsy[,1],rownames(bapsy2))
bapsy2<-bapsy2[match(bapsy[,1],rownames(bapsy2)),] #zmena poradi, aby bylo stejne jako u sekvenci

rownames(bapsy2)
colnames(bapsy2)
write.table(bapsy2,"bapsy.txt",row.names = F,sep=",")
bapsy2<-read.delim("bapsy.txt")
bapsy<-cbind(bapsy,bapsy2)
colnames(bapsy)<-c("1","2","3")
bapsy<-matrix(paste(bapsy$`1`,bapsy$`3`),ncol=1)

baps.nexus<-rbind("#NEXUS","BEGIN DATA;",
                  paste("DIMENSIONS NTAX=",nrow(spec.data), " NCHAR=",nchar(spec.data$nucleotides[1]),";",sep=""),
                  "FORMAT INTERLEAVE DATATYPE=DNA MISSING=? GAP=-;",
                  "MATRIX","[coi]",
                  sekvence,";","end;",
                  
                  "Begin Traits;",
                  paste("Dimensions NTraits=",length(unique(spec.data$baps)),";",sep=""),
                  "Format labels=yes missing=? separator=Comma;",
                  paste("TraitLabels ",toString(unique(spec.data$baps)),";",sep=""),
                  "Matrix",
                  bapsy,
                  ";","end;"
)
write(baps.nexus,paste(spec,"_baps.nex",sep=""))

})

#### Mantel (diversity by distance)
spec.geo.dist<-dist(cbind(spec.data$lon,spec.data$lat))
spec.gen.dist<-dist.dna(baps.seq, model="raw",pairwise.deletion = T,variance = T)
mantel.results<-mantel(spec.gen.dist,spec.geo.dist,permutations = 999,na.rm = T)
sink(paste(spec,"_Mantel_results.txt",sep=""))
print(mantel.results)
sink()
plot(spec.gen.dist,spec.geo.dist)
title(paste(spec," Isolation by distances"))

print(paste("y =",y,", species =",species[y],", out of ",length(species)," species",sep=""))
gc()

dev.off()

}
print(paste("y =",y,", species =",species[y],", out of ",length(species)," species",sep=""))
write.table(nucldiv,paste(grp,"_nucleotide_diversity.txt",sep=""),sep="\t")
dev.off()


save.image()





