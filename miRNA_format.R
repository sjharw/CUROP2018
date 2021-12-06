###########################################################################
        ### Formatting miRNA data for MDS ###
###########################################################################

# Set you working directory to file containing miRNA data
setwd("C:/Users/YourWorkingDirectory")

# Load libraries
library(stringr) # to extract strings from gff file
library(tidyverse) # for dataframe formatting

######### Format miRNA feature data ################################################

# Read in file containing miRNA sequence information
fdata<- read.delim(file.choose(), header=F) # miRNA_fdata.GFF3 file

# Extract transcript ID and add as new column 
transcript_id <- str_match(fdata$V9,"ID=([A-Za-z0-9_]*);")[,2] # V9 is column containing transcript ID
fdata <- cbind(fdata, transcript_id)

# Extract mirBase_ID and add as new column 
mirBase_split<-str_split(fdata$V9, ";")
mirBase_sp<-lapply(mirBase_split, `[[`, 3)
mirBase_id<-sub("Name=", "", mirBase_sp)
fdata<- cbind(fdata, mirBase_id)

# Remove unwanted columns 
fdata$V2<-NULL
fdata$V3<-NULL
fdata$V6<-NULL
fdata$V8<-NULL
fdata$V9<-NULL

# Add columns names
colnames(fdata)[1]<- "chromosome"
colnames(fdata)[2]<- "start"
colnames(fdata)[3]<- "end"
colnames(fdata)[4]<- "strand"

# Save formatted feature data for miRNA
save(fdata, file="mirna_fdata.RData")

rm(list=ls()) # Remove everything from workspace


# Order fdata transcript id by adata
load(file="mirna_fdata.RData") # mirna_fdata
adata <-read.csv(file="miRNA_adata.csv",header = T,row.names = 1) #miRNA_adata
fdata_ordered <- matrix(nrow=nrow(adata), ncol=5)
colnames(fdata_ordered) <- c("mirBase_id", "chromosome","strand", "start", "end")
rownames(fdata_ordered) <- rownames(adata)
for(n in c(1:nrow(fdata_ordered))){
	trasncript_id<- rownames(adata)[n] 
      fdata_ordered[n,"mirBase_id"] <- as.character (fdata[match(trasncript_id,fdata[,"transcript_id"]),"mirBase_id"])
	fdata_ordered[n,"chromosome"] <- as.character (fdata[match(trasncript_id,fdata[,"transcript_id"]),"chromosome"])
	fdata_ordered[n,"strand"] <- as.character (fdata[match(trasncript_id,fdata[,"transcript_id"]),"strand"])
      fdata_ordered[n,"start"] <- as.numeric (fdata[match(trasncript_id,fdata[,"transcript_id"]),"start"])
	fdata_ordered[n,"end"] <- as.numeric (fdata[match(trasncript_id,fdata[,"transcript_id"]),"end"])
}

fdata<- as.data.frame(fdata_ordered)

save(fdata, file="mirna_fdata.RData")

###########################################################################
        ### Add miRNA data to MDS ###
###########################################################################

# Load libraries
library(MultiDataSet)

# Load mRNA assay and feature data
rm(list=ls()) # Remove everything from workspace
load("mirna_fdata.RData") # load miRNA feature data (fdata)
adata <-read.csv(file="miRNA_adata.csv",header = T,row.names = 1) # load miRNA assay data (adata)

# Create expression set data from assay data
eset_m3 <- new("ExpressionSet", exprs = as.matrix(adata))
fData(eset_m3) <- fdata # add feature data
save(eset_m3, file="eset_m3.RData")

# Add fdata and adata for miRNA to existing MDS
load("MDS.RData") # load MDS
MDS<- add_rnaseq(MDS, eset_m3, dataset.name="miRNA", warnings = FALSE, overwrite = TRUE)
save(MDS, file="MDS.RData")

