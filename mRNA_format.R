################################################################################################
                         ### Formatting mRNA data for MDS ###
#################################################################################################

# Set you working directory to file containing miRNA data
setwd("C:/YourWorkingDirectory")

# Load libraries
library(ballgown) # for gff read
library(stringr) # for string extraction
library(tidyverse) # for dataframe formatting
library(data.table) # lists to data frame

############################# Format mRNA assay data ################################################

# Read in mRNA count data
adata <-read.csv(file.choose(), header = TRUE) # mRNA_adata

# Set rownames as transcript_id column
adata <- adata %>% remove_rownames %>% column_to_rownames(var="Transcript_ID")

# Save mRNA assay data
save(adata, file="mrna_adata.RData")


############################# Format mRNA feature data ################################################
## This code generates a dataframe of mRNA Transcript IDs, Gene IDs, Chromosome start and end location ##

# Read in gtf mRNA feature data file
fdata<-gffRead(file.choose(), nrows = -1, verbose = FALSE) # mRNA_fdata GTF file

# Extract Trancript ID
transcript_id <- str_extract(fdata$attributes, "(NM|NR){1}_[0-9]{1,}_*[0-9]*")
fdata<- cbind(fdata, transcript_id)

# Extract Gene ID 
gene_id<-str_split(fdata$attributes, "; ")
gene_id<-lapply(gene_id, `[[`, 1)
gene_id<-sub(".* \"", "", gene_id)
gene_id<-gsub("\"","", gene_id)
fdata<- cbind(fdata, gene_id)

# Remove unwanted columns
fdata$attributes<-NULL
fdata$source <-NULL
fdata$score <-NULL
fdata$frame <-NULL

  ## Extract minimum start value and maximum end value of chromosomes

# Collapse columns based on duplicate transcript IDs
combine<-setDT(fdata)[, lapply(.SD, function(x) toString(na.omit(x))), by = transcript_id]

# Extract min start value for unique transcript id
start<- matrix(nrow=nrow(combine), ncol=1) # create df to add min start values to
start_min<-str_split(as.vector(combine$start), ", ")
for(n in c(1:length(start_min))){
	startval = min(as.numeric(start_min[[n]]))
      start[n,1] = startval
	print(startval)
}

# Extract max end value from combine
end<- matrix(nrow=nrow(combine), ncol=1)
end_max<-str_split(as.vector(combine$end), ", ")
for(n in c(1:length(end_max))){
	endval = max(as.numeric(end_max[[n]]))
      end[n,1] = endval
	print(endval)
}

# Chromosome extraction from combine
# Remove duplicated chr by extracting single chr from column 
chromosome<-str_split(as.vector(combine$seqname), ", ")
chromosome<-sapply(chromosome, `[[`, 1)

# Strand extraction from combine
strand<-str_split(as.vector(combine$strand), ", ")
strand<-sapply(strand, `[[`, 1)

# Extract gene name from combine
gene_name<-str_split(as.vector(combine$gene_id), ", ")
gene_name<-sapply(gene_name, `[[`, 1)

# Extract unique transcript ids from combine
transcript_id<-str_split(as.vector(combine$transcript_id), ", ")
transcript_id<-sapply(transcript_id, `[[`, 1)

# Combine extracted columns into new data frame
# This new df has min start and max end chromosome values
fdata <- data.frame(matrix(unlist(transcript_id), nrow=length(transcript_id), byrow=TRUE))
fdata ['gene_name'] = gene_name
fdata ['strand'] = strand
fdata ['start'] = start
fdata ['end'] = end
fdata ['chromosome'] = chromosome
names(fdata )[1] <- 'transcript_id'


# Save mRNA assay data
save(fdata, file="mrna_fdata.RData")

rm(list=ls()) # Remove everything from workspace

# Load mRNA formatted assay data (adata)
load("mrna_adata.RData")

# Load mRNA formatted assay data (adata)
load("mrna_fdata.RData")

# Order mrna transcript_id in fdata based on order in adata
# bit long-winded but only way I found that works
fdata_ordered <- matrix(nrow=nrow(adata), ncol=5)
colnames(fdata_ordered) <- c("gene_name", "chromosome", "strand", "start", "end")
rownames(fdata_ordered) <- rownames(adata)
for(n in c(1:nrow(fdata_ordered))){
	transcript_id <- rownames(adata)[n]
      fdata_ordered[n,"gene_name"] <- as.character(fdata[match(transcript_id, fdata[,1]),"gene_name"])
	fdata_ordered[n,"chromosome"] <- as.character(fdata[match(transcript_id, fdata[,1]),"chromosome"])
      fdata_ordered[n,"strand"] <- as.character(fdata[match(transcript_id, fdata[,1]),"strand"])
	fdata_ordered[n,"start"] <- as.numeric(fdata[match(transcript_id, fdata[,1]),"start"])
	fdata_ordered[n,"end"] <- as.numeric(fdata[match(transcript_id, fdata[,1]),"end"])
}

# Save formatted feature data for mRNA
save(fdata_ordered, file="mrna_fdata.RData")


################################################################################################
                         ### Add mRNA data to MDS ###
#################################################################################################

# Load libraries
library(MultiDataSet)

# Load mRNA assay and feature data
rm(list=ls()) # Remove everything from workspace
load("mrna_fdata.RData") # load mRNA feature data (fdata)
load("mrna_adata.RData") # load mRNA assay data (adata)
fdata<-as.data.frame(fdata_ordered) # Convert matrix to data frame for MDS


# Create expression set data from assay data
eset_m4 <- new("ExpressionSet", exprs = as.matrix(adata))
fData(eset_m4) <- fdata # add feature data
save(eset_m4, file="eset_m4.RData")

# Create MDS object and add fdata and adata
MDS<-createMultiDataSet()
MDS<- add_rnaseq(MDS, eset_m4, dataset.name="mRNA", warnings = FALSE, overwrite = TRUE)
save(MDS, file="MDS.RData")




