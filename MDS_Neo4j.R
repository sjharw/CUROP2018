###########################################################################
        ### Pushing MDS data to Neo4j ###
###########################################################################

# Load libraries
library(MultiDataSet)
library(RNeo4j)
library(igraph)

# Download and install Neo4j
# Start new project 
# Set password
# Make sure project is active

# Create a connection between R and Neo4j via the RNeo4j package.
graph <- startGraph("http://localhost:7474/db/data", username = "username", password = "password")

# Load MultiDataSet object containing mRNA and microRNA (miRNA)
load("MDS.RData") #MDS

# Extract mRNA and microRNA data fom MDS
fDmi<-fData(MDS[["rnaseq+miRNA"]])
fDm<-fData(MDS[["rnaseq+mRNA"]])


############### Making mRNA transcript and Gene name nodes #######################
############### And creating relationship between them ###########################

# Create Transcript ID nodes for mRNA 
mRNA<-lapply(1:nrow(fDm), function(x){createNode(graph, .label = as.character("mRNA"),
TranscriptID= rownames(fDm)[[x]])})
graph

# Create Gene Name nodes for mRNA 
GeneName<-lapply(1:nrow(fDm), function(x){createNode(graph, .label = as.character("GeneName"),
GeneName= fDm$Gene_name[[x]])})

# Create relationship between gene and transcript
dRelationship <-function(x){ createRel(
GeneName[[x]], "transcript",mRNA[[x]])}
lapply(1:nrow(fDm), dRelationship)

# Merge duplicate Gene Name nodes
cypher(graph, "
MATCH (n3:GeneName)
WITH n3.GeneName as GeneName, collect(n3) AS nodes
WHERE n3.GeneName = n3.GeneName AND size(nodes) >  1
CALL apoc.refactor.mergeNodes(nodes)
YIELD node
MATCH (n3) WHERE NOT labels(n3) DELETE n3
RETURN *
")

# Check its working as expected
# You should see pathway returned in active project opened in Neo4J browser
cypher(graph, "
MATCH (n3:GeneName)
MATCH (n4:mRNA)
MATCH path=(n3)-->(n4)
WHERE n3.GeneName= 'PUSL1'
RETURN n4.TranscriptID
")

############### Making microRNA transcript and mirBase ID nodes #######################
############### And creating relationship between them ###########################

# Create Transcript ID nodes for microRNA
microRNA<-lapply(1:nrow(fDmi), function(x){createNode(graph, .label = as.character("microRNA"),
TranscriptID= rownames(fDmi)[[x]])})

# Create mirBaseID nodes for microRNA
mirBaseID<-lapply(1:nrow(fDmi), function(x){createNode(graph, .label = as.character("mirBaseID"),
mirBaseID= fDmi$mirBase_ID[[x]])})

# Create relationship between transcript and mirBaseID
cRelationship <-function(x){ createRel(
microRNA[[x]], "transcript",mirBaseID[[x]])}
lapply(1:nrow(fDmi), cRelationship)

# Merge duplicate modes
cypher(graph, "
MATCH (n2:mirBaseID)
WITH n2.mirBaseID as mirBaseID, collect(n2) AS nodes
WHERE n2.mirBaseID  = n2.mirBaseID  AND size(nodes) >  1
CALL apoc.refactor.mergeNodes(nodes)
YIELD node
MATCH (n2) WHERE NOT labels(n2) DELETE n2
RETURN *
")

# Check its working
cypher(graph, "
MATCH (n1:microRNA)
MATCH (n2:mirBaseID)
MATCH path=(n1)-->(n2)
WHERE n1.TranscriptID= 'MI0000300'
RETURN n2.mirBaseID
")

############ Making relastionship between microRNA (mirBaseID) and mRNA (GeneName) ################

# Create relastionship between mirBaseID (from microRNA) and GeneName (from mRNA)
cypher(graph, "LOAD CSV WITH HEADERS FROM 'file:///spider.csv' as edges 
MATCH (n2:mirBaseID { mirBaseID: edges.mirBaseID })
MATCH (n3:GeneName { GeneName: edges.GeneName })
CREATE (n2)-[r:TARGET]->(n3)")



############# Extracting infromation from RNeo4j graph ##########################

  ## Within R

# Extract pathway for particular microRNA (mirBaseID) 
cypher(graph, "
MATCH (n1:microRNA)
MATCH (n2:mirBaseID)
MATCH (n3:GeneName)
MATCH (n4:mRNA)
MATCH path=(n1)-->(n2)-->(n3)-->(n4)
WHERE n2.mirBaseID= 'hsa-miR-9-5p'
RETURN n1.TranscriptID, n3.GeneName, n4.TranscriptID
")

  ## Within Neo4J command line

# Return node with certain property 
match (n) 
with n, [x in keys(n) WHERE n[x]='hsa-miR-9-5p'] as doesMatch
where size(doesMatch) > 0
return n

# Returns node with certain property, and its relationships
# Requires path and node location to be known
MATCH (n)
with n, [x in keys(n) WHERE n[x]='hsa-miR-9-5p'] as doesMatch
where size(doesMatch) > 0
MATCH p = ()-->(n)-->()-->() # You need to state location of node in the path
return p

# Returns node with certain property, and its relationships
# Does not require path to be known
MATCH (n)
with n, [x in keys(n) WHERE n[x]='MI0000113'] as doesMatch
where size(doesMatch) > 0
match path = ( (n)-[*]->(m) ) # In this case, you do not need to state the path or the location of the node in the path
return path
# PAth requires a direction, and therefore only returns paths in that direction
# Have yet to find function that returns any path associated with node regardless of direction
