# CUROP2018
Integrating a multi-omic framework called MultiDataSet (MDS) with the graphical database platform Neo4j to demonstrate biological relationships between omics data.

## Project Aim
The aim of this project was to develop a tool to assist integration analysis of transcriptome data by providing an overview of the biological links between datasets such as microRNA and mRNA.

## Technologies
Project is created with: Neo4j 3.5, R version 3.5.2, and RNeo4J library 2018 (you can find repository here: https://github.com/nicolewhite/RNeo4j)

## NOTICE
Please note, this code is for demontration purposes only, as data used during this project is not provided 

## Data
Due to data privacy and protection, the datasets (GTF, GFF, and csv) used during this project cannot be shared

## Files
mRNA_format.R formats mRNA data and uploads it to a new MDS.
miRNA_format.R formats microRNA data and uplaods it to an existing MDS object containing mRNA data.
MDS_Neo4j.R sends data stored in MDS to Neo4j to create a graphical database of transcriptome data, represented as nodes, with biological relationships seen as edges between nodes. Additionally, Neo4j queries are given to extract certain nodes and the associated edges and target nodes.
