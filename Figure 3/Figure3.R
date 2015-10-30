###############################################################################
###############################################################################
## File: Fig3.R                                                              ##
## Author: Sanjeev V Namjoshi (snamjoshi87@utexas.edu)                       ##
## Lab/Advisor: Kimberly Raab-Graham | University of Texas at Austin         ##
## Date: 2/28/14 | Revised for GitHub: 8/10/15                               ## 
## Description: R code used to generate figure 3 for the paper published     ##
##    in Molecular Cellular Proteomics, 2015.                                ##
## Data files: Table S7A protiens (cells A12:C615) and "fullList.csv" which  ##
## is generated in DataProcessing.R and list of mitochontrial genes in       ##
## "mito.csv"                                                                ##
###############################################################################
###############################################################################

### Load libraries ###
library(gplots)

### Import data ###
# Disease lists
disease <- read.csv("diseaseCandidates.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = "")

diseaseList <- list(Epilepsy = disease$Epilepsy,
										Alzheimers = disease$Alzheimers,
										ASD = disease$ASD)

diseaseList <- lapply(diseaseList, na.omit)
diseaseList <- lapply(diseaseList, as.data.frame) # Removes na.action attributes

# Mass spec list
fullList <- read.csv("fullList.csv", header = TRUE)
fullList$Gene <- toupper(fullList$Gene)  # Supp table is in uppercase, easier comparison

# Mitochondrial gene list
mito <- read.csv("mito.csv", header = FALSE, stringsAsFactors = FALSE)
mito <- toupper(mito[ ,1])

### Figure 1A ### 
# Note for all figures that some proteins with unreliable spectral counts (high
# variability in values for the replicates) were removed from the final heatmap.
epilepsy <- diseaseList$'Epilepsy'
epilepsy <- fullList[fullList$Gene %in% epilepsy[ ,1], ] # Filter to isolate genes in our dataset
epilepsy <- epilepsy[!epilepsy$Gene %in% mito, ]
epilepsy <- sapply(epilepsy[ ,2:4], function(x) as.numeric(as.character(x)))
epilepsyMatrix <- data.matrix(epilepsy)

# Labels and out-of-range gene textured boxes added manually in Inkscape later
heatmap.2(
	as.matrix(epilepsyMatrix[ ,c(2,3)]),
	dendrogram = "none",
	Rowv = NULL,
	Colv = NULL,
	col = "heat.colors",
	trace = "none", 
	na.color = "grey",
	density.info = "none",
	key = TRUE,
	keysize = 0.5,
	labCol = FALSE
)

### Figure 1B ###
alz <- diseaseList$'Alzheimers'
alz <- fullList[fullList$Gene %in% alz[ ,1], ] # Filter to isolate genes in our dataset
alz <- alz[!alz$Gene %in% mito, ]
alz <- sapply(alz[ ,2:4], function(x) as.numeric(as.character(x)))
alzMatrix <- data.matrix(alz)

# Labels and out-of-range gene textured boxes added manually in Inkscape later
heatmap.2(
	as.matrix(alzMatrix[ ,c(2,3)]),
	dendrogram = "none",
	Rowv = NULL,
	Colv = NULL,
	col = "heat.colors",
	trace = "none", 
	na.color = "grey",
	density.info = "none",
	key = TRUE,
	keysize = 0.5,
	labCol = FALSE
)	

### Figure 1C ###
asd <- diseaseList$'ASD'
asd <- fullList[fullList$Gene %in% asd[ ,1], ] # Filter to isolate genes in our dataset
asd <- asd[!asd$Gene %in% mito, ]
asd <- sapply(asd[ ,2:4], function(x) as.numeric(as.character(x)))
asdMatrix <- data.matrix(asd)

# Labels and out-of-range gene textured boxes added manually in Inkscape later
heatmap.2(
	as.matrix(asdMatrix[ ,c(2,3)]),
	dendrogram = "none",
	Rowv = NULL,
	Colv = NULL,
	col = "heat.colors",
	trace = "none", 
	na.color = "grey",
	density.info = "none",
	key = TRUE,
	keysize = 0.5,
	labCol = FALSE
)	

### EOF ###