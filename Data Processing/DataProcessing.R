###############################################################################
###############################################################################
## File: DataProcessing.R                                                    ##
## Author: Sanjeev V Namjoshi                                                ##
## Lab/Advisor: Kimberly Raab-Graham | University of Texas at Austin         ##
## Date: 1/28/14 | Revised for GitHub: 8/10/15                               ## 
## Description: R code used to process data obtained from LC-MS/MS for the   ##
## paper titled XXX published in Molecular Cell Proteomics XXX.              ##
## Data files: Numerical data and protein names from Supplemental Table 1A.  ##
## Copy cells B12:T756 from the Excel file to replicate this script exactly. ##
###############################################################################
###############################################################################

### Load packages ###
library(plyr)

### Load data from working directory ###
raw <- read.csv("RawData.csv", header = TRUE, stringsAsFactors = FALSE)

### Cosmetic processing for clarity ###
names(raw) <- c("Protein", "LD1", "LD2", "LD3", 
								           "LR1", "LR2", "LR3",
								           "PD1", "PD2", "PD3",
								           "PR1", "PR2", "PR3",
								           "SD1", "SD2", "SD3",
								           "SR1", "SR2", "SR3")

rownames(raw) <- raw[, "Protein"]
raw <- raw[, -1]
raw[raw == 0] <- NA   # Changing NAs to zero to help with processing downstream

### Split table into separate data frames for rapa/DMSO ###
lysDMSO <- raw[, 1:3]
lysRapa <- raw[, 4:6]
psdDMSO <- raw[, 7:9]
psdRapa <- raw[, 10:12]
solDMSO <- raw[, 13:15]
solRapa <- raw[, 16:18]

### Calculate number of 0 values per protein in each data frame ###

### Then remove proteins from the data frame if they have a 0 values for 2 ### 
### or 3 of 3 replicates
lysDMSO <- lysDMSO[which(apply(lysDMSO, 1, function(x) sum(is.na(x))) <= 1), ]
lysRapa <- lysRapa[which(apply(lysRapa, 1, function(x) sum(is.na(x))) <= 1), ]
psdDMSO <- psdDMSO[which(apply(psdDMSO, 1, function(x) sum(is.na(x))) <= 1), ]
psdRapa <- psdRapa[which(apply(psdRapa, 1, function(x) sum(is.na(x))) <= 1), ]
solDMSO <- solDMSO[which(apply(solDMSO, 1, function(x) sum(is.na(x))) <= 1), ]
solRapa <- solRapa[which(apply(solRapa, 1, function(x) sum(is.na(x))) <= 1), ]

# This produces a filtered protein list:
## lysDMSO = 605 proteins
## lysRapa = 593 proteins
## psdDMSO = 570 proteins
## psdRapa = 570 proteins
## solDMSO = 633 proteins
## solRapa = 641 proteins

### Find row means for all fractions/treatments ###
meanLysDMSO <- rowMeans(lysDMSO, na.rm = TRUE)
meanLysRapa <- rowMeans(lysRapa, na.rm = TRUE)
meanPsdDMSO <- rowMeans(psdDMSO, na.rm = TRUE)
meanPsdRapa <- rowMeans(psdRapa, na.rm = TRUE)
meanSolDMSO <- rowMeans(solDMSO, na.rm = TRUE)
meanSolRapa <- rowMeans(solRapa, na.rm = TRUE)

### Combine rapa and DMSO into one data frame, where NA means the protein ###
### was not found in that fraction/treatment.

##############################################################################
############### The "mergeVec" function will combine the data ################
##############################################################################

mergeVec <- function(x, y) {
	require(plyr)
	
	input <- list(as.data.frame(t(x)), as.data.frame(t(y)))
	input <- t(do.call(rbind.fill, input))
	
	return(input)
}

##############################################################################
##############################################################################
##############################################################################

lysate <- mergeVec(meanLysDMSO, meanLysRapa)
psd <- mergeVec(meanPsdDMSO, meanPsdRapa)
soluble <- mergeVec(meanSolDMSO, meanSolRapa)

### Set the column names and change NA values back to 0 ###
colnames(lysate) <- c("DMSO", "Rapa")
colnames(psd) <- c("DMSO", "Rapa")
colnames(soluble) <- c("DMSO", "Rapa")

lysate[is.na(lysate)] <- 0
psd[is.na(psd)] <- 0
soluble[is.na(soluble)] <- 0

### Compute the fold change for each fraction ###

##############################################################################
### Function to compute log fold changes for fractions and "out-of-range" ###
##############################################################################

foldChange <- function(fraction) {
	output <- as.data.frame(fraction[, 2] / fraction[, 1])
	output <- log2(output)
	
	output[output == Inf, ] <- "HIGH"
	output[output == -Inf, ] <- "LOW"
	
	return(output)
}

##############################################################################
##############################################################################
##############################################################################

lysateFC <- foldChange(lysate)
psdFC <- foldChange(psd)
solubleFC <- foldChange(soluble)

### Merge all lists together ###
fullList <- as.data.frame(mergeVec(lysateFC, psdFC))
fullList <- as.data.frame(mergeVec(fullList, solubleFC))

### Order the list alphabetically, put protein names in first column, and ###
### name the columns ###
fullList <- fullList[ order(row.names(fullList)), ]
fullList <- cbind(row.names(fullList), fullList)
row.names(fullList) <- 1:length(fullList[, 1])
colnames(fullList) <- c("Gene", "Lysate", "PSD", "Soluble")

fullList

### The final dataframe create by this script ("fullList") or some of the 
### precursor dataframes are used for all analysis downstream. ###

### EOF ###