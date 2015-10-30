###############################################################################
###############################################################################
## File: SuppTables.R                                                        ##
## Author: Sanjeev V Namjoshi (snamjoshi87@utexas.edu)                       ##
## Lab/Advisor: Kimberly Raab-Graham | University of Texas at Austin         ##
## Date: 1/29/15 | Revised for GitHub: 8/10/15                               ## 
## Description: R code used to generate supplemental figures S1 and S2 for   ##
##    paper published in Molecular Cellular Proteomics, 2015                 ##
## Data Files: networkstats.csv and fullList.csv                             ##
###############################################################################
###############################################################################

# Table S2

# Load data
# Function reads a CSV file with uneven columns and coverts to a list instead of a dataframe
read.csv.list <- function(file, header) {
	outputList <- apply(read.csv(file, header = header), 2, list)
	outputList <- lapply(outputList, function(x) lapply(x, '[', which(complete.cases(x))))
	
	return(outputList)
}

network <- read.csv.list("networkStats.csv", TRUE)

# Data summary
meanLysate <- unlist(lapply(network, function(x) mean(as.numeric(unlist(x)))))
sdLysate <- unlist(lapply(network, function(x) sd(as.numeric(unlist(x)))))

summaryLysate <- data.frame(mean = round(meanLysate, 3), sd = round(sdLysate, 2))

summaryLysate

# Tests

# DMSO Lysate Analysis
dmsoAnalysis <- lapply(network[c(1,4,7)], function(x) wilcox.test(as.numeric(unlist(x)), alternative = "greater"))

# Rapa Lysate Analysis
rapaAnalysis <- lapply(network[c(2,5,8)], function(x) wilcox.test(as.numeric(unlist(x)), alternative = "greater"))

# Random Clustering Coefficent DMSO
randDmsoCC <- ks.test(as.numeric(unlist(network[[4]])), as.numeric(unlist(network[[6]])), alternative = "two.sided")

# Random Clustering Coefficent Rapa
randRapaCC <- ks.test(as.numeric(unlist(network[[5]])), as.numeric(unlist(network[[6]])), alternative = "two.sided")

# Random Neighborhood Connectivity DMSO
randDmsoNC <- ks.test(as.numeric(unlist(network[[7]])), as.numeric(unlist(network[[9]])), alternative = "two.sided")

# Random Neighborhood Connectivity Rapa
randRapaNC <- ks.test(as.numeric(unlist(network[[8]])), as.numeric(unlist(network[[9]])), alternative = "two.sided")

# Put everything into a table

summary <- data.frame(P.value = c(unname(unlist(sapply(dmsoAnalysis, "[", 3))), unname(unlist(sapply(rapaAnalysis, "[", 3))), randDmsoCC$p.value, randRapaCC$p.value, randDmsoNC$p.value, randRapaNC$p.value),
					 						row.names = c("DMSO Closeness Centrality", "DMSO Clustering Coefficient", "DMSO Neighborhood Connectivity", "Rapa Closeness Centrality", "Rapa Clustering Coefficient", "Rapa Neighborhood Connectivity", "Random DMSO Clustering Coefficent", "Random Rapa Clustering Coefficient", "Random DMSO Neighborhood Connectivity", "Random Rapa Neighborhood Connectivity"))

summary <- cbind(summary, LessThan0.001 = summary$P.value < 0.001)

summary

# Characteristic path length, network density, nodes, and edges calculated from Cytoscape built-in statistics.

# Number of proteins can be found in the DataProcessing.R file.

## Table S4

# Import data
fullList <- read.csv("fullList.csv", header = TRUE, stringsAsFactors = FALSE)

# First get data into correct form. This changes out-of-range genes to NAs for the calculations

oorToNa <- function(geneList) {
	geneList[geneList == "LOW"] <- NA
	geneList[geneList == "HIGH"] <- NA
	geneList[geneList == "<NA>"] <- NA
	
	return(geneList)
}

fullListNumeric <- oorToNa(fullList)

# Calculate
subset(fullListNumeric, as.numeric(PSD) > -0.1 & as.numeric(PSD) < 0.1, c("Gene", "PSD"))   # 76

# Table S5

subset(fullListNumeric, as.numeric(PSD) < -0.1, c("Gene", "PSD"))  # 159
subset(fullListNumeric, as.numeric(PSD) > 0.1, c("Gene", "PSD"))   # 166

# Table S6

as.vector(na.omit(fullList[fullList$PSD == "HIGH", "Gene"]))
as.vector(na.omit(fullList[fullList$PSD == "LOW", "Gene"]))

# A few genes were removed by hand when we realized the spectral count calculations may have been innacurate

# Table S7

# See figure 6 to generate these lists.

# Table S8

# See figure 6 to generate these lists.

# All other figures generated through either Cytoscape or DAVID.















