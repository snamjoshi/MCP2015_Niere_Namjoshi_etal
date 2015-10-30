# Data Transformations ----------------------------------------------------

# Required data:

fullList <- read.csv("fullList.csv", header = TRUE)

### oorToNa()
# Purpose: Converts all instances of "LOW" or "HIGH" in the data to NA
# Usage: To just compare genes but ignore out-of-range genes

oorToNa <- function(geneList) {
	geneList[geneList == "LOW"] <- NA
	geneList[geneList == "HIGH"] <- NA
	geneList[geneList == "<NA>"] <- NA
	
	return(geneList)
}

### oorToValue()
# Purpose: Convert all instances of "LOW" or "HIGH" to -1 and 1 
# Usage: Calculations can now take the out-of-range genes into account
oorToValue <- function(geneList) {
	geneList[geneList == "LOW"] <- -1
	geneList[geneList == "HIGH"] <- 1
	geneList[geneList == "<NA>"] <- NA
	
	return(geneList)
}

### valueToOor()
# Purpose: Reverses the oorToValue function. Changes -1 or 1 back to "LOW" and "HIGH"
# Usage: After calculations, convert back to "LOW" and "HIGH" notation
valueToOor <- function(geneList) {
	geneList[geneList == -1] <- "LOW"
	geneList[geneList == 1] <- "HIGH"
	
	return(geneList)
}

### msDataNumeric()
# Purpose: Changes the msData set fraction values from factor to numeric. Must convert out-of-range values to -1/1 along the way.
# Usage: Data is converted to numeric because factor computations don't work well
# Dependencies: oorToValue()
msDataNumeric <- function(geneList) {
	
	msData <- geneList
	
	# Convert each column to a character
	msData$Lysate <- as.character(msData$Lysate)
	msData$PSD <- as.character(msData$PSD)
	msData$Soluble <- as.character(msData$Soluble)
	
	# Convert out-or-range genes to -1/1
	msData <- oorToValue(msData)
	
	# Convert columns to numeric
	msData$Lysate <- as.numeric(msData$Lysate)
	msData$PSD <- as.numeric(msData$PSD)
	msData$Soluble <- as.numeric(msData$Soluble)
	
	return(msData)
}

# Fold change search tools ------------------------------------------------

ent <- read.csv("entrez.csv", header = TRUE)

# findFoldChange()
# Purpose: Find fold changes for selected set of genes
findFoldChange <- function(data) {
	foldChange <- read.csv("fullList.csv", header = TRUE, na.strings = "<NA>")
	output <- foldChange[foldChange$Gene %in% data, ]
	return(output)
}

findFoldChange(fullList[1:10, 1])

# Arguments: Takes a list of genes (copied from clipboard) and the fullList variable.
# Output: Returns the rows in lysate, PSD, and soluble matching the genes in the input list
findFoldChanges <- function(fullList) {
	search <- readClipboard()
	search <- paste(substring(search,1,1), tolower(substring(search,2)), sep = "")
	data <- fullList[fullList$Gene %in% search,]
	return(data)
}

# Arguments: Takes a list of genes (copied from clipboard) and the fullList variable.
# Output: Returns the genes in the full list matching the input list
matchGenes <- function(fullList) {
	search <- readClipboard()
	search <- paste(substring(search,1,1), tolower(substring(search,2)), sep = "")
	data <- fullList[fullList$Gene %in% search, "Gene", drop = FALSE]
	return(data)
}

# Converts gene name to entrez data then returns gene info from NCBI
GeneSummary <- function(Gene) {
	require("NCBI2R")
	
	entr <- ent[ent$From %in% Gene, "To"]
	entr <- GetGeneInfo(entr)
	return(entr)
}

# Returns mean and sd DMSO values for the three rats for a particular gene. Can specifiy which rats to include in the row mean by modifying the vector passed to dmsoRats.
DMean <- function(fraction, gene, dmsoRats) {
	countMean <- rowMeans(fraction[fraction$Gene == gene, dmsoRats])
	countSD <- sd(fraction[fraction$Gene == gene, dmsoRats])
	data <- data.frame(gene = gene, mean = countMean, SD = countSD)
	return(data)
}

# Returns mean and sd rapa values for the three rats for a particular gene. Can specifiy which rats to include in the row mean by modifying the vector passed to rapaRats.
RMean <- function(fraction, gene, rapaRats) {
	countMean <- rowMeans(fraction[fraction$Gene == gene, rapaRats])
	countSD <- sd(fraction[fraction$Gene == gene, rapaRats])
	data <- data.frame(gene = gene, mean = countMean, SD = countSD)
	return(data)
}

# Returns the fold change (rapa/DMSO) for a particular gene. Can specify which rats to include in the row mean by modifying the vector passed to dmsoRats or rapaRats. Can also specify whether or not to take the log2 for the final value.
foldChange <- function(fraction, gene, dmsoRats, rapaRats, log) {
	countMeanD <- rowMeans(fraction[fraction$Gene == gene, dmsoRats])
	countMeanR <- rowMeans(fraction[fraction$Gene == gene, rapaRats])
	foldChange <- data.frame(gene = gene, foldChange = countMeanR/countMeanD)
	
	if(log == TRUE) {
		foldChange <- data.frame(gene = gene, foldChange = log(countMeanR/countMeanD, 2))
		return(foldChange)
	} else {
		return(foldChange)
	}
}

# Other processing functions ----------------------------------------------

# uniqueProt()
# Purpose: Function that returns list of proteins that are unique to their fraction.
uniqueProt <- function(proteinList){
	proteinList[proteinList == "HIGH"] <- NA
	proteinList[proteinList == "LOW"] <- NA
	
	return(proteinList[which(apply(t(apply(proteinList[,-1], 1, is.na)), 1, sum) == 2), ])
}

uniqueProt(fullList)

### bidirec()
# Purpose: Returns which genes change bidirectionally in the mass spec data 
# Dependencies: run msDataNumeric() first
bidirec <- function(msData) {
	# Assign NA to values between -1 and -1
	msData[ifelse(msData >-0.1 & msData < 0.1, TRUE, FALSE)] <- NA
	
	# Create three new columns   
	match <- data.frame(L_P = msData$Lysate*msData$PSD, P_S = msData$PSD*msData$Soluble, L_S = msData$Lysate*msData$Soluble)
	
	# If the are positive will be FALSE
	match <- ifelse(match > 0, NA, TRUE)
	
	return(match)
}

fullListNumeric <- msDataNumeric(fullList)
cbind(fullList, as.data.frame(bidirec(fullListNumeric)))
