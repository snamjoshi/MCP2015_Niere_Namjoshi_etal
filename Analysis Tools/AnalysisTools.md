# Analysis Tools for Niere and Namjoshi et al. 2015, MCP
Sanjeev V Namjoshi  
Revised for GitHub: October 29, 2015  

This file contains multiple R functions that may be useful in exploratory analysis or cleanup of the mass spectrometry data.

## Data transformation functions


```r
fullList <- read.csv("fullList.csv", header = TRUE)
```

### `oorToNa()`
**Purpose:** Converts all instances of "LOW" or "HIGH" in the data to NA


```r
oorToNa <- function(geneList) {
	geneList[geneList == "LOW"] <- NA
	geneList[geneList == "HIGH"] <- NA
	geneList[geneList == "<NA>"] <- NA
	
	return(geneList)
}

# Example:
oorToNa(fullList)
```

### `oorToValue()`
**Purpose:** Convert all instances of "LOW" or "HIGH" to -1 and 1 


```r
oorToValue <- function(geneList) {
	geneList[geneList == "LOW"] <- -1
	geneList[geneList == "HIGH"] <- 1
	geneList[geneList == "<NA>"] <- NA
	
	return(geneList)
}

# Example:
oorToValue(fullList)
```

### `valueToOor()`
**Purpose:** Reverses the oorToValue function. Changes -1 or 1 back to "LOW" and "HIGH"


```r
valueToOor <- function(geneList) {
	geneList[geneList == -1] <- "LOW"
	geneList[geneList == 1] <- "HIGH"
	
	return(geneList)
}

# Example:
valueToOor(fullList)
```

### msDataNumeric()
**Purpose:** Changes the msData set fraction values from factor to numeric. Must convert out-of-range values to -1/1 along the way.
**Dependencies:** `oorToValue()`


```r
msDataNumeric <- function() {
	
	msData <- read.csv("fullList.csv", header = TRUE)
	
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

# Example
msDataNumeric()
```

## Fold change search tools and gene query tools


```r
ent <- read.csv("entrez.csv", header = TRUE)
```

### `findFoldChange()`
**Purpose:** Find fold changes for selected set of genes


```r
findFoldChange <- function(data) {
	foldChange <- read.csv("fullList.csv", header = TRUE, na.strings = "<NA>")
	output <- foldChange[foldChange$Gene %in% data, ]
	return(output)
}

# Example: Returns fold changes of the first 10 genes
findFoldChange(fullList[1:10, 1])   
```

### `uniqueProt()`
**Purpose:** Function that returns list of proteins that are unique to their fraction.


```r
uniqueProt <- function(proteinList){
	proteinList[proteinList == "HIGH"] <- NA
	proteinList[proteinList == "LOW"] <- NA
	
	return(proteinList[which(apply(t(apply(proteinList[,-1], 1, is.na)), 1, sum) == 2), ])
}

# Example
uniqueProt(fullList)
```

### `findFoldChanges()`
**Purpose:** Returns the rows in lysate, PSD, and soluble matching the genes in the input list (copied from clipboard)


```r
findFoldChanges <- function(fullList) {
	search <- readClipboard()
	search <- paste(substring(search,1,1), tolower(substring(search,2)), sep = "")
	data <- fullList[fullList$Gene %in% search,]
	return(data)
}

# Example: Make sure gene list is copied to clipboard before running
findFoldChanges(fullList)
```

### `matchGenes()`
**Purpose:** Returns the genes in the full list matching the input list (copied from the clipboard). Used to see if a copied set of genes appeared in the mass spec data.


```r
matchGenes <- function(fullList) {
	search <- readClipboard()
	search <- paste(substring(search,1,1), tolower(substring(search,2)), sep = "")
	data <- fullList[fullList$Gene %in% search, "Gene", drop = FALSE]
	return(data)
}
# Example: Make sure gene list is copied to clipboard before running
matchGenes(fullList)
```

### `GeneSummary()`
**Purpose:** Converts gene name to entrez data then returns gene info from NCBI


```r
GeneSummary <- function(Gene) {
	require("NCBI2R")
	
	entr <- ent[ent$From %in% Gene, "To"]
	entr <- GetGeneInfo(entr)
	return(entr)
}

# Example:
GeneSummary("Snap25")
```

### `DMean()`
**Purpose:** Returns mean and sd DMSO values for the three rats for a particular gene. Can specifiy which rats to include in the row mean by modifying the vector passed to dmsoRats. Useful for getting a sense of the accuracy of the spectral count across replicates.


```r
DMean <- function(fraction, gene, dmsoRats) {
	countMean <- rowMeans(fraction[fraction$Gene == gene, dmsoRats])
	countSD <- sd(fraction[fraction$Gene == gene, dmsoRats])
	data <- data.frame(gene = gene, mean = countMean, SD = countSD)
	return(data)
}

# Example:
DMean(fraction = soluble, gene = "Camk2a", dmsoRats = c("SD1", "SD2", "SD3"))
```

### `RMean()`
**Purpose:** Returns mean and sd rapa values for the three rats for a particular gene. Can specifiy which rats to include in the row mean by modifying the vector passed to rapaRats.


```r
RMean <- function(fraction, gene, rapaRats) {
	countMean <- rowMeans(fraction[fraction$Gene == gene, rapaRats])
	countSD <- sd(fraction[fraction$Gene == gene, rapaRats])
	data <- data.frame(gene = gene, mean = countMean, SD = countSD)
	return(data)
}

# Example:
RMean(fraction = soluble, gene = "Camk2a", rapaRats = c("SR1", "SR2", "SR3"))
```

### `foldChange()`
**Purpose:** Returns the fold change (rapa/DMSO) for a particular gene. Can specify which rats to include in the row mean by modifying the vector passed to dmsoRats or rapaRats. Can also specify whether or not to take the log2 for the final value.


```r
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

# Examples:
foldChange(PSD, "Slc17a7", c("PD1", "PD2", "PD3"), c("PR1", "PR2", "PR3"),log = TRUE)
foldChange(lysate, "Slc17a7", c("LD1", "LD2", "LD3"), c("LR1", "LR2", "LR3"),log = TRUE)
foldChange(soluble, "Slc17a7", c("SD1", "SD2", "SD3"), c("SR1", "SR2", "SR3"),log = TRUE)
```

## Other analysis functions

### `uniqueProt()`
**Purpose:** Function that returns list of proteins that are unique to their fraction.

```r
uniqueProt <- function(proteinList){
	proteinList[proteinList == "HIGH"] <- NA
	proteinList[proteinList == "LOW"] <- NA
	
	return(proteinList[which(apply(t(apply(proteinList[,-1], 1, is.na)), 1, sum) == 2), ])
}

# Example:
uniqueProt(fullList)
```

### `bidirec()`
**Purpose:** Returns which genes change bidirectionally in the mass spec data 

**Dependencies:** run `msDataNumeric()` first.


```r
bidirec <- function(msData) {
	# Assign NA to values between -1 and -1
	msData[ifelse(msData >-0.1 & msData < 0.1, TRUE, FALSE)] <- NA
	
	# Create three new columns   
	match <- data.frame(L_P = msData$Lysate*msData$PSD, P_S = msData$PSD*msData$Soluble, L_S = msData$Lysate*msData$Soluble)
	
	# If the are positive will be FALSE
	match <- ifelse(match > 0, NA, TRUE)
	
	return(match)
}

# Examples:
fullListNumeric <- msDataNumeric()
cbind(fullList, as.data.frame(bidirec(fullListNumeric)))
```
