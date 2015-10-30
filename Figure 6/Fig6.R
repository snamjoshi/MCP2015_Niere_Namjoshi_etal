###############################################################################
###############################################################################
## File: Fig6.R                                                              ##
## Author: Sanjeev V Namjoshi (snamjoshi87@utexas.edu)                       ##
## Lab/Advisor: Kimberly Raab-Graham | University of Texas at Austin         ##
## Date: 2/28/14 | Revised for GitHub: 8/10/15                               ## 
## Description: R code used to generate figure 6 for the paper published     ##
##    in Molecular Cellular Proteomics, 2015.                                ##
## Data files: Table S7A protiens (cells A12:C615) and "fullList.csv" which  ##
##    is generated in DataProcessing.R                                       ##
###############################################################################
###############################################################################

library(VennDiagram)

### Import data ###
# Disease lists
disease <- read.csv("diseaseCandidates.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = "")

diseaseList <- list(Epilepsy = disease$Epilepsy,
										Alzheimers = disease$Alzheimers,
										ASD = disease$ASD)

diseaseList <- lapply(diseaseList, na.omit)

# Mass spec list
fullList <- read.csv("fullList.csv", header = TRUE)
fullList$Gene <- toupper(fullList$Gene)  # Supp table is in uppercase, easier comparison

# Epilepsy List
epilepsy <- diseaseList$'Epilepsy'
epilepsy <- fullList[fullList$Gene %in% epilepsy, ] # Filter to isolate genes in our dataset
#epilepsy <- epilepsy[!epilepsy$Gene %in% mito, ]

# Alzheimer's List
alz <- diseaseList$'Alzheimers'
alz <- fullList[fullList$Gene %in% alz, ] # Filter to isolate genes in our dataset
#alz <- alz[!alz$Gene %in% mito, ]

# ASD List
asd <- diseaseList$'ASD'
asd <- fullList[fullList$Gene %in% asd, ] # Filter to isolate genes in our dataset
#asd <- asd[!asd$Gene %in% mito, ]

##############################################################################
######### Creates all possible subsets for a three-way Venn diagram ##########
##############################################################################
#Creates all possible sets for a three-way venn diagram.
VennSet <- function(A, B, C)
{
	unionAB <- union(A, B)
	unionAC <- union(A, C)
	unionBC <- union(B, C)
	uniqueA <- setdiff(A, unionBC)
	uniqueB <- setdiff(B, unionAC)
	uniqueC <- setdiff(C, unionAB)
	intersAB <- setdiff(intersect(A, B), C)
	intersAC <- setdiff(intersect(A, C), B)
	intersBC <- setdiff(intersect(B, C), A)
	intersABC <- intersect(intersect(A, B), intersect(B, C))
	items <- list(A=uniqueA, B=uniqueB, C=uniqueC, AB=intersAB , AC=intersAC , BC=intersBC , ABC=intersABC)
	return(items)
}
##############################################################################
##############################################################################
##############################################################################

### Figure 6A (left) ###
# Find intersections and lengths
diseaseVenn <- VennSet(diseaseList$Epilepsy, diseaseList$Alzheimers, diseaseList$ASD)
lapply(diseaseVenn, length)

grid.newpage()
draw.triple.venn(area1 = 389, area2 = 209, area3 = 603, n12 = 13, n23 = 14, n13 = 59, 
								 n123 = 2, category = c("Epilepsy", "AD", "ASD"), lty = "blank", 
								 fill = c("brown1", "darkorchid2", "lightslateblue"))

### Figure 6B (right) ###
# Import data
OurDiseaseList <- read.csv("Diseases.csv", header = TRUE)
OurDiseaseList.list <- list(A = OurDiseaseList[,1], B = OurDiseaseList[,2], C = OurDiseaseList[,3])

disease.Venn <- VennSet(OurDiseaseList.list$A, OurDiseaseList.list$B, OurDiseaseList.list$C)
lapply(disease.Venn, length)

grid.newpage()
draw.triple.venn(area1 = 48, area2 = 44, area3 = 19, n12 = 3, n23 = 1, n13 = 5, 
								 n123 = 0, category = c("Epilepsy", "AD", "ASD"), lty = "blank", 
								 fill = c("rosybrown2", "plum", "skyblue1"))









