###############################################################################
###############################################################################
## File: Fig2.R                                                              ##
## Author: Sanjeev V Namjoshi (snamjoshi87@utexas.edu)                       ##
## Lab/Advisor: Kimberly Raab-Graham | University of Texas at Austin         ##
## Date: 2/24/14 | Revised for GitHub: 8/10/15                               ## 
## Description: R code used to generate figure 2 for the paper published     ##
##    in Molecular Cellular Proteomics, 2015                                 ##
## Data files: "fullList.csv" variable generated in DataProcessing.R         ##
###############################################################################
###############################################################################

### Load libraries ###
library(ggplot2)
library(reshape2)
library(moments)

### Figure 2d ###
# Import data #
fullList <- read.csv("fullList.csv", header = TRUE)

# Make cosmetic changes and get data into long form #
fullList <- fullList[, -1]
fullList[fullList == "HIGH"] <- NA
fullList[fullList == "LOW"] <- NA
fullList <- apply(fullList, 2, as.numeric)
fullList <- melt(fullList)
fullList <- fullList[, -1]
colnames(fullList) <- c("Fraction", "Expression")

# Subset PSD fraction and graph #
histLongPSD <- subset(fullList, Fraction == "PSD")

binsize <- diff(range(na.omit(histLongPSD$Expression)))/50

a <- ggplot(histLongPSD, aes(x = Expression, fill = Fraction)) +
	geom_histogram(position = "identity", binwidth = binsize) +
	geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
	geom_vline(xintercept = 0.1375, linetype = "solid", color = "red", size = 0.3) +
	geom_vline(xintercept = -0.1375, linetype = "solid", color = "red", size = 0.3) +
	scale_fill_manual(values = c("black")) +
	scale_x_continuous(limits = c(-1.75, 1.75), breaks = seq(-2, 2, by = 0.5)) +
	xlab("Log2 Fold Change") +
	ylab("Number of Proteins") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
				panel.background = element_blank(), axis.line = element_line(color = "black"))

a + 
	theme(text = element_text(size = 15)) +
	theme(legend.position = "none") +
	annotate("rect", xmin = -0.136, xmax = 0.1375, ymin = 0, ymax = Inf, fill = "grey99", alpha = 0.9)

### Statistics in Figure 2d Legend ###
# Import data #
fullList <- read.csv("fullList.csv", header = TRUE)
fullListPSD <- fullList[ ,c(1,3)]
fullListPSD$PSD <- as.numeric(as.character(fullListPSD$PSD))

# Number of proteins calculations (excludes out-of-range)
sum(fullListPSD$PSD > 0, na.rm = TRUE) # 210 proteins greater than 0
sum(fullListPSD$PSD < 0, na.rm = TRUE) # 191 proteins less than 0
sum(complete.cases(fullListPSD$PSD)) # 401 total proteins

## Descriptive statistics ##
mean(fullListPSD$PSD, na.rm = TRUE) # Mean
sd(fullListPSD$PSD, na.rm = TRUE) # Standard deviation
var(fullListPSD$PSD, na.rm = TRUE) # Variance (CI calculated in GraphPad)
skewness(fullListPSD$PSD, na.rm = TRUE) # Skewness
kurtosis(fullListPSD$PSD, na.rm = TRUE) # Kurtosis

## Bartlett's test for comparison of variance ##
# Import Data and process #
fullList <- read.csv("fullList.csv", header = TRUE)
fullList <- fullList[, -1]
fullList[fullList == "HIGH"] <- NA
fullList[fullList == "LOW"] <- NA
fullList <- apply(fullList, 2, as.numeric)
fullList <- melt(fullList)
fullList <- fullList[, -1]
colnames(fullList) <- c("Fraction", "Expression")
fullList[complete.cases(fullList), ]

bartlett.test(Expression~Fraction, fullList) # p < 0.0002 | Variances are different

### Figure 2e ###
# Import data #
fullList <- read.csv("fullList.csv", header = TRUE, stringsAsFactors = FALSE)

nrow(subset(fullList, fullList$PSD == "HIGH", select = c("Gene", "PSD")))
nrow(subset(fullList, fullList$PSD == "LOW", select = c("Gene", "PSD")))

# However, note that for the out-of-range genes we used an even more stringent #
# classification. If either treatment had a 2 zero values we would not classify#
# it as out-of-range because we could not compute the fold-change. This will #
# shorten the list. #

### EOF ###
