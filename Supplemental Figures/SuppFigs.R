###############################################################################
###############################################################################
## File: SuppFigs.R                                                          ##
## Author: Sanjeev V Namjoshi (snamjoshi87@utexas.edu)                       ##
## Lab/Advisor: Kimberly Raab-Graham | University of Texas at Austin         ##
## Date: 1/29/15 | Revised for GitHub: 8/10/15                               ## 
## Description: R code used to generate supplemental figures S1 and S2 for   ##
##    paper published in Molecular Cellular Proteomics, 2015                 ##
## Data Files: Use "RawData.csv" from Supplemental Table 1A, cells B12:T756  ##
###############################################################################
###############################################################################

# Load packages
library(ggplot2)
library(reshape2)

# Load data from working directory
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
head(raw)

##############################################################################
############## The "prop" calculates variance across proteins ################
##############################################################################
prop <- function(cols){
	x <- raw[, cols]
	x <- x[rowSums(x == 0) <= 1, ]
	x[x == 0] <- NA
	
	x.var <- apply(x, 1, function(r) {
		sd(r, na.rm = TRUE) / sqrt(sum(!is.na(r)))
	})
	
	x.prop <- x.var / rowMeans(x, na.rm = TRUE)
	return(x.prop)
}

se <- function(x) sqrt(var(x)/length(x))

totals <- list(LD = prop(1:3),
							 LR = prop(4:6),
							 PD = prop(7:9),
							 PR = prop(10:12),
							 SD = prop(13:15),
							 SR = prop(16:18))


# Create se function for use later
se <- function(x) sqrt(var(x)/length(x))

# Arrange the data for graphing by calculating the mean and standard deviation 
# and putting them into a data frame.
totals.mean <- sapply(totals, mean)
totals.se <- sapply(totals, se)
totals.bar <- data.frame(Mean = totals.mean, SE = totals.se)


# For the publication, we put these values into GraphPad. The graph can be replicated 
# in R with the following code:
ggplot(totals.bar, aes(x = rownames(totals.bar), y = Mean, fill = rownames(totals.bar), color = rownames(totals.bar))) +
	geom_bar(stat = "identity", color = "black", fill = "grey") +
	geom_errorbar(aes(ymin = Mean + 0, ymax = Mean + SE), width = 0.2, color = "black") +
	xlab("Fraction") +
	ylab("Normalized Mean S.E.M.") +
	coord_cartesian(ylim = c(0,0.25)) +
	theme_bw()

# The following is not included in the paper but may be useful. We can plot the 
# probability density functions for each of the fractions and see how the 
# distrubution of standard errors compares.
totals.long <- melt(totals)

ggplot(totals.long, aes(x = value, color = L1)) +
	geom_density() +
	xlab("Distribution of S.E.M.") +
	ylab("Density") +
	theme_bw()

## Correlations for figures S1 and S2

# Here we use the coefficient of determination to see the correlation of the 
# triplicate measurements among all normalized spectral counts by treatment 
# and fraction. We will still use the raw data loaded above for the analysis 
# (storted in the variable raw). 

# The correlations need to be done in sets of three. We can create the pattern 
# sequence we need for each set with the code displayed below. We will store the 
# pattern key-index in the variable m.
v1 <- c(1:18) 
v2 <- v1+1L 
v3 <- c(0,0,3L) 
v2 <- v2-v3
m <- cbind(v1,v2)

# We apply the lm() function to the raw data and extracted r.squared in sets 
# of three. The key specifies the pattern in which the correlations are to be 
# computed. We assign the correlations to the variable fit. 
fit <- lapply(1:length(raw),function(x) summary(lm(raw[,m[x,1]]~raw[,m[x,2]]))$r.squared)

# Now we just have to make some cosmetic adjustments to make fit look tidy.
data.frame(Fraction = c(rep("Lysate", 6), rep("PSD", 6), rep("Soluble", 6)),
					 Treatment = rep(c("DMSO","RAPA"), 3, each = 3),
					 Replicate = rep(c("1 vs 2", "2 vs 3", "3 vs 1"), 6),
					 R2 = round(unlist(fit), 3))

# Additionally, we can also look at a correlation matrix. This does not appear 
# in the publication but may also be useful.

# Find correlations across the entire data.
corDat <- cor(raw)
corDatMelt <- melt(corDat)


# To make the graph more readable we have scaled the colors in the correlation 
# matrix across the range of R squared values (rather than scaled over 0 to 1). 
# This makes the color shading much more readable visually. First we create a 
# midpoint function. Then we find the midpoint of the correlations and set it 
# to `mdpt`. This is then set as the midpoint of the color gradient in the 
# scale_fill_gradient2() function call.

# Function to calculate midpoint for use in color scale
midpoint <- function(x) {
	((max(x) - min(x))/2) + min(x)
}

# Set midpoint of data to variable mdpt
mdpt <- midpoint(corDatMelt$value)

# Plot correlation matrix
ggplot(corDatMelt, aes(x = Var1, y = Var2, fill = value)) +
	geom_tile() +
	scale_fill_gradient2(low = "red", mid = "orange", high = "white", midpoint = mdpt) +
	theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)) +
	xlab("") +
	ylab("")
