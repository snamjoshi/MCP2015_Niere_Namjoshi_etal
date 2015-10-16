# JACoP Analysis
Sanjeev V Namjoshi  
Revised for GitHub: October 16, 2015  

**Code originally written on July 27th, 2015** 

This document includes all the code needed to process the output from the JACoP macro (see ImageJMacros.md). This code was created under the assumption that the JACoP macro is going to output PCC, Manders, and thresholded Manders. The input is the log from the JACoP analysis macro. The titles for the columns will need to be altered to suit your needs.


```r
extract <- function(data){
	dat <- read.table(data, header = FALSE, sep = "", fill = TRUE, stringsAsFactors = FALSE)
	dat <- dat[ , "V1"]
	coef <- dat[grep(pattern = "=", dat)]
	coef <- as.numeric(gsub("r=|M1=|M2=", "",coef))
	coef <- split(coef, ceiling(seq_along(coef)/5))
	coef <- do.call(rbind.data.frame, coef)
	names(coef) <- c("r", "M1", "M2", "M1(T)", "M2(T)")
	coef <- cbind(Value = c("Protein/PSD", "Protein/Synapsin", "PSD/Synapsin"), coef)
	
	return(coef)
}
```
