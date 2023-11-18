# script to perform RMA normalization

library(affy)
library(GEOquery)
library(tidyverse)

# get supplementary files
getGEOSuppFiles("GSE56591")

# untar files
untar("GSE56591/GSE56591_RAW.tar", exdir = 'data/')

# reading in .cel files
raw.data <- ReadAffy(celfile.path = "data/")

# performing RMA normalization
normalized.data <- rma(raw.data)

# get expression estimates
normalized.expr <- as.data.frame(exprs(normalized.data))


# map probe IDs to gene symbols
gse <- getGEO("GSE56591", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE56591_series_matrix.txt.gz@featureData@data
# subset
feature.data <- feature.data[,c(1,11)]

normalized.expr <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data, by = 'ID')

#Loading Necessary Libraries:
  
#  affy: For analyzing Affymetrix GeneChip data.
#GEOquery: For downloading and parsing data from the Gene Expression Omnibus (GEO) database.
#tidyverse: A collection of R packages for data manipulation and visualization.
#Downloading Supplementary Files:
  
# You used getGEOSuppFiles("GSE56591") to download the supplementary files associated with the GEO dataset "GSE56591".
#Extracting the .CEL Files:
  
#You extracted the .CEL files from the downloaded tar archive using the untar function.
#Reading .CEL Files:
  
#You read the .CEL files into R using the ReadAffy function from the affy package.
#RMA Normalization:
  
#You performed RMA (Robust Multi-array Average) normalization on the raw data using the rma function. 
#This method corrects for background noise and normalizes the data across arrays.
#Extracting Expression Estimates:
  
#You extracted the expression estimates from the normalized data and converted it into a data frame.
#Mapping Probe IDs to Gene Symbols:
  
#You fetched the feature data from the GEO dataset to get the mapping between probe IDs and gene symbols.
#You then joined the normalized expression data with the feature data to map probe IDs to their corresponding gene symbols.

head(normalized.expr)


#The resulting normalized.expr data frame contains the following columns:

#ID: The probe ID.
#GSM...: The expression values for each sample (in this case, for different donors and conditions).
#Gene Symbol: The gene symbols corresponding to each probe ID.
#The head of the normalized.expr data frame displays the first few rows of the data, 
#showing the probe IDs, their corresponding expression values for each sample, and the associated gene symbols.



