# Calculating distance between DCs and ECs in spatial cyCIF
# Author: Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 10/08/2023

# Load necessary libraries
library(Matrix) 
library(proxy) 

# Read data from CSV
AllCells <- read.csv("~/sc_phenotype.csv")

# Filter out rows with 'Other' phenotype
AllCells <- AllCells[AllCells$Phenotype != "Other",]

# Process data for ER-0086 sample
ERn1 <- AllCells[AllCells$Sample == "ER-0086", ]
ERn1 <- ERn1[, -1]
DCs <- ERn1[ERn1$Phenotype == "DCs",]
DCs <- DCs[, -3]
colnames(DCs) <- c("x", "y") 
Epi <- ERn1[ERn1$Phenotype == "ECs",]
Epi <- Epi[, -3]
colnames(Epi) <- c("x", "y")

# Calculate distance matrix between DCs and ECs
distance_matrixEpi <- proxy::dist(DCs[, c("x", "y")], Epi[, c("x", "y")]) 

# Prepare distance data
distance_dfEpi <- data.frame(distance = as.vector(distance_matrixEpi))
distance_dfEpi$distance <- as.integer(distance_dfEpi$distance)

# Calculate cumulative frequency distribution of distances
DistEpi <- as.data.frame(table(distance_dfEpi$distance))

# Normalizing for the number of DCs
DistEpi$Freq <- (DistEpi$Freq) / (nrow(DCs))  

# calculating the cumulative
DistEpi$Freq <- cumsum(DistEpi$Freq)

# Prepare data for plotting
DistO <- as.data.frame(1:500)
colnames(DistO) <- "Radius"
DistO$`ERn86Epi` <- 0 
DistO$`ERn86Epi` <- DistEpi$Freq[match(DistO$`Radius`, DistEpi$Var1)]  

# Similar processing for other samples (ER-0058, ER+0073, ER+3019)...

# Handle missing values and write the normalized distance data to a CSV file
DistO[is.na(DistO)] = 0
write.csv(as.data.frame(DistO), file = "~/CumFreq.csv",
          sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
