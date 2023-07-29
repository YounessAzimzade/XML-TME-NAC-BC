# Load the required libraries
library(genefu)
library(dplyr)

# Read the data from a file (replace ".../ScaledExpression.txt" with the actual file path)
x <- read.delim(".../ScaledExpression.txt")


# Read patient info for the corresponding dataset. 
PatientInfo <-  read.delim(".../PatientInfo") 

# Transpose the data frame
x <- as.data.frame(t(x))

# Use the first row as column names and remove it from the data frame
colnames(x) <- x[1,]
x <- x[-1,]

# Convert the data frame to numeric values
x2 <- as.data.frame(sapply(x, as.numeric))

# Set row names using the previous row names from 'x'
rownames(x2) <- rownames(x)

# Scale the data using z-score normalization
x2 <- scale(x2)

# Create a data frame with column names as "EntrezGene.ID" for annotation
annot <- as.data.frame(colnames(x2))
colnames(annot) <- "EntrezGene.ID"

# Load the 'pam50' funtion of genefu 
data(pam50)
str(pam50)
data(pam50.robust)
str(pam50.robust)

# Perform intrinsic clustering prediction using the 'pam50' gene list, 'x2' expression data, and 'annot' for annotation
# Note: Replace 'mapping' with the appropriate mapping object if needed, otherwise, remove that argument.
PAMPred <- intrinsic.cluster.predict(pam50, x2, annot, do.mapping = FALSE, mapping, do.prediction.strength = FALSE, verbose = FALSE)

# Extract the subtype information from the prediction results
subType <- as.data.frame(PAMPred$subtype)

# Assuming you have a data frame named 'PatientInfo' and another named 'Samples',
# match the 'Mixture' column in 'PatientInfo' with row names in 'subType' and assign corresponding subtype values to 'PatientInfo$PAM50'
PatientInfo$PAM50 <- subType$`PAMPred$subtype`[match(PatientInfo$Mixture, rownames(subType))]
