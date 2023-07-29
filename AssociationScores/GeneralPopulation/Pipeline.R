# Pipeline 
# Author: Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 25/07/2023

# Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)


# Function to normalize cell fractions and remove outliers
normalize_data <- function(data) {
  for (i in 1:ncol(data)) {
    data[, i] <- data[, i] / quantile(data[, i], probs = 0.99)
  }
  data[data > 1] <- NA
  return(data)
}

# Function to fit a line to SHAP vs Fraction for each cell type
fit_line <- function(data) {
  result <- data %>% group_by(Feature) %>% do(model = lm(`SHAP value` ~ Fraction, data = .))
  result$Coef <- 0
  result$CI <- 0
  for (i in 1:nrow(result)) {
    xx <- as.data.frame(confint(result[[2]][[i]], level = 0.999))
    result[i, 3] <- result[[2]][[i]][["coefficients"]][["Fraction"]]
    result[i, 4] <- xx[2, 1] * xx[2, 2]
  }
  return(result)
}

# Read data
DiscoveryData <- read.csv("~/DiscoveryData.csv")
ValidationData <- read.csv("~/ValidationData.csv")
DiscoverySHAP <- read.csv("~/DiscoverySHAP.csv")
ValidationSHAP <- read.csv("~/ValidationSHAP.csv")

# Data preprocessing
DiscoveryData <- DiscoveryData[, 4:39]
DiscoveryData <- normalize_data(DiscoveryData)
ValidationData <- ValidationData[, 4:39]
ValidationData <- normalize_data(ValidationData)
DiscoverySHAP <- DiscoverySHAP[, 2:37]
ValidationSHAP <- ValidationSHAP[, 2:37]

# Combine fractions and SHAP values
DiscoveryData.m <- melt(DiscoveryData)
ValidationData.m <- melt(ValidationData)
DiscoverySHAP.m <- melt(DiscoverySHAP)
ValidationSHAP.m <- melt(ValidationSHAP)
AllPat1 <- cbind(DiscoveryData.m, DiscoverySHAP.m$value)
AllPat2 <- cbind(ValidationData.m, ValidationSHAP.m$value)
colnames(AllPat1) <- c("Feature", "Fraction", "SHAP value")
colnames(AllPat2) <- c("Feature", "Fraction", "SHAP value")
AllPat1 <- na.omit(AllPat1)
AllPat2 <- na.omit(AllPat2)

# Fit lines for Discovery and Validation cohorts
fitDiscv <- fit_line(AllPat1)
fitValid <- fit_line(AllPat2)

# Combining two fits
fitCombined <- merge(fitValid, fitDiscv, by = "Feature", suffixes = c("_Valid", "_Discv"))
 
# Select cell types that show a clear association in both cohorts
fitCombined <- subset(fitCombined, Coef_Discv*Coef_Valid > 0 & CI_Valid > 0 & CI_Discv > 0)

# All data put together
AllPat <- rbind(AllPat1, AllPat2)

# Select only the cell types that pass the validation
AllPat <- AllPat[AllPat$Feature %in% fitCombined$Feature, ]

# A new fitting is performed on discovery and validation data
fit =AllPat%>%group_by(Feature)%>%do(model=lm(`SHAP value`~Fraction,data=.)) 
fit $Coef <- 0  
fit $LI <- 0 
fit $HI <- 0 
for (i in 1:nrow(fit)) {
  xx <- as.data.frame(confint(fit [[2]][[i]],level=0.999))
  fit [i,3]<-fit[[2]][[i]][["coefficients"]][["Fraction"]] 
  fit [i,4]<-xx[2,1]
  fit [i,5]<-xx[2,2]}
 
# Preparations to plot Fig2
fit <- fit[, -2]

fit$Sign <- sign(fit$Coef)

fit <- fit[order(-fit$Coef), ]

fit$Feature <- as.character(fit$Feature)
fit$Feature <- factor(fit$Feature, levels = fit$Feature)
colnames(fit) <- c("Cell Type", "Coef", "LI", "HI", "Sign")

# Remove   major cell types.  By keeping them, we get suppl Fig 13 b. 
fit <- subset(fit, `Cell Type` != "CAFs" & `Cell Type` != "PVLs" & `Cell Type` != "TCells" & 
                `Cell Type` != "Normal.Epi" & `Cell Type` != "BCells" & `Cell Type` != "Endothelials" & 
                `Cell Type` != "Myeloids" & `Cell Type` != "Cancer.Cells")

# Normalize and calculate association score  
fit$`Association Score`<-fit$Coef/(quantile(fit$Coef,probs=0.95))  
fit$`Association Score`[fit$`Association Score`>1] <- 1
fit$`Association Score`[fit$`Association Score`<  -1] <- -1 

# Update the CI scales
fit$LI <- fit$LI*(fit$`Association Score`/fit$Coef)  
fit$HI <- fit$HI*(fit$`Association Score`/fit$Coef) 

# Plot Fig2
pdf("~/Fig2.pdf", width = 20, height = 10.5)
ggplot(fit)+
  geom_bar(aes(x=`Cell Type`, y=`Association Score`, fill= as.factor(Sign)),stat="identity",width = 0.8)+ 
  geom_errorbar( aes(x=`Cell Type`, ymin=LI, ymax=HI), width=0.3, colour="black", alpha=0.9, size=1)+
  theme_bw()+theme(text=element_text(size=30))+theme(axis.text.x=element_text(angle=90))+
  theme(legend.position="none")+
  scale_fill_manual(values=c("yellow2","purple"))+
  scale_x_discrete(position = "bottom") +
  ylim(-1.15,1.15)+geom_segment(x=0.55,y=1.05,xend=0.55,yend=-1.05,lineend="round", 
                                linejoin="round",size=1,arrow=arrow(length=unit(0.15,"inches"),ends="both"),
                                colour="black")+ annotate(geom="text", size =8,x=0.75,y=1.15, label="pCR",color="black")+ 
  annotate(geom="text",size=8,x=0.75, y=-1.13,label="RD",color="black")+
  theme(axis.title.x=element_text(margin=margin(t=9)))
dev.off()               
    

write.csv(fit,file="~/AssociationScoreAll.csv",sep="\t",
          row.names=FALSE,quote=FALSE)              
                 