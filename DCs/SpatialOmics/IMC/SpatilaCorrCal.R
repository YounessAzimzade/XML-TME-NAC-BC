# Calculating the spatial correlation between CD11c, CD16 and ER, HER2, CXCL12
# Author: Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 10/08/2023
 
library(reshape2)  
require(dplyr) 


DataAll <-read.csv("~/Metabric/Data.csv")

SingleCells<-read.csv("~/METABRIC/SingleCells.csv") 

SingleCells$HER2 <- SingleCells$HER2..3B5.+SingleCells$HER2..D8F12
SingleCells$DNA<-SingleCells$DNA1 +SingleCells$DNA2

SingleCells <-SingleCells%>%relocate(HER2,.before=ER) 
SingleCells <-SingleCells%>%relocate(DNA,.before=ER) 

SingleCells=SingleCells[,-5:-11]
 
SingleCells=SingleCells[,-c(17,24,44,45,48)]

 

SingleCells[,5:41] <- log(SingleCells[,5:41]+1)
SingleCells2 <-  SingleCells[, 5:41]

for(i in 1:ncol(SingleCells2)){ 
  SingleCells2[,i]<-SingleCells2[,i]/(quantile(SingleCells2[,i],probs=0.99))
}


SingleCells2$metabric_id <- SingleCells$metabric_id 

SingleCells2 <-SingleCells2 %>% relocate(metabric_id, .before = Histone.H3 ) 

colnames(SingleCells2)


func <- function(SingleCells2)
{
  return(data.frame(COR = cor(SingleCells2[,2:38] ,  SingleCells2[,2:38] )))
}
 
CorAll <- ddply(SingleCells2, .(metabric_id), func) 

colnames(CorAll) <- colnames(SingleCells2)

Ps <- list(colnames(SingleCells2[,-1]))
Ps2 <- rep(list(Ps), 718)
Ps2 <- as.data.frame(unlist(Ps2))

CorAll$Protein  <- Ps2[,1]
CorAll <-CorAll %>% relocate(Protein, .before = Histone.H3 ) 

CorAll$PAM50  <- 0
CorAll$ERStatus  <- 0

CorAll$metabric_id<-gsub('-','.',CorAll$metabric_id)

CorAll$PAM50<-DataAll$PAM50[match(CorAll$metabric_id,DataAll$Mixture)] 
CorAll$ERStatus<-DataAll$ER[match(CorAll$metabric_id,DataAll$Mixture)] 



CorAll <-CorAll[CorAll$Protein=="CD11c"|CorAll$Protein=="CD16" ,]

CorAll <- CorAll[, c("metabric_id", "Protein","ER","HER2","CXCL12","ERStatus")] 
CorAll.m <- melt(CorAll)
CorAll.m <- na.omit(CorAll.m)

colnames(CorAll.m) <- c("metabric_id", "Protein","ERStatus","ECs Protein", "Correlation")

write.csv(CorAll.m,file="~/SpatialCorrelation.csv"
          ,sep="\t", row.names = FALSE,quote = FALSE) 






