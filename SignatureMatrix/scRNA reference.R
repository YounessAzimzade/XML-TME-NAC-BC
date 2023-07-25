# Sampling sc IDs data from metadata to create signature matrix (SM)
# Youness Azimzade
# email: youness.azimzade@gmail.com
# 25/07/2023


MetaData<-read.csv("~/metadata.csv")


# removing cell types that have a very small fraction or have a finer resolution 
MetaData <- MetaData[MetaData$celltype_minor!="Cycling_Myeloid" &
                        MetaData$celltype_minor!="Cycling T-cells"&
                        MetaData$celltype_minor!="Cycling PVL",]


# Creating new annotation that combine gene module annotation with  celltype_minor annotation: 
MetaData$Cell_Type <- ifelse(MetaData$gene_module != "no_gene_module", MetaData$gene_module,
                             MetaData$celltype_minor)

 
MetaData$Cell_Type[grepl("^\\d+$", MetaData$Cell_Type)] <- paste0( "GenMod",MetaData$Cell_Type[grepl("^\\d+$", MetaData$Cell_Type)])

# selecting cell ID and cell type
MetaData <-MetaData[,c(1,11)]



MySample  <- MetaData[sample(nrow(MetaData)),]

SM_cellID <- MySample[1:10000,]
PS_cellID <- MySample[10000:nrow(MySample),]
write.table(SM_cellID,file="~/SM_cellID.txt"
            ,sep="\t",row.names=FALSE,quote=FALSE)
 
write.table(PS_cellID,file="~/PS_cellID.txt"
            ,sep="\t",row.names=FALSE,quote=FALSE)

###############################################################
# Retrieving  scRNA seq data for SM_cellID.txt (10000 cells sampled before) 
# from matrix.mtx to create scRNA seq reference matrix  
# Youness Azimzade
# email: youness.azimzade@gmail.com
# 25/07/2023

library(Matrix)

x<-readMM("~/matrix.mtx")
MetaData<-read.csv("~/metadata.csv")
colnames(x)<- MetaData$X

SM_cellID<-read.delim("~/SM_cellID.txt")
Data1<-x[,c(SM_cellID$X)] 
Data1 <- as.matrix(Data1)
xg<-read.table("~/genes.tsv")
CellType<- as.matrix(t(SM_cellID$Cell_Type))
rownames(Data1) <- xg$V1
Data1 <- rbind(CellType, Data1)
write.table(Data1,file="~/scRNAReference.txt"
            ,sep="\t",row.names = TRUE, col.names = FALSE, quote = FALSE)
 