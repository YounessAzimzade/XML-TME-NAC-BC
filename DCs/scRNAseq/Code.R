library(dplyr) 
library(ggplot2)
library(ggpubr) 

AllPat <-read.csv("E:/Dropbox (UiO)/Validation/SI/DCs/scRNAseq/DCsSubsetsinER.csv")

  
subSteFreq <- AllPat%>% count(DCsubset, ERStatus, Patient, sort = FALSE)
colnames(subSteFreq) <- c("DCSubset", "ERStatus","Patient", "Frequency")

pdf("E:/Dropbox (UiO)/Validation/SI/DCs/scRNAseq/FigS12a.pdf",width=9,height=4)
ggplot( subSteFreq,aes(x=`ERStatus`,y=Frequency,col=`ERStatus`))+
  geom_violin()+ geom_boxplot(width=0.1,alpha=0.2) +
  facet_wrap(~DCSubset,ncol=4)+  
  scale_color_manual(values=c("blue","red"))+ 
  stat_compare_means(label="p",size=6,hide.ns=FALSE,label.y=50)+ 
  theme_bw()+theme(axis.text.x=element_text(angle=90))+
  theme(text=element_text(size=20))+
  theme(legend.position="none")
dev.off() 




### Normalizing  
MyAv<- subSteFreq %>%
  group_by(Patient) %>%
  summarize(sum_value = sum(Frequency))
 
subSteFreq$Ave  <- MyAv$sum_value[match(subSteFreq$Patient, MyAv$Patient)] 
subSteFreq$Fraction <- subSteFreq$Frequency/subSteFreq$Ave


pdf("E:/Dropbox (UiO)/Validation/SI/DCs/scRNAseq/FigS12b.pdf",width=9,height=4)
ggplot( subSteFreq,aes(x=`ERStatus`,y=Fraction,col=`ERStatus`))+
  geom_violin()+ geom_boxplot(width=0.1,alpha=0.2) +
  facet_wrap(~DCSubset,ncol=4)+  
  scale_color_manual(values=c("blue","red"))+ 
  stat_compare_means(label="p",size=6,hide.ns=FALSE,label.y=0.8)+ 
  theme_bw()+theme(axis.text.x=element_text(angle=90))+
  theme(text=element_text(size=20))+
  theme(legend.position="none")
dev.off() 
 