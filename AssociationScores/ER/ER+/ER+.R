library(pheatmap) 
library(ggplot2)
library(reshape2)  
library(ggpubr) 
library(wesanderson)
library(dplyr)
library(lme4) 
library(plyr)
library(broom)
 
 

Data1<-read.csv("~/SHAP values/ER/ER+/DiscoveryData")
Data2<-read.csv("~/SHAP values/ER/ER+/ValidationData.csv")

Shap1<-read.csv("~/SHAP values/ER/ER+/DiscoverySHAP.csv")
Shap2<-read.csv("~/SHAP values/ER/ER+/ValidationSHAP.csv")


Data1<-Data1[,4:39]
for(i in 1:ncol(Data1)){ 
  Data1[,i]<- Data1[,i]/(quantile(Data1[,i],probs=0.99))
}
Data1[Data1>1] <-NA 


Data2<-Data2[,4:39]
for(i in 1:ncol(Data2)){ 
  Data2[,i]<- Data2[,i]/(quantile(Data2[,i],probs=0.99))
}
Data2[Data2>1] <-NA 

Shap1<- Shap1[,2:37]
Shap2<- Shap2[,2:37]

Data1.m <- melt(Data1) 
Data2.m <- melt(Data2) 

Shap1.m <- melt(Shap1) 
Shap2.m <- melt(Shap2) 

AllPat1 <- cbind(Data1.m,Shap1.m$value)
AllPat2 <- cbind(Data2.m,Shap2.m$value)

colnames(AllPat1) <- c("Feature", "Fraction", "SHAP value") 
colnames(AllPat2) <- c("Feature", "Fraction", "SHAP value") 

AllPat1 <- na.omit(AllPat1)
AllPat2 <- na.omit(AllPat2)


fitt1=AllPat1%>%group_by(Feature)%>%do(model=lm(`SHAP value` ~ Fraction, data = .)) 
fitt1$Coef<-0
fitt1$CI<-0 
for(i in 1:nrow(fitt1)){ 
  xx <- as.data.frame(confint(fitt1[[2]][[i]], level = 0.999) )
  fitt1[i,3]<-fitt1[[2]][[i]][["coefficients"]][["Fraction"]]
  fitt1[i,4]<-xx[2,1]*xx[2,2]
}


fitt2=AllPat2%>%group_by(Feature)%>%do(model=lm(`SHAP value`~Fraction,data=.)) 
fitt2$Coef<-0
fitt2$CI<-0 
for(i in 1:nrow(fitt2)) {
  xx <- as.data.frame(confint(fitt2[[2]][[i]],level=0.999))
  fitt2[i,3]<-fitt2[[2]][[i]][["coefficients"]][["Fraction"]] 
  fitt2[i,4]<-xx[2,1]*xx[2,2]
}

SXM <- as.data.frame(fitt2[,1])

SXM$Coef1<-fitt1$Coef [match(SXM[,1],fitt1$Feature)]
SXM$CI1<-fitt1$CI [match(SXM[,1],fitt1$Feature)]

SXM$Coef2<-fitt2$Coef [match(SXM[,1],fitt2$Feature)]
SXM$CI2<-fitt2$CI [match(SXM[,1],fitt2$Feature)]


SXM$RS <-SXM$Coef1*SXM$Coef2 

colnames(SXM)<-c("Feature","Coef1","CI1","Coef2","CI2","RS")


SXM<-subset(SXM,SXM$RS>0&SXM$CI1>0&SXM$CI2>0)

Data <- rbind(Data1, Data2)
Shap <- rbind(Shap1, Shap2)

AllPat <- rbind(AllPat1, AllPat2)

AllPat <- AllPat[ AllPat$Feature %in%SXM$Feature,]
 




fitt =AllPat%>%group_by(Feature)%>%do(model=lm(`SHAP value`~Fraction,data=.)) 
fitt $Slope <- 0  
fitt $LI <- 0 
fitt $HI <- 0 
for (i in 1:nrow(fitt)) {
  xx <- as.data.frame(confint(fitt [[2]][[i]],level=0.999))
  fitt [i,3]<-fitt[[2]][[i]][["coefficients"]][["Fraction"]] 
  fitt [i,4]<-xx[2,1]
  fitt [i,5]<-xx[2,2]}

fitt  <- fitt [,-2]

fitt$Sign<-sign(fitt$Slope)


fitt<-fitt[order(-fitt$Slope),]

fitt$Feature<-as.character(fitt$Feature)
fitt$Feature<-factor(fitt$Feature,levels=fitt$Feature)
colnames(fitt)<-c("Cell Type","Slope","LI", "HI", "Sign") 


fitt<-subset(fitt,fitt$`Cell Type`!="CAFs" & fitt$`Cell Type`!="PVLs" & fitt$
               `Cell Type`!="TCells" & fitt$`Cell Type`!="Normal.Epi"& fitt$`Cell Type`!="BCells"&
               fitt$`Cell Type`!="Endothelials"& fitt$`Cell Type`!="Myeloids"& fitt$`Cell Type`!="Cancer.Cells")

fitt$`Association Score`<-fitt$Slope/(quantile(fitt$Slope,probs=0.95)) 
fitt$`Association Score`[fitt$`Association Score`>1] <- 1
fitt$`Association Score`[fitt$`Association Score`<  -1] <- -1 

fitt$LI <- fitt$LI*(fitt$`Association Score`/fitt$Slope) 
fitt$HI <- fitt$HI*(fitt$`Association Score`/fitt$Slope) 



pdf("E:/Dropbox (UiO)/Validation/0Code/ER/ER+/Analysis/Fig3a.pdf",width=16,height=10)
ggplot(fitt)+
  geom_bar(aes(x=`Cell Type`, y=`Association Score`, fill= as.factor(Sign)),stat="identity",width = 0.8)+ 
  geom_errorbar( aes(x=`Cell Type`, ymin=LI, ymax=HI), width=0.3, colour="black", alpha=0.8, size=0.8)+
  theme_bw()+theme(text=element_text(size=30))+theme(axis.text.x=element_text(angle=90))+
  theme(legend.position="none")+
  scale_fill_manual(values=c("yellow2","purple"))+
  ylim(-1.15,1.15)+geom_segment(x=0.55,y=1.05,xend=0.55,yend=-1.05,lineend="round", 
                                linejoin="round",size=1,arrow=arrow(length=unit(0.15,"inches"),ends="both"),
                                colour="black")+ annotate(geom="text", size =8,x=0.9,y=1.15, label="pCR",color="black")+ 
  annotate(geom="text",size=8,x=0.8, y=-1.15,label="RD",color="black")+
  theme(axis.title.x=element_text(margin=margin(t=9)))
dev.off() 
 

fittERp <- fitt 
fittERp<-cbind( "ER+", fittERp)
colnames(fittERp) <- c("ER","Cell Type","Slope","LI", "HI", "Sign", "`Association Score`")

 
 