setwd("E:/Dropbox (UiO)/Validation/0Code/PAM50/Basal/Analysis") 


library(pheatmap) 
library(ggplot2)
library(reshape2)  
library(ggpubr) 
library(wesanderson)

library(dplyr)
library(lme4) 
library(plyr)
library(broom)
 
 

Data1<-read.csv(".../PAM50/Basal/DiscoveryData.csv")
Data2<-read.csv(".../PAM50/Basal/ValidationData")

Shap1<-read.csv(".../PAM50/Basal/DiscoverySHAP.csv")
Shap2<-read.csv(".../PAM50/Basal/ValidationSHAP.csv")

Data1<-Data1[,4:39]
for(i in 1:ncol(Data1)){ 
  Data1[,i]<-Data1[,i]/(quantile(Data1[,i],probs=0.99))}
Data1[Data1>1] <-NA 


Data2<-Data2[,4:39]
for(i in 1:ncol(Data2)){ 
  Data2[,i]<-Data2[,i]/(quantile(Data2[,i],probs=0.99))}
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

fit1=AllPat1%>%group_by(Feature)%>%do(model=lm(`SHAP value` ~ Fraction, data = .)) 
fit1$Coef<-0
fit1$CI<-0 
for(i in 1:nrow(fit1)){ 
  xx <- as.data.frame(confint(fit1[[2]][[i]], level = 0.999) )
  fit1[i,3]<-fit1[[2]][[i]][["coefficients"]][["Fraction"]]
  fit1[i,4]<-xx[2,1]*xx[2,2]
}


fit2=AllPat2%>%group_by(Feature)%>%do(model=lm(`SHAP value`~Fraction,data=.)) 
fit2$Coef<-0
fit2$CI<-0 
for(i in 1:nrow(fit2)) {
  xx <- as.data.frame(confint(fit2[[2]][[i]],level=0.999))
  fit2[i,3]<-fit2[[2]][[i]][["coefficients"]][["Fraction"]] 
  fit2[i,4]<-xx[2,1]*xx[2,2]
}

SXM <- as.data.frame(fit2[,1])

SXM$Coef1<-fit1$Coef [match(SXM[,1],fit1$Feature)]
SXM$CI1<-fit1$CI [match(SXM[,1],fit1$Feature)]

SXM$Coef2<-fit2$Coef [match(SXM[,1],fit2$Feature)]
SXM$CI2<-fit2$CI [match(SXM[,1],fit2$Feature)]


SXM$RS <-SXM$Coef1*SXM$Coef2 

colnames(SXM)<-c("Feature","Coef1","CI1","Coef2","CI2","RS")


SXM<-subset(SXM,SXM$RS>0&SXM$CI1>0&SXM$CI2>0)

Data <- rbind(Data1, Data2)
Shap <- rbind(Shap1, Shap2)

AllPat <- rbind(AllPat1, AllPat2)


AllPat <- AllPat[ AllPat$Feature %in%SXM$Feature,]
  


fit =AllPat%>%group_by(Feature)%>%do(model=lm(`SHAP value`~Fraction,data=.)) 
fit $Slope <- 0  
fit $LI <- 0 
fit $HI <- 0 
for (i in 1:nrow(fit)) {
  xx <- as.data.frame(confint(fit [[2]][[i]],level=0.999))
  fit [i,3]<-fit[[2]][[i]][["coefficients"]][["Fraction"]] 
  fit [i,4]<-xx[2,1]
  fit [i,5]<-xx[2,2]}

fit  <- fit [,-2]

fit$Sign<-sign(fit$Slope)


fit<-fit[order(-fit$Slope),]

fit$Feature<-as.character(fit$Feature)
fit$Feature<-factor(fit$Feature,levels=fit$Feature)


fit$`Association Score`<-fit$Slope/(quantile(abs(fit$Slope),probs=0.95)) 
fit$`Association Score`[fit$`Association Score`>1] <- 1
fit$`Association Score`[fit$`Association Score`<  -1] <- -1 

fit$LI <- fit$LI*(fit$`Association Score`/fit$Slope) 
fit$HI <- fit$HI*(fit$`Association Score`/fit$Slope) 

colnames(fit)<-c("Cell Type","Slope","LI", "HI", "Sign", "Association Score") 


fit<-subset(fit,fit$`Cell Type`!="CAFs" & fit$`Cell Type`!="PVLs" & fit$
               `Cell Type`!="TCells" & fit$`Cell Type`!="Normal.Epi"& fit$`Cell Type`!="BCells"&
               fit$`Cell Type`!="Endothelials"& fit$`Cell Type`!="Myeloids"& fit$`Cell Type`!="Cancer.Cells")

pdf("E:/Dropbox (UiO)/Validation/0Code/PAM50/Basal/Analysis/Fig4c.pdf",width=16,height=10)
ggplot(fit)+
  geom_bar(aes(x=`Cell Type`, y=`Association Score` , fill= as.factor(Sign)),stat="identity",width = 0.8)+ 
  geom_errorbar( aes(x=`Cell Type`, ymin=LI, ymax=HI), width=0.3, colour="black", alpha=0.8, size=0.8)+
  theme_bw()+theme(text=element_text(size=30))+theme(axis.text.x=element_text(angle=90))+
  theme(legend.position="none")+
  scale_fill_manual(values=c("yellow2","purple"))+
  ylim(-1.3,1.3)+geom_segment(x=0.55,y=1.05,xend=0.55,yend=-1.05,lineend="round", 
                              linejoin="round",size=1,arrow=arrow(length=unit(0.15,"inches"),ends="both"),
                              colour="black")+ annotate(geom="text", size =8,x=0.9,y=1.3, label="pCR",color="black")+ 
  annotate(geom="text",size=8,x=0.8, y=-1.3,label="RD",color="black")+
  theme(axis.title.x=element_text(margin=margin(t=9)))
dev.off() 








  

fitBasal <- fit 

fitBasal<-cbind( "Basal", fitBasal)
colnames(fitBasal) <- c("PAM50","Cell Type","Slope","LI", "HI", "Sign", "`Association Score`")
 

 