###   Code for Fig 3c 
library(pheatmap) 
library(ggplot2)
library(reshape2)  
library(ggpubr)
library(wesanderson)
  
fitAll<-read.csv("~/AssociationScoreAll.csv")
fitAll<-cbind("All",fitAll)
colnames(fitAll)<-c("Data","Cell Type","Slope","LI", "HI", "Sign", "Association Score")
fitAll$`Cell Type` <- as.character(fitAll$`Cell Type`)
fitAll$`Cell Type` <- factor(fitAll$`Cell Type`,levels=fitAll$`Cell Type`)

fitER<- read.csv("E:/Dropbox (UiO)/Validation/0Code/ER/Fit2.csv")
fitER <- rbind(fitERp, fitERn)

colnames(fitER) <- c("Data","Cell Type","Slope","LI", "HI", "Sign", "Association Score")

fitAllx <- rbind(fitAll, fitER)

fitAllx$Data <-  as.character(fitAllx$Data)
fitAllx$Data <- factor(fitAllx$Data,levels= c("All", "ER+", "ER-"))


fitAllx$Sign[fitAllx$Sign==1] <- "Association Score > 0"
fitAllx$Sign[fitAllx$Sign==-1] <- "Association Score < 0"


pdf("E:/Dropbox (UiO)/Validation/0Code/ER/Fig3c.pdf",width=25,height=8.2)
ggplot(fitAllx,aes(`Cell Type`, forcats::fct_rev(as.factor(Data)),fill=as.factor(Sign),size=abs(`Association Score`)))+
  geom_point(shape=21,stroke=0) +
  # geom_hline(yintercept = seq(.5, 4.5, 1), size = .9) +
 scale_x_discrete(position = "top") +
  scale_radius(range=c(10,28)) + 
  # scale_fill_viridis(low = "#F8766D", high = "#00BFC4", limits = c(-1, 1)) +
  scale_fill_manual(values=c("yellow2","purple"))+
  theme_bw() +theme(text=element_text(size =30))+
  theme(axis.text.x=element_text(angle = 90))+
  theme(legend.position = "top", 
        # panel.grid.major = element_blank(),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18))+
  guides(size=guide_legend(override.aes=list(fill = NA,color="black",stroke=.1), 
                           label.position="bottom",
                           title.position="right",order=1),
         fill = guide_legend(override.aes = list(size = 25),
                             title = NULL))+
  labs(size="Area = |Association Score|" ,fill="Association Score",x=NULL,y=NULL)
dev.off()  




AllPat<-read.csv("~/AllSamples.csv") 
 
for(i in 3:30){ 
  AllPat[,i]<- AllPat[,i]/(quantile(AllPat[,i],probs=0.9))
}

AllPat.m<-melt(AllPat[,c(3:30,39)])
colnames(AllPat.m)<-c("ER","Cell Type","Fraction")
AllPat.m[AllPat.m=="Positive"] <-  "ER+"
AllPat.m[AllPat.m=="Negative"] <-  "ER-"

AllPat.m$ER <-  as.character(AllPat.m$ER)
AllPat.m$ER <- factor(AllPat.m$ER,levels= c("ER-", "ER+"))

AllPat.m <- subset(AllPat.m, AllPat.m$`Cell Type` %in% fitAllx$`Cell Type` )


CT <- as.data.frame(fitAllx$`Cell Type`) 
CT<-as.data.frame(CT[!duplicated(CT[,1]),])

AllPat.m$`Cell Type` <-  as.character(AllPat.m$`Cell Type`)
AllPat.m$`Cell Type` <- factor(AllPat.m$`Cell Type`,levels= CT[,1])

  
colnames(AllPat.m) <- c("ER", "Cell Type","fraction in samples of each subtype")


####  Lower panel of Fig 3 
dev.off()
pdf("E:/Dropbox (UiO)/Validation/0Code/ER/Fig3c lowe panel.pdf",width=25,height=2.5)   
ggplot(AllPat.m,aes(x=ER,y=`fraction in samples of each subtype`,col=ER))+theme_bw()+
    geom_boxplot(width=0.95,alpha=1)+ #geom_violin(width=0.9,alpha=0.5)+
  facet_wrap(~`Cell Type`,nrow=1)+   #theme_minimal() +
  theme(strip.background=element_blank(),strip.text.x=element_blank())+   
  theme(text=element_text(size=30))+
  scale_y_continuous(expand=expansion(mult=c(0.002,-0.985)))+
  theme(legend.position="none")+coord_flip()+theme(axis.title.y=element_blank(),axis.text.x=element_blank())+
  scale_color_manual(values=c("blue","red"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank(),axis.line=element_line(colour="black"))
dev.off() 