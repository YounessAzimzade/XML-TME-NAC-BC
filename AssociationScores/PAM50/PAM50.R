library(pheatmap) 
library(ggplot2)
library(reshape2)  
library(ggpubr)
library(wesanderson)
   

fitAll<-read.csv(".../AssociationScoreAll.csv")
fitAll<-cbind("All",fitAll)
colnames(fitAll)<-c("Data","Cell Type","Slope","LI", "HI", "Sign", "Association Score")
fitAll$`Cell Type`<-as.character(fitAll$`Cell Type`)
fitAll$`Cell Type`<-factor(fitAll$`Cell Type`,levels=fitAll$`Cell Type`)

 

fitPAM<-read.csv(".../AssociationScorePAM50.csv")
colnames(fitPAM)<-c("Data","Cell Type","Slope","LI","HI","Sign","Association Score")
 
fitAllx<-rbind(fitAll,fitPAM)  #,


fitAllx$Data<- as.character(fitAllx$Data)
fitAllx$Data <- factor(fitAllx$Data,levels=c("All","ER+", "ER-", "LumA", "LumB", "Basal", "Her2"))

fitAllx$Sign[fitAllx$Sign==1] <- "Association Score > 0"
fitAllx$Sign[fitAllx$Sign==-1] <- "Association Score < 0"

pdf(".../Fig4e.pdf",width=25,height=9.7)
ggplot(fitAllx,aes(`Cell Type`, forcats::fct_rev(as.factor(Data)),fill=as.factor(Sign),size=abs(`Association Score`)))+
  geom_point(shape=21,stroke=0)+
  # geom_hline(yintercept = seq(.5, 4.5, 1), size = .9) +
   scale_x_discrete(position = "top") +
  scale_radius(range=c(12,29)) + 
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


AllPat<-read.csv("E:/Dropbox (UiO)/Validation/All2007x.csv") 
AllPat <- AllPat[AllPat$PAM50!="Normal",]
for(i in 3:30){ 
  AllPat[,i]<- AllPat[,i]/(quantile(AllPat[,i],probs=0.9))
}

AllPat.m<-melt(AllPat[,c(3:30,40)])
#AllPat.m$value [AllPat.m$value > 1.0] <-  1 
colnames(AllPat.m)<-c("PAM50","Cell Type","Fraction") 

AllPat.m$PAM50 <- as.character(AllPat.m$PAM50)
AllPat.m$PAM50  <- factor(AllPat.m$PAM50 , levels=c( "Her2", "Basal","LumB", "LumA"))

AllPat.m <- subset(AllPat.m, AllPat.m$`Cell Type` %in% fitAllx$`Cell Type` )

CT <- as.data.frame(fitAllx$`Cell Type`) 
CT<-as.data.frame(CT[!duplicated(CT[,1]),])

AllPat.m$`Cell Type` <-  as.character(AllPat.m$`Cell Type`)
AllPat.m$`Cell Type` <- factor(AllPat.m$`Cell Type`,levels= CT[,1])

colnames(AllPat.m) <- c("PAM50", "Cell Type","fraction in samples of each subtype")


dev.off()
pdf("E:/Dropbox (UiO)/Validation/0Code/PAM50/Fig4e lower pane;.pdf",width=25,height=4)
ggplot(AllPat.m,aes(x=PAM50,y=`fraction in samples of each subtype`,col=PAM50))+theme_bw()+
 geom_boxplot(width=0.95,alpha=1) + # geom_violin()+
  facet_wrap(~`Cell Type`,nrow=1)+   #theme_minimal() +
  scale_y_continuous(expand=expansion(mult=c(0.002,-0.981)))+
  theme(strip.background=element_blank(),strip.text.x=element_blank())+   
  theme(text=element_text(size=30))+ #ylim(-0.035,0.4)+
  theme(legend.position="none")+coord_flip()+ theme(axis.title.y=element_blank(),axis.text.x=element_blank())+
  scale_color_manual(values=c("lightpink","red","cyan","darkblue"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank(),axis.line=element_line(colour="black"))
dev.off() 