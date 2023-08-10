# Figs 5f S15a,b
# Author: Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 10/08/2023

library(ggplot2)
library(reshape2)  
library(ggpubr)     
library(rstatix)

CorAll<-read.csv(".../SpatialCorrelation.csv")

CorAll$ERStatus<-as.character(CorAll$ERStatus)
CorAll$ERStatus<-factor(CorAll$ERStatus,levels= c("Positive","Negative")) 
 
stat.test <- CorAll%>%
  group_by(`ECs Protein`,Protein) %>%
  t_test(Correlation~ERStatus) %>%
  adjust_pvalue(method="bonferroni") %>%
  add_significance("p.adj")
stat.test<-stat.test%>%
  add_xy_position(fun="mean_se",x="ERStatus") 

stat.test2<-stat.test

stat.test2$statistic2<-stat.test2$statistic/abs(stat.test2$statistic)
stat.test2$`-Log(p.adj)*Direction`<- -log10(stat.test2$p.adj)*stat.test2$statistic2
stat.test2<-stat.test2[order(-stat.test2$`-Log(p.adj)*Direction`),]

stat.CD11c <- stat.test2[,c(1,2,18)] 
colnames(stat.CD11c) <- c("Protein", "ECs Protein" ,"-Log(p.adj)*Direction" )

stat.CD11c$`EC Protein`<-as.character(stat.CD11c$`EC Protein`)
stat.CD11c$`EC Protein`<-factor(stat.CD11c$`EC Protein`,levels=c("ER", "HER2","CXCL12"))

pdf(".../Fig5f.pdf",width=14,height=6) 
ggplot(stat.CD11c,aes(x=`ECs Protein`,y=`-Log(p.adj)*Direction`,
                       fill=as.factor(sign(`-Log(p.adj)*Direction`))))+
  geom_bar(stat="identity",width=0.6)+ 
    facet_wrap(~`Protein2`,nrow=1)+  
  scale_fill_manual(values=c("yellow2","purple"))+
  theme_bw()+#theme(axis.text.x = element_text(angle = 90)) +
  theme(text=element_text(size=23))+
  ylim(-1,20)+
  geom_segment(x=0.65,y=18,xend=0.65,yend=0,lineend="round", 
               linejoin="round",size=1,arrow=arrow(length=unit(0.15,"inches"),ends="both"),
               colour="black")+
  annotate(geom="text",size=8,x=0.7,y=19.5,label="ER+",color="black")+ 
  annotate(geom="text",size=8,x=0.7,y=-1,label="ER-",color="black")+
  theme(axis.title.x=element_text(margin=margin(t=9)))+
  theme(legend.position="none")
dev.off()



 
CorAll$`ECs Protein`<-as.character(CorAll$`ECs Protein`)
CorAll$`ECs Protein`<-factor(CorAll$`ECs Protein`,levels= c( "ER", "HER2", "CXCL12")) 

stat.test <- subset(CorAll,CorAll$Protein=="CD11c") %>%
  group_by(`ECs Protein`, Protein) %>%
  t_test(Correlation ~ ERStatus) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test <- stat.test %>%
  add_xy_position(fun = "mean_se", x ="ERStatus") 


dev.off()
pdf(".../FigS15a.pdf",width=10,height=4.5) 
ggplot( subset(CorAll,CorAll$Protein=="CD11c"),
        aes(x=ERStatus,y=Correlation,col=ERStatus))+
  geom_violin()+ 
  geom_boxplot(width=0.1, alpha=0.4) +
  facet_wrap(~ `ECs Protein`,ncol =8)+  
  scale_color_manual(values=c("red",   "blue"))+
  #stat_summary(fun=mean,geom="point",shape=23,size=2)+ 
  # stat_compare_means( label = "p",size=8 ,hide.ns=TRUE,label.y=0.28)  +
  stat_pvalue_manual(stat.test,label="p.adj",hide.ns=TRUE,y.position=0.7,size=7 )+
  # scale_y_continuous(expand=expansion(mult=c(0.03,-0.3))) +
  theme_bw()+theme(axis.text.x=element_text(angle=90)) +
  theme(text=element_text(size=23)) +
  theme(legend.position="none")
dev.off() 


stat.test <- subset(CorAll,CorAll$Protein=="CD16") %>%
  group_by(`ECs Protein`, Protein) %>%
  t_test(Correlation ~ ERStatus) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test <- stat.test %>%
  add_xy_position(fun = "mean_se", x ="ERStatus") 


dev.off()
pdf(".../FigS15b.pdf",width=10,height=4.5) 
ggplot( subset(CorAll,CorAll$Protein=="CD16"),
        aes(x=ERStatus,y=Correlation,col=ERStatus))+
  geom_violin()+ 
  geom_boxplot(width=0.1, alpha=0.4) +
  facet_wrap(~ `ECs Protein`,ncol =8)+  
  scale_color_manual(values=c("red",   "blue"))+
  #stat_summary(fun=mean,geom="point",shape=23,size=2)+ 
  # stat_compare_means( label = "p",size=8 ,hide.ns=TRUE,label.y=0.28)  +
  stat_pvalue_manual(stat.test,label="p.adj",hide.ns=TRUE,y.position=0.6,size=7 )+
  # scale_y_continuous(expand=expansion(mult=c(0.03,-0.3))) +
  theme_bw()+theme(axis.text.x=element_text(angle=90)) +
  theme(text=element_text(size=23)) +
  theme(legend.position="none")
dev.off()


 