rm(list=ls())
library(ggplot2)
library(plyr)
library(ggpubr)
library(dplyr)
library(reshape2)
library(stringr)
group = read.csv("../Analysis/A01.1-add-meta/non-melanoma-cohort-sample_10_10.csv",check.names = F)
group = group[order(group$days.sample.collect),]
group = group[order(match(group$patient_rename,unique(group$patient_rename))),]
group = group[order(match(group$clinical.response,c("Responder","Nonresponder","Stable_Disease","Other"))),]
group = group[order(match(group$Tumor,unique(group$Tumor))),]

ID = unique(group$patient_rename)
rownames(group)=group$sample_id

first_sample = NULL
for(i in ID){
  tmp = group[group$patient_rename == i,]
  tmp = tmp[tmp$days.sample.collect == min(tmp$days.sample.collect),]
  first_sample = rbind(first_sample,tmp)
}

first_sample = as.data.frame(first_sample[first_sample$days.sample.collect < 180,])


first_sample = first_sample[first_sample$Tumor == "SCC",]
squamous = group[group$Tumor == "SCC",]

squamous$time_group = NA

for (i in 1:nrow(squamous)){
  if (squamous[i,"days.sample.collect"] < 180){
    squamous[i,"time_group"] = "<180 days"
  }
  else if (squamous[i,"days.sample.collect"] >=180 & squamous[i,"days.sample.collect"]< 360){
    squamous[i,"time_group"] = "180-360 days"
  }
  else if (squamous[i,"days.sample.collect"] >= 360){
    squamous[i,"time_group"] = ">=360 days"
  }}

color2 = c("#FC8002","#4995C6","#369F2D") 
alpha_diversity=read.csv("../Analysis/A01.1-add-meta/jhmi-sears-non-melanoma.alpha-diversity.txt",sep="\t",row.names = 1,check.names = F)
res_squamous = squamous[squamous$clinical.response == "Responder",]
alpha_diversity=alpha_diversity[rownames(res_squamous),]
alpha_diversity$Treatment_Duration = factor(res_squamous$time_group,levels=c("<180 days","180-360 days",">=360 days"))
alpha_diversity$ICI_response = factor(res_squamous$clinical.response,c( "Responder","Nonresponder"))

format = theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.text=element_text(size=12,colour="black"),text=element_text(size=12,colour="black"),
               axis.title = element_text(size=12,colour="black"),plot.title = element_text(size = 12),
               #axis.text.x = element_text(angle = 90, hjust = 1),
               axis.ticks.x = element_blank(),legend.position = "none",
               axis.text.x = element_blank())

p1 = ggplot(alpha_diversity, aes(x =Treatment_Duration, y = shannon,fill=Treatment_Duration))+scale_fill_manual(values=color2)+ 
  geom_boxplot(width=0.4,show.legend = FALSE) + 
  geom_point(position = position_jitter(0.1))+
  stat_summary(fun=median,geom='crossbar',size=0.6,width=0.2)+
  stat_compare_means(comparisons =list(c("<180 days","180-360 days"),c("<180 days",">=360 days"),c(">=360 days","180-360 days")),method = 'wilcox.test',size=4)+
  coord_cartesian(ylim = c(min(alpha_diversity$shannon),max(alpha_diversity$shannon)*1.3))+ 
  theme_bw()+ggtitle("Shannon Index")+theme(axis.title.x = element_blank())+format+ylab("")

p2 = ggplot(alpha_diversity, aes(x =Treatment_Duration, y = simpson_reciprocal,fill=Treatment_Duration))+scale_fill_manual(values=color2)+ 
  geom_boxplot(width=0.4,show.legend = FALSE) + 
  geom_point(position = position_jitter(0.1))+
  stat_summary(fun=median,geom='crossbar',size=0.6,width=0.2)+
  stat_compare_means(comparisons =list(c("<180 days","180-360 days"),c("<180 days",">=360 days"),c(">=360 days","180-360 days")),method = 'wilcox.test',size=4)+
  coord_cartesian(ylim = c(min(alpha_diversity$simpson_reciprocal),max(alpha_diversity$simpson_reciprocal)*1.3))+ 
  theme_bw()+ggtitle("Simpson_reciprocal")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+format
#ggsave(p2,file="./squamous/baseline_alpha_diversity_simpson_reciprocal.pdf", width=4, height=3.5)


p3 = ggplot(alpha_diversity, aes(x =Treatment_Duration, y = chao1,fill=Treatment_Duration))+scale_fill_manual(values=color2)+ 
  geom_boxplot(width=0.4,show.legend = FALSE) + 
  geom_point(position = position_jitter(0.1))+
  stat_summary(fun=median,geom='crossbar',size=0.6,width=0.2)+
  stat_compare_means(comparisons =list(c("<180 days","180-360 days"),c("<180 days",">=360 days"),c(">=360 days","180-360 days")),method = 'wilcox.test',size=4)+
  coord_cartesian(ylim = c(min(alpha_diversity$chao1),max(alpha_diversity$chao1)*1.3))+ 
  theme_bw()+ggtitle("Chao1")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+format
#ggsave(p3,file="./squamous/baseline_alpha_diversity_chao1.pdf", width=4, height=3.5)



p4 = ggplot(alpha_diversity, aes(x =Treatment_Duration, y = observed_species,fill=Treatment_Duration))+scale_fill_manual(values=color2)+ 
  geom_boxplot(width=0.4) + 
  geom_point(position = position_jitter(0.1))+
  stat_summary(fun=median,geom='crossbar',size=0.6,width=0.2)+
  stat_compare_means(comparisons =list(c("<180 days","180-360 days"),c("<180 days",">=360 days"),c(">=360 days","180-360 days")),method = 'wilcox.test',size=4)+
  coord_cartesian(ylim = c(min(alpha_diversity$observed_species),max(alpha_diversity$observed_species)*1.3))+ 
  theme_bw()+ggtitle("Observed_species")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(size = 12),
        axis.text=element_text(size=12,colour="black"),text=element_text(size=12,colour="black"),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title =element_text(size=12,colour="black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))

#labels = "A"                                        # Labels of the scatter plot



combine = ggarrange(p1,p2, p3,p4,ncol = 4,widths  = c(1, 1, 1.1,2))