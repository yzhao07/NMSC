rm(list=ls())
library(ggplot2)
library(vegan)
library(permute)
library(lattice)
library(plyr)
library(getopt)
library(ggfortify)
library(ade4)
library(dplyr)
library(reshape2)
group = read.csv("../Analysis/A01.1-add-meta/non-melanoma-cohort-sample_10_10.csv")
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
color=c("Responder" = "#A94643", "Nonresponder"="#3E4A7B","Stable_Disease"="#79ADD6","Other"="grey")
i="genus"
########
dat=read.csv(paste0("../Analysis/A01.1-add-meta/jhmi-sears-non-melanoma.",i,".txt"),sep="\t",row.names = 1,check.names = F)
dat=dat[rownames(first_sample),]
data=as.data.frame(dat[,10:ncol(dat)])
data=data[,which(colSums(data)>0.01)]
tab.dist<-vegdist(data,method="bray")
data$group = first_sample$clinical.response
set.seed(123)
data.adonis<- adonis2(tab.dist~group,data,permutations = 999)
pr = data.adonis$`Pr(>F)`[1]
pcoa<- dudi.pco(tab.dist, scan = FALSE,nf=10)
#$eig means how much is each component among 10 components 
pcoa_eig <- (pcoa$eig)/ sum(pcoa$eig)
#get 10 components
plotdata <- data.frame({pcoa$li})
plotdata$names <- rownames(plotdata)
#name 10 components
names(plotdata)[1:3] <- c('PCoA1', 'PCoA2','PCoA3')
plotdata$response = as.factor(first_sample$clinical.response)
centroids <- aggregate(cbind(PCoA1,PCoA2)~response, data=plotdata, mean)
plotdata     <- merge(plotdata,centroids,by="response",suffixes=c("",".centroid"))
pcall <- ggplot(plotdata, aes(x=PCoA1,y=PCoA2,fill = response)) +
  geom_point(aes(x=PCoA1,y=PCoA2,fill=response,color=response), size=4) +
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, fill=response,color=response), pch=21, size=2) +
  geom_segment(aes(x=PCoA1.centroid, y=PCoA2.centroid, xend=PCoA1, yend=PCoA2, color=response), alpha=0.8) +
  scale_fill_manual(values=color) +
  scale_colour_manual(values=color) +
  xlab(paste("PCoA Axis 1 (", round(100 * pcoa_eig[1], 2), "%)", sep="")) +
  ylab(paste("PCoA Axis 2 (", round(100 * pcoa_eig[2], 2), "%)", sep="")) +
  ggtitle("All Tumor")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size=14,colour="black"),legend.position = "none",
        axis.text=element_text(size=14,colour="black"),text=element_text(size=14,colour="black"),
        aspect.ratio=1) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
#################################
dat=read.csv(paste0("../analysis/A01.1-add-meta/jhmi-sears-non-melanoma.",i,".txt"),sep="\t",row.names = 1,check.names = F)
first_sample_bcc = first_sample[first_sample$Tumor == "BCC",]
data=as.data.frame(dat[,10:ncol(dat)])
data=data[,which(colSums(data)>0.01)]
data=data[rownames(first_sample_bcc),]
tab.dist<-vegdist(data,method="bray")
data$group = first_sample_bcc$clinical.response
set.seed(123)
data.adonis<- adonis2(tab.dist~group,data,permutations = 999)
pr = data.adonis$`Pr(>F)`[1]

pcoa<- dudi.pco(tab.dist, scan = FALSE,nf=10)
#$eig means how much is each component among 10 components 
pcoa_eig <- (pcoa$eig)/ sum(pcoa$eig)
#get 10 components
plotdata <- data.frame({pcoa$li})
plotdata$names <- rownames(plotdata)
#name 10 components
names(plotdata)[1:3] <- c('PCoA1', 'PCoA2','PCoA3')
plotdata$response = as.factor(first_sample_bcc$clinical.response)
centroids <- aggregate(cbind(PCoA1,PCoA2)~response, data=plotdata, mean)
plotdata     <- merge(plotdata,centroids,by="response",suffixes=c("",".centroid"))
pc1 <- ggplot(plotdata, aes(x=PCoA1,y=PCoA2,fill = response)) +
  geom_point(aes(x=PCoA1,y=PCoA2,fill=response,color=response), size=4) +
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, fill=response,color=response), pch=21, size=2) +
  geom_segment(aes(x=PCoA1.centroid, y=PCoA2.centroid, xend=PCoA1, yend=PCoA2, color=response), alpha=0.8) +
  scale_fill_manual(values=color) +
  scale_colour_manual(values=color) +
  xlab(paste("PCoA Axis 1 (", round(100 * pcoa_eig[1], 2), "%)", sep="")) +
  ylab(paste("PCoA Axis 2 (", round(100 * pcoa_eig[2], 2), "%)", sep="")) +
  ggtitle("BCC")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size=14,colour="black"),legend.position = "none",
        axis.text=element_text(size=14,colour="black"),text=element_text(size=14,colour="black"),
        aspect.ratio=1) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
##########################################
dat=read.csv(paste0("../analysis/A01.1-add-meta/jhmi-sears-non-melanoma.",i,".txt"),sep="\t",row.names = 1,check.names = F)
first_sample_mcc = first_sample[first_sample$Tumor == "MCC",]
data=as.data.frame(dat[,10:ncol(dat)])
data=data[,which(colSums(data)>0.01)]
data=data[rownames(first_sample_mcc),]
tab.dist<-vegdist(data,method="bray")
data$group = first_sample_mcc$clinical.response
set.seed(123)
data.adonis<- adonis2(tab.dist~group,data,permutations = 999)
pr = data.adonis$`Pr(>F)`[1]

pcoa<- dudi.pco(tab.dist, scan = FALSE,nf=10)
#$eig means how much is each component among 10 components 
pcoa_eig <- (pcoa$eig)/ sum(pcoa$eig)
#get 10 components
plotdata <- data.frame({pcoa$li})
plotdata$names <- rownames(plotdata)
#name 10 components
names(plotdata)[1:3] <- c('PCoA1', 'PCoA2','PCoA3')
plotdata$response = as.factor(first_sample_mcc$clinical.response)

centroids <- aggregate(cbind(PCoA1,PCoA2)~response, data=plotdata, mean)
plotdata     <- merge(plotdata,centroids,by="response",suffixes=c("",".centroid"))

pc3 <- ggplot(plotdata, aes(x=PCoA1,y=PCoA2,fill = response)) +
  geom_point(aes(x=PCoA1,y=PCoA2,fill=response,color=response), size=4) +
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, fill=response,color=response), pch=21, size=2) +
  geom_segment(aes(x=PCoA1.centroid, y=PCoA2.centroid, xend=PCoA1, yend=PCoA2, color=response), alpha=0.8) +
  scale_fill_manual(values=color) +
  scale_colour_manual(values=color) +
  xlab(paste("PCoA Axis 1 (", round(100 * pcoa_eig[1], 2), "%)", sep="")) +
  ylab(paste("PCoA Axis 2 (", round(100 * pcoa_eig[2], 2), "%)", sep="")) +
  ggtitle("MCC")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=14,colour="black"),
        axis.text=element_text(size=14,colour="black"),text=element_text(size=14,colour="black"),
        aspect.ratio=1) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
##########################################
dat=read.csv(paste0("../analysis/A01.1-add-meta/jhmi-sears-non-melanoma.",i,".txt"),sep="\t",row.names = 1,check.names = F)
first_sample_scc = first_sample[first_sample$Tumor == "SCC",]
data=as.data.frame(dat[,10:ncol(dat)])
data=data[,which(colSums(data)>0.01)]
data=data[rownames(first_sample_scc),]
tab.dist<-vegdist(data,method="bray")
data$group = first_sample_scc$clinical.response
set.seed(123)
data.adonis<- adonis2(tab.dist~group,data,permutations = 999)
pr = data.adonis$`Pr(>F)`[1]

pcoa<- dudi.pco(tab.dist, scan = FALSE,nf=10)
#$eig means how much is each component among 10 components 
pcoa_eig <- (pcoa$eig)/ sum(pcoa$eig)
#get 10 components
plotdata <- data.frame({pcoa$li})
plotdata$names <- rownames(plotdata)
#name 10 components
names(plotdata)[1:3] <- c('PCoA1', 'PCoA2','PCoA3')
plotdata$response = as.factor(first_sample_scc$clinical.response)

centroids <- aggregate(cbind(PCoA1,PCoA2)~response, data=plotdata, mean)
plotdata     <- merge(plotdata,centroids,by="response",suffixes=c("",".centroid"))
pc2 <- ggplot(plotdata, aes(x=PCoA1,y=PCoA2,fill = response)) +
  geom_point(aes(x=PCoA1,y=PCoA2,fill=response,color=response), size=4) +
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, fill=response,color=response), pch=21, size=2) +
  geom_segment(aes(x=PCoA1.centroid, y=PCoA2.centroid, xend=PCoA1, yend=PCoA2, color=response), alpha=0.8) +
  scale_fill_manual(values=color) +
  scale_colour_manual(values=color) +
  xlab(paste("PCoA Axis 1 (", round(100 * pcoa_eig[1], 2), "%)", sep="")) +
  ylab(paste("PCoA Axis 2 (", round(100 * pcoa_eig[2], 2), "%)", sep="")) +
  ggtitle("CSCC")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=14,colour="black"),
        axis.text=element_text(size=14,colour="black"),text=element_text(size=14,colour="black"),
        aspect.ratio=1) +
  guides(fill = guide_legend(override.aes=list(shape=21)))

supplement = ggarrange(pcall,pc1,pc3,pc2,nrow = 1,ncol =4,widths=c(1,1,1,1.48) )