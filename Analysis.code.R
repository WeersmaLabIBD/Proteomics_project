# this is the whole script for all analysis, including:

# 1) data check, PCA, Umap, clustering 2) multiple variables regression to assess factors
# 3) clinics outcome and host genestics correlation
# 4) CD/UC prediction
# 5) Mendelien Randomization
# 6) eQTL + mbQTL integration
# 7) pQTL affect on cell type enrichment (predicted from RNA-seq)

library(ggplot2)
library(ggridges)
library(dplyr)
library(corrplot)
library(pheatmap)
library(vegan)
library(reshape)
library(foreach)
library(randomcoloR)
library(UpSetR)
library(ggsci)
library(crayon)
library(RColorBrewer)
library(doParallel)
library(scales)
library(glmnet)
library(broom)
library(dplyr)
library(tidyr)
library(stringr)
library(gCMAP)
library(cosinor)
library(ggpubr)
library(cosinor2)
library(reshape2)
library(umapr)
library(tidyverse)
library(uwot)
source("Microbiome.function.R")


# ======================================  data check ======================================
# Proteomics, data import and rename 
pro_all=read.table("data.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
coupling=read.table("Coupling.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
pheno_all=read.table("Olink_proteomics_phenotype_v4.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
surgery_data=read.table("NewSurgeryData!.txt",sep = "\t",header = T,stringsAsFactors = F)
pheno_all=merge(pheno_all,surgery_data,by.x="UMCG1000IBD-ID",by.y="ID",all=T)
aa=pheno_all[,"DateofPlasmaSample",drop=F]
bb=str_split_fixed(aa$`DateofPlasmaSample`, "/", 3)
bb=apply(bb,2,function(x){
  x=as.numeric(x)
})
colnames(bb)=c("year","month","day")
bb=as.data.frame(bb)
bb$day=bb$day+bb$month*30
pheno_all=cbind(bb,pheno_all)

pheno_all$anti_TNF=NA
pheno_all$anti_TNF[pheno_all$IFX_use_all==2]=1
pheno_all$anti_TNF[pheno_all$IFX_use_all==4]=1
pheno_all$anti_TNF[pheno_all$IFX_use_all==5]=1
pheno_all$anti_TNF[pheno_all$IFX_use_all==7]=1
pheno_all$anti_TNF[pheno_all$ADA_use_all==2]=1
pheno_all$anti_TNF[pheno_all$ADA_use_all==4]=1
pheno_all$anti_TNF[pheno_all$ADA_use_all==5]=1
pheno_all$anti_TNF[pheno_all$ADA_use_all==7]=1
pheno_all$anti_TNF[pheno_all$TNF_other_all==2]=1
pheno_all$anti_TNF[is.na(pheno_all$anti_TNF)]=0

pheno_all$Vedolizumab=NA
pheno_all$Vedolizumab[pheno_all$VEDO_use_all==2]=1
pheno_all$Vedolizumab[pheno_all$VEDO_use_all!=2]=0
pheno_all$Vedolizumab[is.na(pheno_all$VEDO_use_all)]=0

pheno_all$Ustekinumab=NA
pheno_all$Ustekinumab[pheno_all$UST_use_all==2]=1
pheno_all$Ustekinumab[pheno_all$UST_use_all==4]=1
pheno_all$Ustekinumab[pheno_all$UST_use_all==5]=1
pheno_all$Ustekinumab[pheno_all$UST_use_all!=2]=0
pheno_all$Ustekinumab[is.na(pheno_all$UST_use_all)]=0

#pheno_all$Season=NA
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-1-" ]=1
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-2-" ]=2
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-3-" ]=3
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-4-" ]=4
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-5-" ]=5
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-6-" ]=6
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-7-" ]=7
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-8-" ]=8
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-9-" ]=9
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-10-" ]=10
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-11-" ]=11
#pheno_all$Season[pheno_all$`Date of Plasma Sample` %like% "-12-" ]=12
#pheno_all$Season=as.factor(pheno_all$Season)

qc_all=read.table("samples.QC.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
pro_all=pro_all[order(pro_all$ProID),]
coupling=coupling[order(coupling$`2D code`),]
qc_all=qc_all[order(qc_all$sample),]
rownames(pro_all)=coupling$`Participant ID`
pro_all$ProID=NULL
rownames(qc_all)=coupling$`Participant ID`

# distribution plot
pro_all_plot=CompositionTable(pro_all,92)
pro_all_plot$Level=factor(pro_all_plot$Level,levels = rev(unique(pro_all_plot$Level)))
ggplot(pro_all_plot, aes(x = Relative, y = Level,fill = Level)) + 
  geom_density_ridges(scale = 3, rel_min_height = 0.01, gradient_lwd = 1.)+
  theme_ridges(font_size = 5, grid = TRUE) + 
  theme(axis.title.y = element_blank())+guides(fill=F)
ggsave("Plot/Pro.distribution.pdf",width = 8,height = 8)
count=data.frame(Pro=unique(pro_all_plot$Level))
count$Valied_number=NA
for(i in 1:nrow(count)){
  pro=as.character(count$Pro)[i]
  tmp=pro_all[,pro]
  tmp=(length(which(!is.na(tmp))))
  count$Valied_number[i]=tmp
}
ggplot(count, aes(x=Pro, y=Valied_number,fill=Pro)) +
  geom_bar(stat='identity') +
  coord_flip()+guides(fill=F)
ggsave("Plot/Barplot.pdf",width = 3,height = 12)

# remove QC failed/warning samples perform PCA, based on manhattan or euclidian distance matrix
pro_all_distance=pro_all
pro_all_distance=pro_all_distance[rownames(pro_all_distance)!="UMCGIBD01528",]
pro_all_distance=pro_all_distance[rownames(pro_all_distance)!="UMCGIBD00980",]
pro_all_distance=pro_all_distance[rownames(pro_all_distance)!="UMCGIBD00727",]
pro_all_distance=pro_all_distance[rownames(pro_all_distance)!="UMCGIBD00092",]
pro_all_distance[is.na(pro_all_distance)]=0
beta_diversity=vegdist(pro_all_distance,method = "manhattan")
PCoAList=do_PCoA(beta_diversity)
pcoas=PCoAList$Coordinates
pcoas=merge(pcoas,pheno_all,by.x="Sample",by.y="UMCG1000IBD-ID",all=T)
pcoas=merge(pcoas,qc_all,by.x="Sample","row.names")
pcoas$qc_warning[pcoas$qc_warning1=="TRUE"]="QC.warned"
pcoas$qc_warning[pcoas$qc_warning2=="TRUE"]="QC.warned"
pcoas$qc_warning[is.na(pcoas$qc_warning)]="QC.passed"
ggplot (pcoas, aes(PCoA1,PCoA2,color=qc_warning)) + 
  geom_point() + theme_bw() + 
  xlab(label = paste("PCoA1",PCoAList$Variance1,sep = " ")) + ylab(label = paste("PCoA2",PCoAList$Variance2,sep = " ")) + 
  theme_bw()+scale_color_jama()
ggsave("./Plot/QC.pdf")

# check per protein
limit=read.table("proteins.txt",sep = "\t",header = T,stringsAsFactors = F)

summary_protein=matrix(nrow = nrow(limit),ncol = 4)
summary_protein=as.data.frame(summary_protein)
colnames(summary_protein)=c("protein","LOD","sample_above_limit","sample_size")
for(i in 1:nrow(limit)){
  tmp.protein=limit$protein[i]
  summary_protein$protein[i]=tmp.protein
  tmp.data=pro_all_distance[,tmp.protein,drop=F]
  tmp.data=tmp.data[tmp.data>limit$LOD[i]]
  summary_protein$sample_above_limit[i]=length(tmp.data)
  summary_protein$sample_size[i]=nrow(pro_all_distance)
  summary_protein$LOD=limit$LOD[i]
}
summary_protein$rate=summary_protein$sample_above_limit/summary_protein$sample_size
write.table(summary_protein,"tables/QC.perProtein.txt",sep = "\t",row.names = F,quote = F)

# Umap clustering
set.seed(10)
umap= umap(pro_all_distance[,c(1:90,92)], init = "spca")
rownames(umap)=rownames(pro_all_distance)
umap=as.data.frame(umap)
umap=merge(umap,pheno_all,by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
umap=merge(umap,qc_all,by.x="Row.names","row.names")
umap$qc_warning[umap$qc_warning1=="TRUE"]="QC.warned"
umap$qc_warning[umap$qc_warning2=="TRUE"]="QC.warned"
umap$qc_warning[is.na(umap$qc_warning)]="QC.passed"
umap$PSC=as.factor(umap$PSC)
umap$anti_TNF=as.factor(umap$anti_TNF)
umap$Vedolizumab=as.factor(umap$Vedolizumab)
umap$Ustekinumab=as.factor(umap$Ustekinumab)
ggplot (umap, aes(V1,V2,color=qc_warning)) + 
  geom_point() + theme_bw() + 
  theme_bw()+scale_color_jama()
ggplot (umap, aes(V1,V2,color=PSC)) + 
  geom_point() + theme_bw() + 
  theme_bw()+scale_color_jama()

umap=umap[umap$`IBD type` %in% c("CD","UC"),]
ggplot (umap, aes(V1,V2,color=`IBD type`)) + 
  geom_point(size=3) + theme_bw() + 
  scale_color_npg()
ggsave("Umap.IBDtype.pdf")

umap_CD=umap[umap$`IBD type`=="CD",]
umap_UC=umap[umap$`IBD type`=="UC",]
p1=ggplot (umap_CD, aes(V1,V2,color=anti_TNF)) + 
  geom_point(size=3) + theme_bw() + 
  theme_bw()+scale_color_lancet()+ theme(legend.position = 'top')
p2=ggplot (umap_UC, aes(V1,V2,color=anti_TNF)) + 
  geom_point(size=3) + theme_bw() + 
  theme_bw()+scale_color_lancet()+ theme(legend.position = 'top')
ggarrange(p1,p2)
ggsave("Plot/Anti_tnf.pdf",height = 5,width = 12)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(200))
umap_CD=merge(umap_CD,pro_all_distance,by.x="Row.names",by.y="row.names",all=F)
umap_UC=merge(umap_UC,pro_all_distance,by.x="Row.names",by.y="row.names",all=F)
p1=ggplot (umap_CD, aes(V1,V2)) + geom_point(aes(colour =TNF),size=3) + theme_bw()+sc+
  theme(legend.position = 'top')
p2=ggplot (umap_UC, aes(V1,V2)) + geom_point(aes(colour =TNF),size=3) + theme_bw()+sc+
  theme(legend.position = 'top')
ggarrange(p1,p2)
ggsave("Plot/protein.tnf.pdf",height = 5,width = 12)

write.table(umap_CD[,c("Row.names","anti_TNF")],"Umap.CD.cluster2.txt",row.names = F,quote = F,sep = "\t")
write.table(umap_UC[,c("Row.names","anti_TNF")],"Umap.UC.cluster2.txt",row.names = F,quote = F,sep = "\t")
write.table(umap_CD,"umap_CD.txt",sep = "\t",row.names = F,quote = F)
write.table(umap_UC,"umap_UC.txt",sep = "\t",row.names = F,quote = F)

qc_warning=qc_all[qc_all$qc_warning1=="TRUE" | qc_all$qc_warning2=="TRUE",]
pro_all_distance=pro_all_distance[!rownames(pro_all_distance) %in% rownames(qc_warning),]
pro_all_distance=pro_all_distance[,c(1:90,92),drop=F]

#pro_all_distance$group=NA
#pro_all_distance$group[rownames(pro_all_distance)=="UMCGIBD01668"]=1
#pro_all_distance$group[rownames(pro_all_distance)=="UMCGIBD00441"]=1
#pro_all_distance$group[rownames(pro_all_distance)=="UMCGIBD00143"]=1
#pro_all_distance$group[rownames(pro_all_distance)=="UMCGIBD00246"]=1
#pro_all_distance$group[rownames(pro_all_distance)=="UMCGIBD00571"]=1
#pro_all_distance$group[is.na(pro_all_distance$group)]=0
#pro_all_distance$group=as.factor(pro_all_distance$group)
#ggplot(pro_all_distance, aes(x=group, y=TNF)) + 
#  geom_boxplot()
#ggsave("Plot/TNFcheck.pdf")
#pro_all_distance$group=NULL
#write.table(pro_all_distance[,c("TNF","group"),drop=F],file = "TNFcheck.txt",row.names = T,quote = F)


# ======================================  factor assess ======================================
# alignment proteins and covariate
pheno_all_distance=pheno_all
pheno_all_distance=pheno_all_distance[pheno_all_distance$`UMCG1000IBD-ID` %in% rownames(pro_all_distance),]

pro_cd=pro_all_distance[rownames(pro_all_distance) %in% pheno_all_distance$`UMCG1000IBD-ID`[pheno_all_distance$`IBD type`=="CD"],]
pro_uc=pro_all_distance[rownames(pro_all_distance) %in% pheno_all_distance$`UMCG1000IBD-ID`[pheno_all_distance$`IBD type`=="UC"],]
pro_covariate=pheno_all_distance[,c("Age","Gender","Plasma Storage Time","BMI","Current smoking","Aminosalicylates","Thiopurines","Steroids",
                                    "Calcineurin_inhibitors","Methotrexaat","Mycofenolatemofetil","oac","antibiotics",
                                    "anti_TNF","History.of.colectomy","History.of.ileocecal.resection")]
rownames(pro_covariate)=pheno_all_distance$`UMCG1000IBD-ID`
pro_covariate$Gender[pro_covariate$Gender=="Female"]=1
pro_covariate$Gender[pro_covariate$Gender=="Male"]=2
pro_covariate$Gender=as.numeric(pro_covariate$Gender)
pro_covariate$Aminosalicylates[is.na(pro_covariate$Aminosalicylates)]=0
pro_covariate$Thiopurines[is.na(pro_covariate$Thiopurines)]=0
pro_covariate$Steroids[is.na(pro_covariate$Steroids)]=0
pro_covariate$Calcineurin_inhibitors[is.na(pro_covariate$Calcineurin_inhibitors)]=0
pro_covariate$Methotrexaat[is.na(pro_covariate$Methotrexaat)]=0
pro_covariate$Mycofenolatemofetil[is.na(pro_covariate$Mycofenolatemofetil)]=0
pro_covariate$oac[is.na(pro_covariate$oac)]=0
pro_covariate$antibiotics[is.na(pro_covariate$antibiotics)]=0

write.table(pro_covariate,"Phenotype.afterQC.txt",sep = "\t",row.names = T,quote = F)

# try clustering samples based on proteins
heatmap_plot=pro_all_distance
heatmap_plot=cbind(heatmap_plot,pheno_all_distance[,c("Gender","IBD type","IBD Medication","PSC"),drop=F])
heatmap_plot=heatmap_plot[heatmap_plot$`IBD type`=="CD" | heatmap_plot$`IBD type`=="UC",]
cols <- colorRampPalette(brewer.pal(6,name="PuOr"))(12)
brks <- seq(-1,1,length.out=12)  
pheatmap(heatmap_plot[,1:92],cluster_cols = T, cluster_rows = T,scale = "column",
         show_rownames=F, show_colnames=T, 
         breaks=brks, border_color=F, color=cols)

# CD and UC compare without any covariate, just cohort compare
cohort_compare=matrix(nrow = ncol(pro_all_distance),ncol = 7)
cohort_compare=as.data.frame(cohort_compare)
colnames(cohort_compare)=c("Protein","CD.mean","CD.sd","UC.mean","UC.sd","Log(fold)","Pvalue")
for(i in 1:ncol(pro_all_distance)){
  
  pro=colnames(pro_all_distance)[i]
  tmp.cd=pro_cd[,i,drop=F]
  tmp.uc=pro_uc[,i,drop=F]
  mm=wilcox.test(tmp.cd[,1],tmp.uc[,1])
  tmp.cd.mean=mean(tmp.cd[,1])
  tmp.uc.mean=mean(tmp.uc[,1])
  tmp.cd.sd=sd(tmp.cd[,1])
  tmp.uc.sd=sd(tmp.uc[,1])
  fold.change=(tmp.cd.mean-tmp.uc.mean)
  cohort_compare$Protein[i]=pro
  cohort_compare$CD.mean[i]=tmp.cd.mean
  cohort_compare$UC.mean[i]=tmp.uc.mean
  cohort_compare$CD.sd[i]=tmp.cd.sd
  cohort_compare$UC.sd[i]=tmp.uc.sd
  cohort_compare$`Log(fold)`[i]=fold.change
  cohort_compare$Pvalue[i]=mm$p.value
  
}
cohort_compare$FDR=p.adjust(cohort_compare$Pvalue)
cohort_compare$threshold=NA
cohort_compare$threshold[cohort_compare$FDR<0.1]="Significant"
cohort_compare$threshold[cohort_compare$FDR>0.1]="Non-Significant"
ggplot(cohort_compare) +
  geom_point(aes(x=`Log(fold)`, y=-log10(Pvalue), fill=threshold),
             shape = 21, color = "black", size = 2) +
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +geom_hline(yintercept=-log10(9.11e-04),linetype="dashed", color = "lightblue")+
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +xlim(-1.4,1.4)+theme_bw()+scale_fill_manual(values = c("grey90","darkred"))
ggsave("Plot/Volcano.pdf")
write.table(cohort_compare,"tables/CD.UC.proteins.compare.txt",sep = "\t",quote = F,row.names = F)

# multivariable model
pro_covariate1=pheno_all_distance[,c("IBD type","Age","Gender","Plasma Storage Time","BMI","Current smoking","Aminosalicylates","Thiopurines","Steroids",
                                    "Calcineurin_inhibitors","Methotrexaat","Mycofenolatemofetil","oac","antibiotics",
                                    "anti_TNF","History.of.colectomy","History.of.ileocecal.resection"),drop=F]
rownames(pro_covariate1)=pheno_all_distance$`UMCG1000IBD-ID`
pro_covariate1$Gender[pro_covariate1$Gender=="Female"]=1
pro_covariate1$Gender[pro_covariate1$Gender=="Male"]=2
pro_covariate1$Gender=as.numeric(pro_covariate1$Gender)
pro_covariate1$Aminosalicylates[is.na(pro_covariate1$Aminosalicylates)]=0
pro_covariate1$Thiopurines[is.na(pro_covariate1$Thiopurines)]=0
pro_covariate1$Steroids[is.na(pro_covariate1$Steroids)]=0
pro_covariate1$Calcineurin_inhibitors[is.na(pro_covariate1$Calcineurin_inhibitors)]=0
pro_covariate1$Methotrexaat[is.na(pro_covariate1$Methotrexaat)]=0
pro_covariate1$Mycofenolatemofetil[is.na(pro_covariate1$Mycofenolatemofetil)]=0
pro_covariate1$oac[is.na(pro_covariate1$oac)]=0
pro_covariate1$antibiotics[is.na(pro_covariate1$antibiotics)]=0
pro_covariate1=pro_covariate1[pro_covariate1$`IBD type`!="DR",]
pro_covariate1=pro_covariate1[pro_covariate1$`IBD type`!="MC",]
pro_covariate1=pro_covariate1[pro_covariate1$`IBD type`!="IBDU",]
pro_covariate1$`IBD type`[pro_covariate1$`IBD type`=="CD"]=0
pro_covariate1$`IBD type`[pro_covariate1$`IBD type`=="UC"]=1
pro_covariate1$`IBD type`=as.numeric(pro_covariate1$`IBD type`)

glm_combine=foreach(i=1:ncol(pro_all_distance),.combine = rbind) %do%  {
  pro=colnames(pro_all_distance)[i]
  tmp.pro=pro_all_distance[,i,drop=F]
  tmp.cov=pro_covariate1[rownames(pro_covariate1) %in% rownames(tmp.pro),]
  tmp.pro=tmp.pro[rownames(tmp.pro) %in% rownames(pro_covariate1),,drop=F]
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  
  tmp.glm=glm(tmp.pro[,1]~.,data=tmp.cov,family = gaussian)
  tmp.coef=as.data.frame(summary(tmp.glm)$coef)
  tmp.coef$FDR=p.adjust(tmp.coef$`Pr(>|t|)`)
  tmp.coef$CIup=tmp.coef$Estimate+1.96*tmp.coef$`Std. Error`
  tmp.coef$CIdown=tmp.coef$Estimate-1.96*tmp.coef$`Std. Error`
  tmp.av=anova(tmp.glm)
  tmp.av$Explain=NA
  for(j in 2:nrow(tmp.av)){
    tmp.av$Explain[j]=(tmp.av$`Resid. Dev`[j-1]-tmp.av$`Resid. Dev`[j])/tmp.av$`Resid. Dev`[1]
  }
  tmp.coef=tmp.coef[tmp.coef$FDR<0.1,]
  tmp.av=merge(tmp.av,tmp.coef,by="row.names",all=F)
  
  if(nrow(tmp.av)==0){
    return.string=data.frame(MultiVariate=NA,MultiVariate.FDR=NA,MultiVariate.explain=NA,Protein=pro,CIup=NA,CIdonw=NA,Estimate=NA)
  }else{
    return.string=data.frame(MultiVariate=tmp.av$Row.names,MultiVariate.FDR=tmp.av$FDR,
                             MultiVariate.explain=tmp.av$Explain,Protein=pro,CIup=tmp.av$CIup,CIdonw=tmp.av$CIdown,Estimate=tmp.av$Estimate)
  }
  
}

glm_cd = foreach(i=1:ncol(pro_cd),.combine = rbind) %do%  {
  pro=colnames(pro_cd)[i]
  tmp.pro=pro_cd[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
    
    tmp.glm=glm(tmp.pro[,1]~.,data=tmp.cov,family = gaussian)
    tmp.coef=as.data.frame(summary(tmp.glm)$coef)
    tmp.coef$FDR=p.adjust(tmp.coef$`Pr(>|t|)`)
    tmp.coef$CIup=tmp.coef$Estimate+1.96*tmp.coef$`Std. Error`
    tmp.coef$CIdown=tmp.coef$Estimate-1.96*tmp.coef$`Std. Error`
    tmp.av=anova(tmp.glm)
    tmp.av$Explain=NA
    for(j in 2:nrow(tmp.av)){
      tmp.av$Explain[j]=(tmp.av$`Resid. Dev`[j-1]-tmp.av$`Resid. Dev`[j])/tmp.av$`Resid. Dev`[1]
    }
    tmp.coef=tmp.coef[tmp.coef$FDR<0.1,]
    tmp.av=merge(tmp.av,tmp.coef,by="row.names",all=F)
    
    if(nrow(tmp.av)==0){
      return.string=data.frame(MultiVariate=NA,MultiVariate.FDR=NA,MultiVariate.explain=NA,Protein=pro,CIup=NA,CIdonw=NA,Estimate=NA)
    }else{
      return.string=data.frame(MultiVariate=tmp.av$Row.names,MultiVariate.FDR=tmp.av$FDR,
                               MultiVariate.explain=tmp.av$Explain,Protein=pro,CIup=tmp.av$CIup,CIdonw=tmp.av$CIdown,Estimate=tmp.av$Estimate)
    }

}
glm_uc = foreach(i=1:ncol(pro_uc),.combine = rbind) %do%  {
  pro=colnames(pro_uc)[i]
  tmp.pro=pro_uc[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  
  tmp.glm=glm(tmp.pro[,1]~.,data=tmp.cov,family = gaussian)
  tmp.coef=as.data.frame(summary(tmp.glm)$coef)
  tmp.coef$FDR=p.adjust(tmp.coef$`Pr(>|t|)`)
  tmp.coef$CIup=tmp.coef$Estimate+1.96*tmp.coef$`Std. Error`
  tmp.coef$CIdown=tmp.coef$Estimate-1.96*tmp.coef$`Std. Error`
  tmp.av=anova(tmp.glm)
  tmp.av$Explain=NA
  for(j in 2:nrow(tmp.av)){
    tmp.av$Explain[j]=(tmp.av$`Resid. Dev`[j-1]-tmp.av$`Resid. Dev`[j])/tmp.av$`Resid. Dev`[1]
  }
  tmp.coef=tmp.coef[tmp.coef$FDR<0.1,]
  tmp.av=merge(tmp.av,tmp.coef,by="row.names",all=F)
  
  if(nrow(tmp.av)==0){
    return.string=data.frame(MultiVariate=NA,MultiVariate.FDR=NA,MultiVariate.explain=NA,Protein=pro,CIup=NA,CIdonw=NA,Estimate=NA)
  }else{
    return.string=data.frame(MultiVariate=tmp.av$Row.names,MultiVariate.FDR=tmp.av$FDR,
                             MultiVariate.explain=tmp.av$Explain,Protein=pro,CIup=tmp.av$CIup,CIdonw=tmp.av$CIdown,Estimate=tmp.av$Estimate)
  }
  
}
glm_cd_sig=glm_cd[glm_cd$MultiVariate.FDR<0.1,]
glm_uc_sig=glm_uc[glm_uc$MultiVariate.FDR<0.1,]
glm_cd_sig=na.omit(glm_cd_sig)
glm_uc_sig=na.omit(glm_uc_sig)

glm_cd_sig_pro=as.character(unique(glm_cd_sig$Protein))
glm_uc_sig_pro=as.character(unique(glm_uc_sig$Protein))
glm_cd_sig_fact=data.frame(Protein=NA,Group="CD",SignificantFactor=NA)
for(i in glm_cd_sig_pro){
  tmp=glm_cd_sig[glm_cd_sig$Protein==i,,drop=F]
  tmp.facotr=paste(tmp$MultiVariate,collapse  = "&")
  tmp.data=data.frame(Protein=i,Group="CD",SignificantFactor=tmp.facotr)
  glm_cd_sig_fact=rbind(glm_cd_sig_fact,tmp.data)
}
glm_cd_sig_fact=na.omit(glm_cd_sig_fact)
glm_uc_sig_fact=data.frame(Protein=NA,Group="UC",SignificantFactor=NA)
for(i in glm_uc_sig_pro){
  tmp=glm_uc_sig[glm_uc_sig$Protein==i,,drop=F]
  tmp.facotr=paste(tmp$MultiVariate,collapse  = "&")
  tmp.data=data.frame(Protein=i,Group="uc",SignificantFactor=tmp.facotr)
  glm_uc_sig_fact=rbind(glm_uc_sig_fact,tmp.data)
}
glm_uc_sig_fact=na.omit(glm_uc_sig_fact)

glm_sig_fact=rbind(glm_cd_sig_fact,glm_uc_sig_fact)
write.table(glm_sig_fact,"tables/Factor.perProtein.txt",sep = "\t",row.names = F,quote = F)

glm_cd$MultiVariate[is.na(glm_cd$MultiVariate)]="Age"
glm_uc$MultiVariate[is.na(glm_uc$MultiVariate)]="Age"
glm_cd[is.na(glm_cd)]=0
glm_uc[is.na(glm_uc)]=0
glm_cd$MultiVariate=gsub("`","",glm_cd$MultiVariate)
glm_uc$MultiVariate=gsub("`","",glm_uc$MultiVariate)
category_order=colnames(pro_covariate)
glm_cd$MultiVariate<-factor(glm_cd$MultiVariate,levels = rev(category_order))
glm_uc$MultiVariate<-factor(glm_uc$MultiVariate,levels = rev(category_order))
sample_order=glm_cd[order(glm_cd$MultiVariate.explain,decreasing = T),]
sample_order=sample_order[sample_order$MultiVariate=="Age",]
sample_order=as.character(sample_order$Protein)
sample_order=append(sample_order,colnames(pro_all_distance)[!colnames(pro_all_distance) %in% sample_order])

assign_color=read.table("Assigin_color.txt",sep = "\t",header = T,stringsAsFactors = F,quote = "",comment.char = "")
palette=assign_color$color
ggplot() +
  geom_bar(data=glm_cd,aes(x=Protein,y=MultiVariate.explain,fill=MultiVariate),stat="identity",colour="black") +
  scale_x_discrete(limits = sample_order)+
  scale_fill_manual(breaks = assign_color$category,values = assign_color$color)+
  geom_bar(data=glm_uc,aes(x=Protein,y=-MultiVariate.explain,fill=MultiVariate), stat="identity",colour="black") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 8),
        axis.line.y =  element_line(colour = 'white'))+
  geom_hline(yintercept = 0,  size=1,colour = 'white')
ggsave("Plot/Multiplevariation.pdf",height = 8,width = 14)

glm_combine$MultiVariate[is.na(glm_combine$MultiVariate)]="Age"
glm_combine[is.na(glm_combine)]=0
glm_combine$MultiVariate=gsub("`","",glm_combine$MultiVariate)
glm_combine$MultiVariate<-factor(glm_combine$MultiVariate,levels = rev(category_order))
sample_order=glm_combine[order(glm_combine$MultiVariate.explain,decreasing = T),]
sample_order=sample_order[sample_order$MultiVariate=="Age",]
sample_order=as.character(sample_order$Protein)
sample_order=append(sample_order,colnames(pro_all_distance)[!colnames(pro_all_distance) %in% sample_order])

ggplot() +
  geom_bar(data=glm_combine,aes(x=Protein,y=MultiVariate.explain,fill=MultiVariate),stat="identity",colour="black") +
  scale_x_discrete(limits = sample_order)+
  scale_fill_manual(breaks = assign_color$category,values = assign_color$color)+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 8),
        axis.line.y =  element_line(colour = 'white'))+
  geom_hline(yintercept = 0,  size=1,colour = 'white')
ggsave("Plot/Multiplevariation.combine.pdf",height = 8,width = 14)

matrix_cd=matrix(ncol = 0,nrow = 91)
matrix_cd=as.data.frame(matrix_cd)
rownames(matrix_cd)=colnames(pro_all_distance)
matrix_cd=matrix_cd[order(rownames(matrix_cd)),,drop=F]
for(i in 1:ncol(pro_covariate)){
  cov=colnames(pro_covariate)[i]
  tmp=data.frame(Protein=colnames(pro_all_distance))
  tmp$R2=NULL
  tmp=merge(tmp,glm_cd[glm_cd$MultiVariate==cov,c("Protein","MultiVariate.explain")],all=T)
  tmp$MultiVariate.explain[tmp$MultiVariate.explain!=0]=1
  tmp[is.na(tmp)]=0
  rownames(tmp)=tmp$Protein
  tmp$Protein=NULL
  colnames(tmp)=cov
  tmp=tmp[order(rownames(tmp)),,drop=F]
  matrix_cd=cbind(matrix_cd,tmp)
}
pdf("./Plot/Upset.plot.CD.pdf",width = 15,height = 6)
upset(matrix_cd, sets = colnames(matrix_cd), sets.bar.color = "#56B4E9",point.size=3,
      order.by = "freq",  keep.order = TRUE)
dev.off()
matrix_uc=matrix(ncol = 0,nrow = 91)
matrix_uc=as.data.frame(matrix_uc)
rownames(matrix_uc)=colnames(pro_all_distance)
matrix_uc=matrix_uc[order(rownames(matrix_uc)),,drop=F]
for(i in 1:ncol(pro_covariate)){
  cov=colnames(pro_covariate)[i]
  tmp=data.frame(Protein=colnames(pro_all_distance))
  tmp$R2=NULL
  tmp=merge(tmp,glm_uc[glm_uc$MultiVariate==cov,c("Protein","MultiVariate.explain")],all=T)
  tmp$MultiVariate.explain[tmp$MultiVariate.explain!=0]=1
  tmp[is.na(tmp)]=0
  rownames(tmp)=tmp$Protein
  tmp$Protein=NULL
  colnames(tmp)=cov
  tmp=tmp[order(rownames(tmp)),,drop=F]
  matrix_uc=cbind(matrix_uc,tmp)
}
pdf("./Plot/Upset.plot.UC.pdf",width = 15,height = 6)
upset(matrix_uc, sets = colnames(matrix_uc), sets.bar.color = "#56B4E9",point.size=3,
      order.by = "freq",  keep.order = TRUE)
dev.off()

write.table(glm_cd,file = "tables/glm.cd.txt",sep = "\t",row.names = F,quote = F)
write.table(glm_uc,file = "tables/glm.uc.txt",sep = "\t",row.names = F,quote = F)

glm_cd_all = foreach(i=1:ncol(pro_cd),.combine = rbind) %do%  {
  pro=colnames(pro_cd)[i]
  tmp.pro=pro_cd[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  
  tmp.glm=glm(tmp.pro[,1]~.,data=tmp.cov,family = gaussian)
  tmp.coef=as.data.frame(summary(tmp.glm)$coef)
  tmp.coef$CIup=tmp.coef$Estimate+1.96*tmp.coef$`Std. Error`
  tmp.coef$CIdown=tmp.coef$Estimate-1.96*tmp.coef$`Std. Error`
  tmp.coef$FDR=p.adjust(tmp.coef$`Pr(>|t|)`)
  tmp.coef=tmp.coef[-1,]
 
  if(nrow(tmp.coef)==0){
    return.string=data.frame(MultiVariate=NA,MultiVariate.P=NA,Protein=pro,CIup=NA,CIdonw=NA,Estimate=NA,FDR=NA)
  }else{
    return.string=data.frame(MultiVariate=rownames(tmp.coef),MultiVariate.P=tmp.coef$`Pr(>|t|)`,
                             Protein=pro,CIup=tmp.coef$CIup,CIdonw=tmp.coef$CIdown,Estimate=tmp.coef$Estimate,FDR=tmp.coef$FDR)
  }
}
glm_uc_all = foreach(i=1:ncol(pro_uc),.combine = rbind) %do%  {
  pro=colnames(pro_uc)[i]
  tmp.pro=pro_uc[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  
  tmp.glm=glm(tmp.pro[,1]~.,data=tmp.cov,family = gaussian)
  tmp.coef=as.data.frame(summary(tmp.glm)$coef)
  tmp.coef$CIup=tmp.coef$Estimate+1.96*tmp.coef$`Std. Error`
  tmp.coef$CIdown=tmp.coef$Estimate-1.96*tmp.coef$`Std. Error`
  tmp.coef$FDR=p.adjust(tmp.coef$`Pr(>|t|)`)
  tmp.coef=tmp.coef[-1,]
  
  if(nrow(tmp.coef)==0){
    return.string=data.frame(MultiVariate=NA,MultiVariate.P=NA,Protein=pro,CIup=NA,CIdonw=NA,Estimate=NA,FDR=NA)
  }else{
    return.string=data.frame(MultiVariate=rownames(tmp.coef),MultiVariate.P=tmp.coef$`Pr(>|t|)`,
                             Protein=pro,CIup=tmp.coef$CIup,CIdonw=tmp.coef$CIdown,Estimate=tmp.coef$Estimate,FDR=tmp.coef$FDR)
  }
}

glm_cd_all$Disease="CD"
glm_uc_all$Disease="UC"
combine=rbind(glm_cd_all,glm_uc_all)

write.table(combine,file = "tables/glm.all.txt",sep = "\t",row.names = F,quote = F)

combine_gender=combine[combine$MultiVariate=="Gender",]
gender=unique(union(glm_cd$Protein[glm_cd$MultiVariate=="Gender"],glm_uc$Protein[glm_uc$MultiVariate=="Gender"]))
combine_gender=combine_gender[combine_gender$Protein %in% gender,]
combine_gender$Protein=as.character(combine_gender$Protein)
combine_gender=combine_gender[order(combine_gender$Estimate,decreasing = T),]
sample_order=as.character(combine_gender[combine_gender$Disease=="CD",]$Protein)
combine_gender$Protein=factor(combine_gender$Protein,levels = rev(sample_order))
combine_gender$Sig=NA
combine_gender$Sig[combine_gender$FDR<0.1]="Significant"
combine_gender$Sig[combine_gender$FDR>0.1]="Significant-non"
ggplot(combine_gender, aes(x=Protein, y=Estimate, ymin=CIdonw, ymax=CIup,col=Disease,fill=Disease,linetype=Sig)) + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_pointrange(aes(col=Disease),position=position_dodge(width = 0.5))+
  coord_flip() +scale_color_jama()+scale_fill_jama()+theme_bw()
ggsave("Plot/gender.effect.proteins.pdf",width = 6,height = 8)

combine_age=combine[combine$MultiVariate=="Age",]
age=unique(union(glm_cd$Protein[glm_cd$MultiVariate=="Age" & glm_cd$MultiVariate.explain!=0],glm_uc$Protein[glm_uc$MultiVariate=="Age" & glm_uc$MultiVariate.explain!=0]))
combine_age=combine_age[combine_age$Protein %in% age,]
combine_age$Protein=as.character(combine_age$Protein)
combine_age=combine_age[order(combine_age$Estimate,decreasing = T),]
sample_order=as.character(combine_age[combine_age$Disease=="CD",]$Protein)
combine_age$Protein=factor(combine_age$Protein,levels = rev(sample_order))
combine_age$Sig=NA
combine_age$Sig[combine_age$FDR<0.1]="Significant"
combine_age$Sig[combine_age$FDR>0.1]="Significant-non"
ggplot(combine_age, aes(x=Protein, y=Estimate, ymin=CIdonw, ymax=CIup,col=Disease,fill=Disease,linetype=Sig)) + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_pointrange(aes(col=Disease),position=position_dodge(width = 0.5))+
  coord_flip() +scale_color_jama()+scale_fill_lancet()+theme_bw()
ggsave("Plot/Age.effect.proteins.pdf",width = 6,height = 8)

combine_Thiopurines=combine[combine$MultiVariate=="Thiopurines",]
Thiopurines=unique(union(glm_cd$Protein[glm_cd$MultiVariate=="Thiopurines" & glm_cd$MultiVariate.explain!=0],glm_uc$Protein[glm_uc$MultiVariate=="Thiopurines" & glm_uc$MultiVariate.explain!=0]))
combine_Thiopurines=combine_Thiopurines[combine_Thiopurines$Protein %in% Thiopurines,]
combine_Thiopurines$Protein=as.character(combine_Thiopurines$Protein)
combine_Thiopurines=combine_Thiopurines[order(combine_Thiopurines$Estimate,decreasing = T),]
sample_order=as.character(combine_Thiopurines[combine_Thiopurines$Disease=="CD",]$Protein)
combine_Thiopurines$Protein=factor(combine_Thiopurines$Protein,levels = rev(sample_order))
combine_Thiopurines$Sig=NA
combine_Thiopurines$Sig[combine_Thiopurines$FDR<0.1]="Significant"
combine_Thiopurines$Sig[combine_Thiopurines$FDR>0.1]="Significant-non"
ggplot(combine_Thiopurines, aes(x=Protein, y=Estimate, ymin=CIdonw, ymax=CIup,col=Disease,fill=Disease,linetype=Sig)) + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_pointrange(aes(col=Disease),position=position_dodge(width = 0.5))+
  coord_flip() +scale_color_jama()+scale_fill_jama()+theme_bw()
ggsave("Plot/Thiopurines.effect.proteins.pdf",width = 6,height = 8)

combine_anti_TNF=combine[combine$MultiVariate=="anti_TNF",]
anti_TNF=unique(union(glm_cd$Protein[glm_cd$MultiVariate=="anti_TNF" & glm_cd$MultiVariate.explain!=0],glm_uc$Protein[glm_uc$MultiVariate=="anti_TNF" & glm_uc$MultiVariate.explain!=0]))
combine_anti_TNF=combine_anti_TNF[combine_anti_TNF$Protein %in% anti_TNF,]
combine_anti_TNF$Protein=as.character(combine_anti_TNF$Protein)
combine_anti_TNF=combine_anti_TNF[order(combine_anti_TNF$Estimate,decreasing = T),]
sample_order=as.character(combine_anti_TNF[combine_anti_TNF$Disease=="CD",]$Protein)
combine_anti_TNF$Protein=factor(combine_anti_TNF$Protein,levels = rev(sample_order))
combine_anti_TNF$Sig=NA
combine_anti_TNF$Sig[combine_anti_TNF$FDR<0.1]="Significant"
combine_anti_TNF$Sig[combine_anti_TNF$FDR>0.1]="Significant-non"
ggplot(combine_anti_TNF, aes(x=Protein, y=Estimate, ymin=CIdonw, ymax=CIup,col=Disease,fill=Disease,linetype=Sig)) + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_pointrange(aes(col=Disease),position=position_dodge(width = 0.5))+
  coord_flip() +scale_color_jama()+scale_fill_jama()+theme_bw()
ggsave("Plot/anti_TNF.effect.proteins.pdf",width = 6,height = 8)

combine_Steroids=combine[combine$MultiVariate=="Steroids",]
Steroids=unique(union(glm_cd$Protein[glm_cd$MultiVariate=="Steroids" & glm_cd$MultiVariate.explain!=0],glm_uc$Protein[glm_uc$MultiVariate=="Steroids" & glm_uc$MultiVariate.explain!=0]))
combine_Steroids=combine_Steroids[combine_Steroids$Protein %in% Steroids,]
combine_Steroids$Protein=as.character(combine_Steroids$Protein)
combine_Steroids=combine_Steroids[order(combine_Steroids$Estimate,decreasing = T),]
sample_order=as.character(combine_Steroids[combine_Steroids$Disease=="CD",]$Protein)
combine_Steroids$Protein=factor(combine_Steroids$Protein,levels = rev(sample_order))
combine_Steroids$Sig=NA
combine_Steroids$Sig[combine_Steroids$FDR<0.1]="Significant"
combine_Steroids$Sig[combine_Steroids$FDR>0.1]="Significant-non"
ggplot(combine_Steroids, aes(x=Protein, y=Estimate, ymin=CIdonw, ymax=CIup,col=Disease,fill=Disease,linetype=Sig)) + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_pointrange(aes(col=Disease),position=position_dodge(width = 0.5))+
  coord_flip() +scale_color_jama()+scale_fill_jama()+theme_bw()
ggsave("Plot/Steroids.effect.proteins.pdf",width = 6,height = 8)

medication=colnames(pro_covariate)[6:14]
medication=c("Aminosalicylates","Thiopurines","Steroids","`TNF-antagonists`","Calcineurin inhibitors",
             "Methotrexate","Mycophenolate mofetil","Ustekinumab","Vedolizumab")
for(i in c("Aminosalicylates","Thiopurines","Steroids","`TNF-antagonists`","Calcineurin inhibitors",
           "Methotrexate","Mycophenolate mofetil","Ustekinumab","Vedolizumab")){
  tmp=glm_cd_all
}
tmp_cd_medication=glm_cd_all[glm_cd_all$MultiVariate %in% medication,c("MultiVariate","Protein","Estimate")]
tmp_cd_medication=dcast(tmp_cd_medication,Protein ~ MultiVariate, value.var="Estimate")
rownames(tmp_cd_medication)=tmp_cd_medication$Protein
tmp_cd_medication$Protein=NULL
clustering=dist((as.matrix(tmp_cd_medication)),method = "euclidean")
clustering=hclust(clustering, method = "complete", members = NULL)
plot(clustering)
#clustering=dist(as.matrix(tmp_cd_medication),method = "euclidean")
#clustering=hclust(clustering, method = "complete", members = NULL)
pdf("Plot/Protein.medication.cluster.pdf",width = 30,height = 30)
plot(clustering,hang = -1, cex = 0.6)
dev.off()

sample_order=rownames(tmp_cd_medication)[clustering$order]
tmp_cd_medication=glm_cd_all[glm_cd_all$MultiVariate %in% medication,c("MultiVariate","Protein","Estimate","FDR")]
tmp_cd_medication$Protein=factor(tmp_cd_medication$Protein,levels = sample_order)
tmp_cd_medication$MultiVariate=factor(tmp_cd_medication$MultiVariate,levels = c("Methotrexate","Aminosalicylates",
                                                                                "Thiopurines","`TNF-antagonists`",
                                                                                "Steroids","Vedolizumab","Ustekinumab"))
tmp_cd_medication$color=NA
tmp_cd_medication$color[tmp_cd_medication$Estimate<0]="Decrease"
tmp_cd_medication$color[tmp_cd_medication$Estimate>0]="Increase"
tmp_cd_medication$color[tmp_cd_medication$FDR>0.05]="NotSig"
ggplot(tmp_cd_medication, aes(x=Protein, y=MultiVariate, fill=color)) + 
  geom_tile(aes(fill=color), colour="white")+
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y= "") + 
  coord_fixed(ratio=1)


# ======================================  clinical outcome/genetics correlation ======================================
# correct for all significant covariates
corrected_cd = foreach(i=1:ncol(pro_cd),.combine = cbind) %do%  {
  pro=colnames(pro_cd)[i]
  cat(green("Correction ===>",pro,"\n"))
  tmp.pro=pro_cd[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.sig.cov=tmp.cov[,colnames(tmp.cov) %in% glm_cd$MultiVariate[glm_cd$Protein==pro & glm_cd$MultiVariate.explain>0],drop=F]

  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.sig.cov=tmp.sig.cov[order(rownames(tmp.sig.cov)),,drop=F]
  
  if(ncol(tmp.sig.cov)==0){
    
    return.string=as.data.frame(tmp.pro)
    
  }else{
    tmp.sig.cov=as.data.frame(apply(tmp.sig.cov,2,function(x){
      x[is.na(x)]=median(x[!is.na(x)])
      return(x)
    }))
    
    x.resid = resid(lm(tmp.pro[,1] ~ .,data = tmp.sig.cov))
    tmp.pro[,1]=x.resid
    return.string=as.data.frame(tmp.pro)
  }
}
corrected_uc = foreach(i=1:ncol(pro_uc),.combine = cbind) %do%  {
  pro=colnames(pro_uc)[i]
  cat(green("Correction ===>",pro,"\n"))
  tmp.pro=pro_uc[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.sig.cov=tmp.cov[,colnames(tmp.cov) %in% glm_uc$MultiVariate[glm_uc$Protein==pro & glm_uc$MultiVariate.explain>0],drop=F]
  
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.sig.cov=tmp.sig.cov[order(rownames(tmp.sig.cov)),,drop=F]
  
  if(ncol(tmp.sig.cov)==0){
    
    return.string=as.data.frame(tmp.pro)
    
  }else{
    tmp.sig.cov=as.data.frame(apply(tmp.sig.cov,2,function(x){
      x[is.na(x)]=median(x[!is.na(x)])
      return(x)
    }))
    
    x.resid = resid(lm(tmp.pro[,1] ~ .,data = tmp.sig.cov))
    tmp.pro[,1]=x.resid
    return.string=as.data.frame(tmp.pro)
  }
}
corrected_cd = as.data.frame(t(corrected_cd))
corrected_cd = cbind(rownames(corrected_cd),corrected_cd)
colnames(corrected_cd)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(corrected_cd),
                   Gene = rownames(corrected_cd),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"

write.table(corrected_cd, file = "CorrectedData/CD_numeric.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "CorrectedData/CD_numeric.txt.annot",sep="\t",row.names=F,quote = F)
corrected_uc = as.data.frame(t(corrected_uc))
corrected_uc = cbind(rownames(corrected_uc),corrected_uc)
colnames(corrected_uc)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(corrected_uc),
                   Gene = rownames(corrected_uc),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"

write.table(corrected_uc, file = "CorrectedData/UC_numeric.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "CorrectedData/UC_numeric.txt.annot",sep="\t",row.names=F,quote = F)

#write.table(rownames(pro_cd),file = "CD.list.txt",row.names = F,quote = F)
#write.table(rownames(pro_uc),file = "UC.list.txt",row.names = F,quote = F)

# HBI, SCCAI, CRP,  montreal classification 
corrected_cd=as.data.frame(t(corrected_cd))
corrected_uc=as.data.frame(t(corrected_uc))
corrected_cd=corrected_cd[-1,]
corrected_uc=corrected_uc[-1,]

hbi_cd = foreach(i=1:ncol(corrected_cd),.combine = rbind) %do%  {
  
  pro=colnames(corrected_cd)[i]
  tmp.pro=corrected_cd[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","HBI"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method = "spearman")
  
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,R=tmp$estimate)
}
hbi_cd$FDR=p.adjust(hbi_cd$Pvalue)

tmp=pro_all_distance[,colnames(pro_all_distance) %in% as.character(hbi_cd$Protein[hbi_cd$FDR<0.1]),drop=F]
tmp=merge(tmp,pheno_all_distance[pheno_all_distance$`IBD type`=="CD",c("UMCG1000IBD-ID","HBI")],by.x = "row.names",by.y = "UMCG1000IBD-ID")
rownames(tmp)=tmp$Row.names
tmp$Row.names=NULL
tmp_plot=data.frame(Value=NA,Protein=NA,HBI=NA)
for(i in 1:7){
  tmp.data=data.frame(row.names = rownames(tmp),Value=tmp[,i],Protein=colnames(tmp)[i],HBI=tmp$HBI)
  tmp_plot=rbind(tmp_plot,tmp.data)
}
tmp_plot=tmp_plot[-1,]
tmp_plot=na.omit(tmp_plot)
tmp_plot$HBI=as.numeric(tmp_plot$HBI)
ggplot(tmp_plot, aes(x = HBI, y = Value))+
  geom_smooth(aes(color = Protein, fill = Protein), method = lm) +
  scale_color_manual(values = pal_rickandmorty("schwifty")(11))+
  scale_fill_manual(values = pal_rickandmorty("schwifty")(11))+ylab("")+theme_bw()
ggsave("Plot/HBI.elevenGenes.pdf",width = 6,height = 3)

hbi_quartile=pheno_all_distance[pheno_all_distance$`IBD type`=="CD",c("UMCG1000IBD-ID","HBI")]
hbi_quartile=na.omit(hbi_quartile)
quantile(hbi_quartile$HBI)
hbi_quartile$Group=NA
hbi_quartile$Group[hbi_quartile$HBI==0]=1
hbi_quartile$Group[hbi_quartile$HBI<=2 & hbi_quartile$HBI>0]=2
hbi_quartile$Group[hbi_quartile$HBI<=5 & hbi_quartile$HBI>2]=3
hbi_quartile$Group[hbi_quartile$HBI<=30 & hbi_quartile$HBI>5]=4

hbi_cd = foreach(i=1:ncol(corrected_cd),.combine = rbind) %do%  {
  
  pro=colnames(corrected_cd)[i]
  tmp.pro=corrected_cd[,i,drop=F]
  tmp.pro=merge(tmp.pro,hbi_quartile,by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp=lm(tmp.pro[,2]~tmp.pro[,4])
  tmp=summary(tmp)$coef
  return.string=data.frame(Protein=pro,Pvalue=tmp[2,4],R=tmp[2,1])
}
hbi_cd$FDR=p.adjust(hbi_cd$Pvalue)

ggplot(tmp_plot, aes(tmp_plot$HBI, tmp_plot$Value,color=Protein,fill=Protein)) +
  geom_point(shape = 21, 
             color = "black", size = 2)+
  geom_smooth(method = lm)+
  facet_wrap( .~ Protein, scales="free_y",nrow = 1)+scale_color_jama()+
  theme_bw()+ theme(legend.position="bottom")+guides(color=guide_legend(nrow = 1))
ggsave("Plot/HBI.pdf",width = 12,height = 4)
write.table(hbi_cd,"tables/HBI.txt",sep = "\t",row.names = F,quote = F)

SCCAI_uc = foreach(i=1:ncol(corrected_uc),.combine = rbind) %do%  {
  
  pro=colnames(corrected_uc)[i]
  tmp.pro=corrected_uc[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","SCCAI"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method = "spearman")
  
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,R=tmp$estimate)
}
SCCAI_uc$FDR=p.adjust(SCCAI_uc$Pvalue)

tmp=pro_all_distance[,colnames(pro_all_distance) %in% as.character(SCCAI_uc$Protein[SCCAI_uc$FDR<0.1]),drop=F]
tmp=merge(tmp,pheno_all_distance[pheno_all_distance$`IBD type`=="UC",c("UMCG1000IBD-ID","SCCAI")],by.x = "row.names",by.y = "UMCG1000IBD-ID")
rownames(tmp)=tmp$Row.names
tmp$Row.names=NULL
tmp_plot=data.frame(Value=NA,Protein=NA,SCCAI=NA)
for(i in 1:7){
  tmp.data=data.frame(row.names = rownames(tmp),Value=tmp[,i],Protein=colnames(tmp)[i],SCCAI=tmp$SCCAI)
  tmp_plot=rbind(tmp_plot,tmp.data)
}
tmp_plot=tmp_plot[-1,]
tmp_plot=na.omit(tmp_plot)
tmp_plot$SCCAI=as.numeric(tmp_plot$SCCAI)
ggplot(tmp_plot, aes(x = SCCAI, y = Value))+
  geom_smooth(aes(color = Protein, fill = Protein), method = lm) +
  scale_color_npg()+
  scale_fill_npg()+ylab("")+theme_bw()
ggsave("Plot/SCCAI.eightGenes.pdf",width = 6,height = 3)

ggplot(tmp_plot, aes(tmp_plot$SCCAI, tmp_plot$Value,color=Protein,fill=Protein)) +
  geom_point(shape = 21, 
             color = "black", size = 2)+
  geom_smooth(method = lm)+
  facet_wrap( .~ Protein, scales="free_y",nrow = 1)+scale_color_jama()+
  theme_bw()+ theme(legend.position="bottom")+guides(color=guide_legend(nrow = 1))
ggsave("Plot/SCCAI.pdf",width = 12,height = 4)
write.table(SCCAI_uc,"tables/SCCAI.txt",sep = "\t",row.names = F,quote = F)

Mon_age_cd = foreach(i=1:ncol(corrected_cd),.combine = rbind) %do%  {
  
  pro=colnames(corrected_cd)[i]
  tmp.pro=corrected_cd[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","Montreal Age"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method = "spearman")
  
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,R=tmp$estimate)
}
Mon_age_cd$FDR=p.adjust(Mon_age_cd$Pvalue)
write.table(Mon_age_cd,file = "tables/Mon_age.CD.txt",row.names = F,quote = F,sep = "\t")

Mon_age_uc = foreach(i=1:ncol(corrected_uc),.combine = rbind) %do%  {
  
  pro=colnames(corrected_uc)[i]
  tmp.pro=corrected_uc[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","Montreal Age"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method = "spearman")
  
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,R=tmp$estimate)
}
Mon_age_uc$FDR=p.adjust(Mon_age_uc$Pvalue)
write.table(Mon_age_uc,file = "tables/Mon_age.UC.txt",row.names = F,quote = F,sep = "\t")

pheno_all_distance$`Montreal Location`[pheno_all_distance$`Montreal Location`==4]=0
pheno_all_distance$`Montreal Location`[pheno_all_distance$`Montreal Location`==5]=1
pheno_all_distance$`Montreal Location`[pheno_all_distance$`Montreal Location`==6]=2
Mon_locate_cd = foreach(i=1:ncol(corrected_cd),.combine = rbind) %do%  {
  
  pro=colnames(corrected_cd)[i]
  tmp.pro=corrected_cd[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","Montreal Location"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp.pro[,3]=as.factor(tmp.pro[,3])
  tmp=kruskal.test(tmp.pro[,2]~tmp.pro[,3])
  
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value)
}
Mon_locate_cd$FDR=p.adjust(Mon_locate_cd$Pvalue)
tmp=pro_all_distance[,colnames(pro_all_distance) %in% as.character(Mon_locate_cd$Protein[Mon_locate_cd$FDR<0.1]),drop=F]
tmp=merge(tmp,pheno_all_distance[pheno_all_distance$`IBD type`=="CD",c("UMCG1000IBD-ID","Montreal Location")],by.x = "row.names",by.y = "UMCG1000IBD-ID")
rownames(tmp)=tmp$Row.names
tmp$Row.names=NULL
tmp_plot=data.frame(Value=NA,Protein=NA,Location=NA)
for(i in 1:1){
  tmp.data=data.frame(row.names = rownames(tmp),Value=tmp[,i],Protein=colnames(tmp)[i],Location=tmp$`Montreal Location`)
  tmp_plot=rbind(tmp_plot,tmp.data)
}
tmp_plot=tmp_plot[-1,]
tmp_plot=na.omit(tmp_plot)
tmp_plot=tmp_plot[tmp_plot$Location!=3,]
tmp_plot$Location=as.factor(tmp_plot$Location)
ggplot(tmp_plot, aes(Location, Value,color=Location,fill=Location)) +
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+theme_bw()+scale_fill_jama()
ggsave("Plot/MontrealLocation.pdf",width = 8,height = 5)
write.table(Mon_locate_cd,file = "tables/Mon_locate.CD.txt",row.names = F,quote = F,sep = "\t")

wilcox.test(tmp$`FGF-19`[tmp$`Montreal Location`=="0"],tmp$`FGF-19`[tmp$`Montreal Location`=="1"])

Mon_behave_cd = foreach(i=1:ncol(corrected_cd),.combine = rbind) %do%  {
  
  pro=colnames(corrected_cd)[i]
  tmp.pro=corrected_cd[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","Montreal Behavior"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp.pro[,3]=as.factor(tmp.pro[,3])
  tmp=kruskal.test(tmp.pro[,2]~tmp.pro[,3])
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value)
}
Mon_behave_cd$FDR=p.adjust(Mon_behave_cd$Pvalue)
tmp=pro_all_distance[,colnames(pro_all_distance) %in% as.character(Mon_behave_cd$Protein[Mon_behave_cd$FDR<0.1]),drop=F]
tmp=merge(tmp,pheno_all_distance[pheno_all_distance$`IBD type`=="CD",c("UMCG1000IBD-ID","Montreal Behavior")],by.x = "row.names",by.y = "UMCG1000IBD-ID")
rownames(tmp)=tmp$Row.names
tmp$Row.names=NULL
tmp_plot=data.frame(Value=NA,Protein=NA,CD_behave=NA)
for(i in 1){
  tmp.data=data.frame(row.names = rownames(tmp),Value=tmp[,i],Protein=colnames(tmp)[i],CD_behave=tmp$`Montreal Behavior`)
  tmp_plot=rbind(tmp_plot,tmp.data)
}
tmp_plot=tmp_plot[-1,]
tmp_plot=na.omit(tmp_plot)
tmp_plot$CD_behave=as.factor(tmp_plot$CD_behave)
ggplot(tmp_plot, aes(CD_behave, Value,color=CD_behave,fill=CD_behave)) +
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+theme_bw()+scale_fill_jama()
ggsave("Plot/CD_behave.oneGenes.pdf",width = 12,height = 5)
write.table(Mon_behave_cd,file = "tables/Mon_behave.CD.txt",row.names = F,quote = F,sep = "\t")

wilcox.test(tmp$`FGF-19`[tmp$`Montreal Behavior`=="0"],tmp$`FGF-19`[tmp$`Montreal Behavior`=="2"])

Mon_exten_uc = foreach(i=1:ncol(corrected_uc),.combine = rbind) %do%  {
  
  pro=colnames(corrected_uc)[i]
  tmp.pro=corrected_uc[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","Montreal Extension"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp.pro[,3]=as.factor(tmp.pro[,3])
  tmp=kruskal.test(tmp.pro[,2]~tmp.pro[,3])
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.valu)
}
Mon_exten_uc$FDR=p.adjust(Mon_exten_uc$Pvalue)
tmp=pro_all_distance[,colnames(pro_all_distance) %in% as.character(Mon_exten_uc$Protein[Mon_exten_uc$FDR<0.1]),drop=F]
tmp=merge(tmp,pheno_all_distance[pheno_all_distance$`IBD type`=="UC",c("UMCG1000IBD-ID","Montreal Extension")],by.x = "row.names",by.y = "UMCG1000IBD-ID")
rownames(tmp)=tmp$Row.names
tmp$Row.names=NULL
tmp_plot=data.frame(Value=NA,Protein=NA,UC_extension=NA)
for(i in 1:8){
  tmp.data=data.frame(row.names = rownames(tmp),Value=tmp[,i],Protein=colnames(tmp)[i],UC_extension=tmp$`Montreal Extension`)
  tmp_plot=rbind(tmp_plot,tmp.data)
}
tmp_plot=tmp_plot[-1,]
tmp_plot=na.omit(tmp_plot)
tmp_plot$UC_extension=as.factor(tmp_plot$UC_extension)
ggplot(tmp_plot, aes(UC_extension, Value,color=UC_extension,fill=UC_extension)) +
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+theme_bw()+scale_fill_jama()+facet_wrap( .~ Protein, scales="free_y",ncol = 2)
ggsave("Plot/UC_extension.pdf",width = 12,height = 12)
write.table(Mon_exten_uc,file = "tables/Mon_exten.UC.txt",row.names = F,quote = F,sep = "\t")

Mon_sever_uc = foreach(i=1:ncol(corrected_uc),.combine = rbind) %do%  {
  
  pro=colnames(corrected_uc)[i]
  tmp.pro=corrected_uc[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","Montreal Severity"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method="spearman")
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,estimate=tmp$estimate)
}
Mon_sever_uc$FDR=p.adjust(Mon_sever_uc$Pvalue)

write.table(Mon_sever_uc,file = "tables/Mon_sever.UC.txt",row.names = F,quote = F,sep = "\t")

tmp=pro_all_distance[,colnames(pro_all_distance) %in% as.character(Mon_sever_uc$Protein[Mon_sever_uc$FDR<0.1]),drop=F]
tmp=merge(tmp,pheno_all_distance[pheno_all_distance$`IBD type`=="UC",c("UMCG1000IBD-ID","Montreal Severity")],by.x = "row.names",by.y = "UMCG1000IBD-ID")
rownames(tmp)=tmp$Row.names
tmp$Row.names=NULL
tmp_plot=data.frame(Value=NA,Protein=NA,UC_severity=NA)
for(i in 1:3){
  tmp.data=data.frame(row.names = rownames(tmp),Value=tmp[,i],Protein=colnames(tmp)[i],UC_severity=tmp$`Montreal Severity`)
  tmp_plot=rbind(tmp_plot,tmp.data)
}
tmp_plot=tmp_plot[-1,]
tmp_plot=na.omit(tmp_plot)
tmp_plot$UC_severity=as.numeric(tmp_plot$UC_severity)
ggplot(tmp_plot, aes(UC_severity, Value,color=Protein,fill=Protein)) +
  geom_jitter(shape = 21, 
             color = "black", size = 2)+
  geom_smooth(method = lm)+
  facet_wrap( .~ Protein, scales="free_y")+scale_color_jama()+
  theme_bw()
ggsave("Plot/UC_severitys.pdf",width = 6,height = 3)

write.table(Mon_behave_cd,"tables/Mon_behave.cd.txt",sep = "\t",row.names = F,quote = F)
write.table(Mon_locate_cd,"tables/Mon_locate.cd.txt",sep = "\t",row.names = F,quote = F)
write.table(Mon_exten_uc,"tables/Mon_exten.uc.txt",sep = "\t",row.names = F,quote = F)
write.table(Mon_sever_uc,"tables/Mon_sever.uc.txt",sep = "\t",row.names = F,quote = F)

FCal_cd = foreach(i=1:ncol(corrected_cd),.combine = rbind) %do%  {
  
  pro=colnames(corrected_cd)[i]
  tmp.pro=corrected_cd[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","Fecal calprotectin"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp.pro=tmp.pro[tmp.pro$`Fecal calprotectin`!="Aangevraagd",]
  tmp.pro=tmp.pro[tmp.pro$`Fecal calprotectin`!="<40",]
  tmp.pro=tmp.pro[tmp.pro$`Fecal calprotectin`!="<100",]
  tmp.pro[,3]=as.numeric(as.character(tmp.pro[,3]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method="spearman")
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,estimate=tmp$estimate)
}
FCal_cd$FDR=p.adjust(FCal_cd$Pvalue)

FCal_uc = foreach(i=1:ncol(corrected_uc),.combine = rbind) %do%  {
  
  pro=colnames(corrected_uc)[i]
  tmp.pro=corrected_uc[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","Fecal calprotectin"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp.pro=tmp.pro[tmp.pro$`Fecal calprotectin`!="Aangevraagd",]
  tmp.pro=tmp.pro[tmp.pro$`Fecal calprotectin`!="<40",]
  tmp.pro=tmp.pro[tmp.pro$`Fecal calprotectin`!="<100",]
  tmp.pro[,3]=as.numeric(as.character(tmp.pro[,3]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method="spearman")
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,estimate=tmp$estimate)
}
FCal_uc$FDR=p.adjust(FCal_uc$Pvalue)

CRP_cd = foreach(i=1:ncol(corrected_cd),.combine = rbind) %do%  {
  
  pro=colnames(corrected_cd)[i]
  tmp.pro=corrected_cd[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","CRP"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro=tmp.pro[!tmp.pro$CRP=="Aangevraagd",]
  tmp.pro$CRP[tmp.pro$CRP=="<5"]=0
  tmp.pro$CRP[tmp.pro$CRP=="<0,3"]=0
  tmp.pro$CRP=as.numeric(tmp.pro$CRP)
  tmp.pro$CRP[tmp.pro$CRP<=5]=0
  tmp.pro$CRP[tmp.pro$CRP>5]=1
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp.pro[,3]=as.numeric(as.character(tmp.pro[,3]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method = "spearman")
  
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,R=tmp$estimate)
}
CRP_cd$FDR=p.adjust(CRP_cd$Pvalue)
CRP_uc = foreach(i=1:ncol(corrected_uc),.combine = rbind) %do%  {
  
  pro=colnames(corrected_uc)[i]
  tmp.pro=corrected_uc[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","CRP"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro=tmp.pro[!tmp.pro$CRP=="Aangevraagd",]
  tmp.pro$CRP[tmp.pro$CRP=="<5"]=0
  tmp.pro$CRP[tmp.pro$CRP=="<0,3"]=0
  tmp.pro$CRP=as.numeric(tmp.pro$CRP)
  tmp.pro$CRP[tmp.pro$CRP<=5]=0
  tmp.pro$CRP[tmp.pro$CRP>5]=1
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp.pro[,3]=as.numeric(as.character(tmp.pro[,3]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method = "spearman")
  
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,R=tmp$estimate)
}
CRP_uc$FDR=p.adjust(CRP_uc$Pvalue)

CRP_cd$Type="CD"
CRP_uc$Type="UC"
CRP_all=rbind(CRP_cd,CRP_uc)
feature=union(CRP_cd$Protein[CRP_cd$FDR<0.1],CRP_uc$Protein[CRP_uc$FDR<0.1])
CRP_all=CRP_all[CRP_all$Protein %in% feature,]

hclust_plot=data.frame(CD=CRP_cd$R,UC=CRP_uc$R,row.names = CRP_cd$Protein)
hclust_plot=hclust_plot[rownames(hclust_plot) %in% feature,]
cluster=hclust(dist(hclust_plot),method = "complete")
pdf("Plot/CRP.cluster.pdf",width = 10,height = 5)
plot(cluster,hang = -1)
dev.off()
protein_order=rownames(hclust_plot)[cluster$order]
CRP_all$Protein=factor(CRP_all$Protein,levels = protein_order)

ggplot(CRP_all, aes(x=Protein, y=Type, fill=R)) + 
  geom_tile(aes(fill=R), colour="white") + 
  scale_fill_gradient2(low="#456BB3", high = "#F26A55", mid = "grey", midpoint = 0) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y= "") + 
  coord_fixed(ratio=1)
ggsave("Plot/CRP.correlation.pdf",width = 9,height = 5)

tmp=pro_all_distance[,c("IL6","CSF-1","SCF","DNER"),drop=F]
tmp=merge(tmp,pheno_all_distance[,c("UMCG1000IBD-ID","CRP")],by.x = "row.names",by.y = "UMCG1000IBD-ID")
rownames(tmp)=tmp$Row.names
tmp$Row.names=NULL
tmp_plot=data.frame(Value=NA,Protein=NA,CRP=NA)
for(i in 1:4){
  tmp.data=data.frame(row.names = rownames(tmp),Value=tmp[,i],Protein=colnames(tmp)[i],CRP=tmp$CRP)
  tmp_plot=rbind(tmp_plot,tmp.data)
}
tmp_plot=tmp_plot[-1,]
tmp_plot=tmp_plot[!tmp_plot$CRP=="Aangevraagd",]
tmp_plot$CRP[tmp_plot$CRP=="<5"]=0
tmp_plot$CRP[tmp_plot$CRP=="<0,3"]=0
tmp_plot$CRP=as.numeric(tmp_plot$CRP)
tmp_plot$CRP[tmp_plot$CRP<=5]=0
tmp_plot$CRP[tmp_plot$CRP>5]=1
tmp_plot=na.omit(tmp_plot)
tmp_plot$CRP=as.factor(tmp_plot$CRP)
tmp_plot$Protein=factor(tmp_plot$Protein,levels = c("IL6","CSF-1","SCF","DNER"))
ggplot(tmp_plot, aes((CRP), Value,color=Protein,fill=Protein)) +
  geom_point(alpha = 0.5, position = "jitter",size=1) +
  geom_boxplot(alpha = 0)+
  facet_wrap( .~ Protein, scales="free_y")+scale_color_jama()+
  theme_bw()
ggsave("./Plot/CRP.Four.genes.pdf")
write.table(CRP_all,"tables/CRP.txt",sep = "\t",quote = F,row.names = F)

season_CD = foreach(i=1:ncol(corrected_cd),.combine = rbind) %do%  {
  i=colnames(corrected_cd)[i]
  tmp.pro=corrected_cd[,i,drop=F]
  tmp=merge(pheno_all_distance[,c("UMCG1000IBD-ID","month")],tmp.pro,by.x="UMCG1000IBD-ID",by.y="row.names",all = F)
  rownames(tmp)=tmp$`UMCG1000IBD-ID`
  tmp$`UMCG1000IBD-ID`=NULL
  tmp=apply(tmp,2,function(x){
    x=as.numeric(as.character(x))
  })
  tmp=as.data.frame(tmp)
  fit <- cosinor.lm(tmp[,2] ~ time(month), data = tmp, period = 12)
  Pvalue=cosinor.detect(fit)[4]
  mm=fit$coefficients
  return.string = data.frame(protein = i, amp=mm[2],acr=mm[3],P=Pvalue)
}
season_CD$FDR=p.adjust(season_CD$P)

tmp.pro=pro_all_distance[,"IL6",drop=F]
tmp=merge(pheno_all_distance[,c("UMCG1000IBD-ID","month")],tmp.pro,by.x="UMCG1000IBD-ID",by.y="row.names",all = F)
rownames(tmp)=tmp$`UMCG1000IBD-ID`
tmp$`UMCG1000IBD-ID`=NULL
mm=rownames(tmp)
tmp=apply(tmp,2,function(x){
  x=as.numeric(as.character(x))
})
tmp=as.data.frame(tmp)
fit <- cosinor.lm(tmp[,2] ~ time(month), data = tmp, period = 12)
rownames(tmp)=mm
tmp=merge(tmp,pheno_all_distance[,c("UMCG1000IBD-ID","year")],by.x = "row.names",by.y = "UMCG1000IBD-ID")
tmp=tmp[tmp$year %in% c(2011,2012,2013),]
tmp$year=as.factor(tmp$year)
ggplot.cosinor.lm(fit)+
  geom_point(data=tmp, aes(y=`IL6`, x=month,color=year))+facet_grid(year ~. )+theme_bw()+scale_color_jama()
ggsave("Plot/Season.pdf",width = 10,height = 8)


# ======================================  calssification  ======================================
library(caret)
library(pROC)
outcome=pheno_all_distance[pheno_all_distance$`IBD type`=="CD" | pheno_all_distance$`IBD type`=="UC",,drop=F]
predictors=pro_all_distance[rownames(pro_all_distance) %in% outcome$`UMCG1000IBD-ID`,]
covariate=outcome[,c("UMCG1000IBD-ID","Age","Gender","BMI","Current smoking"),drop=F]
outcome=outcome[,c("UMCG1000IBD-ID","IBD type"),drop=F]
predictors=merge(predictors,covariate,by.x="row.names",by.y = "UMCG1000IBD-ID",all=F)
rownames(predictors)=predictors$Row.names
predictors$Row.names=NULL
rownames(outcome)=outcome$`UMCG1000IBD-ID`
outcome$`UMCG1000IBD-ID`=NULL
predictors=predictors[order(rownames(predictors)),]
outcome=outcome[order(rownames(outcome)),,drop=F]
outcome[outcome=="CD"]=0
outcome[outcome=="UC"]=1
outcome$`IBD type`=as.numeric(outcome$`IBD type`)

set.seed(12)
inTrain <- createDataPartition(y=outcome$`IBD type`,p=0.8,list=F)
trainSet_predict <- predictors[inTrain,]
trainSet_predict$Gender[trainSet_predict$Gender=="Female"]=0
trainSet_predict$Gender[trainSet_predict$Gender=="Male"]=1
trainSet_predict$Gender=as.numeric(trainSet_predict$Gender)
trainSet_predict = as.data.frame(apply(trainSet_predict,2,function(x){
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
}))
trainSet_outcome <- outcome[inTrain,,drop=F]
trainSet=cbind(trainSet_outcome,trainSet_predict)

testSet_predict <- predictors[-inTrain,]
testSet_predict$Gender[testSet_predict$Gender=="Female"]=0
testSet_predict$Gender[testSet_predict$Gender=="Male"]=1
testSet_predict$Gender=as.numeric(testSet_predict$Gender)
testSet_predict = as.data.frame(apply(testSet_predict,2,function(x){
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
}))
testSet_outcome <- outcome[-inTrain,,drop=F]
testSet=cbind(testSet_outcome,testSet_predict)

# elastic net, 10 fold, "Age","Gender","BMI","Current smoking"
registerDoParallel(cores = 4)
a <- seq(0.1, 0.9, 0.05)
response=as.matrix(trainSet_outcome$`IBD type`)
feature=as.matrix(trainSet_predict[,c("Age","Gender","BMI","Current smoking")])
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(feature, response, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(feature, response, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
selection=as.matrix(coef(md3))
selection=selection[selection[,1]!=0,,drop=F]
print(selection[-1,])
new_response=as.matrix(testSet_outcome$`IBD type`)
new_feature=as.matrix(testSet_predict[,c("Age","Gender","BMI","Current smoking")])
modelroc1_train=roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")))
modelroc_plot1_train=data.frame(FPR=1-modelroc1_train$specificities,TPR=modelroc1_train$sensitivities)
modelroc1_test=roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")))
modelroc_plot1_test=data.frame(FPR=1-modelroc1_test$specificities,TPR=modelroc1_test$sensitivities)
roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )
roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )

# elastic net, 10 fold, all features
registerDoParallel(cores = 4)
a <- seq(0.1, 0.9, 0.05)
response=as.matrix(trainSet_outcome$`IBD type`)
feature=as.matrix(trainSet_predict)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(feature, response, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(feature, response, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
selection=as.matrix(coef(md3))
selection=selection[selection[,1]!=0,,drop=F]
print(selection[-1,])
new_response=as.matrix(testSet_outcome$`IBD type`)
new_feature=as.matrix(testSet_predict)
modelroc2_train=roc(as.numeric(response), as.numeric(predict(md3, feature, type = "response")))
modelroc_plot2_train=data.frame(FPR=1-modelroc2_train$specificities,TPR=modelroc2_train$sensitivities)
modelroc2_test=roc(as.numeric(new_response), as.numeric(predict(md3, new_feature, type = "response")))
modelroc_plot2_test=data.frame(FPR=1-modelroc2_test$specificities,TPR=modelroc2_test$sensitivities)
roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )
roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )

modelroc_plot1_train$Set="Age+Gender+Smoking"
modelroc_plot2_train$Set="Age+Gender+Smoking+Proteins"
modelroc_plot_train=rbind(modelroc_plot1_train,modelroc_plot2_train)
ggplot(modelroc_plot_train,aes(FPR,TPR,color=Set))+geom_line(size = 1, alpha = 0.7)+geom_segment(aes(x = 1, y = 1, xend = 0,yend = 0))+
  labs(title= "ROC curve", 
       x = "False Positive Rate (1-Specificity)", 
       y = "True Positive Rate (Sensitivity)")+theme_minimal()+scale_color_nejm()
ggsave("Plot/Train.model.pdf")
modelroc_plot1_test$Set="Age+Gender+Smoking"
modelroc_plot2_test$Set="Age+Gender+Smoking+Proteins"
modelroc_plot_test=rbind(modelroc_plot1_test,modelroc_plot2_test)
ggplot(modelroc_plot_test,aes(FPR,TPR,color=Set))+geom_line(size = 1, alpha = 0.7)+geom_segment(aes(x = 1, y = 1, xend = 0,yend = 0))+
  labs(title= "ROC curve", 
       x = "False Positive Rate (1-Specificity)", 
       y = "True Positive Rate (Sensitivity)")+theme_minimal()+scale_color_nejm()
ggsave("Plot/Test.model.pdf")


# ======================================  MR  ======================================
library(TwoSampleMR)
library(MRInstruments)

cis_pqtl=read.table("pQTLs/Cis.pQTL.all.txt",header = T,sep = "\t",stringsAsFactors = F)
trans_pqtl=read.table("pQTLs/Trans.pQTL.all.txt",header = T,sep = "\t",stringsAsFactors = F)
pqtl=rbind(cis_pqtl,trans_pqtl)

effect_allele_cis=read.table("pQTLs/Effect_allele.cis.txt",sep = "\t",stringsAsFactors = F,header = T)
effect_allele_trans=read.table("pQTLs/Effect_allele.trans.txt",sep = "\t",stringsAsFactors = F,header = T)
effect_allele=rbind(effect_allele_cis,effect_allele_trans)
effect_allele=effect_allele[!duplicated(effect_allele$SNPName_GWAS1),]
pqtl=merge(pqtl,effect_allele,by.x="SNP",by.y="SNPName_GWAS1",all=F)
colnames(pqtl)=c("SNP","Phenotye","chr","pos","N","P","beta","Q","I","effect_allele")

# estimate standard error (library(gCMAP))
pqtl$se=NULL
pqtl$se=abs(pqtl$beta/zScores(pqtl$P,tails = 2))

# define CCL25 as exposure, for example
tmp_exposure=pqtl[pqtl$Protein=="CCL25",]
tmp_exposure=tmp_exposure[tmp_exposure$SNP %in% c("rs2032887","rs485073"),]
IBD_exp_dat <- format_data(tmp_exposure, type="exposure")

# outcome
ao <- available_outcomes()
chd_out_dat <- extract_outcome_data(
  snps = IBD_exp_dat$SNP,
  outcomes = 'ebi-a-GCST004131'
)

# harmonize
dat <- harmonise_data(
  exposure_dat = IBD_exp_dat,
  outcome_dat = chd_out_dat
)
dat$mr_keep="TRUE"
dat$mr_keep=as.logical(dat$mr_keep)

# run analysis
res <- mr(dat,method_list=c("mr_egger_regression", "mr_ivw","mr_wald_ratio"))
res_single <- mr_singlesnp(dat, single_method="mr_meta_fixed")

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
res_loo <- mr_leaveoneout(dat)


# for all
proteins=unique(pqtl$Phenotye)
uniq_pqtl=pqtl[1,]
for(i in proteins){
  tmp=pqtl[pqtl$Phenotye==i,]
  tmp=tmp[order(tmp$P),]
  tmp=tmp[1,]
  uniq_pqtl=rbind(uniq_pqtl,tmp)
}
uniq_pqtl=uniq_pqtl[-1,]
uniq_pqtl=rbind(uniq_pqtl,pqtl[pqtl$SNP=="rs62124660",])
uniq_pqtl=rbind(uniq_pqtl,pqtl[pqtl$SNP=="rs602662",])
uniq_pqtl=rbind(uniq_pqtl,pqtl[pqtl$SNP=="rs2735477",])
uniq_pqtl=rbind(uniq_pqtl,pqtl[pqtl$SNP=="rs2071442",])
uniq_pqtl=rbind(uniq_pqtl,pqtl[pqtl$SNP=="rs6065729",])
uniq_pqtl=rbind(uniq_pqtl,pqtl[pqtl$SNP=="rs80277428",])

MR_CD = foreach(i=1:nrow(uniq_pqtl),.combine = rbind) %do%  {
  tmp=uniq_pqtl[i,]
  IBD_exp_dat <- format_data(tmp, type="exposure")
  chd_out_dat <- extract_outcome_data(
    snps = IBD_exp_dat$SNP,
    outcomes = 'ebi-a-GCST004132'
  )
  if(!is.null(chd_out_dat)){
    dat <- harmonise_data(
      exposure_dat = IBD_exp_dat,
      outcome_dat = chd_out_dat
    )
    dat$mr_keep="TRUE"
    dat$mr_keep=as.logical(dat$mr_keep)
    res_single <- mr_singlesnp(dat, single_method="mr_meta_fixed")
    dat$WaldRatioPvalue=res_single$p[1]
    dat$WaldRatioBeta=res_single$b[1]
    dat$WaldRatioSE=res_single$se[1]
    
    dat$Protein=uniq_pqtl$Phenotye[i]
    cat(yellow("Number",i,"MR processing","\n"))
    
    return.string = dat
  }else{cat(green(uniq_pqtl$Phenotye[i],"no MR ===================","\n"))}

}
MR_UC = foreach(i=1:nrow(uniq_pqtl),.combine = rbind) %do%  {
  tmp=uniq_pqtl[i,]
  IBD_exp_dat <- format_data(tmp, type="exposure")
  chd_out_dat <- extract_outcome_data(
    snps = IBD_exp_dat$SNP,
    outcomes = 'ebi-a-GCST004133'
  )
  if(!is.null(chd_out_dat)){
    dat <- harmonise_data(
      exposure_dat = IBD_exp_dat,
      outcome_dat = chd_out_dat
    )
    dat$mr_keep="TRUE"
    dat$mr_keep=as.logical(dat$mr_keep)
    res_single <- mr_singlesnp(dat, single_method="mr_meta_fixed")
    dat$WaldRatioPvalue=res_single$p[1]
    dat$WaldRatioBeta=res_single$b[1]
    dat$WaldRatioSE=res_single$se[1]
    
    dat$Protein=uniq_pqtl$Phenotye[i]
    cat(yellow("Number",i,"MR processing","\n"))
    
    return.string = dat
  }else{cat(green(uniq_pqtl$Phenotye[i],"no MR ===================","\n"))}
  
}
MR_IBD = foreach(i=1:nrow(uniq_pqtl),.combine = rbind) %do%  {
  tmp=uniq_pqtl[i,]
  IBD_exp_dat <- format_data(tmp, type="exposure")
  chd_out_dat <- extract_outcome_data(
    snps = IBD_exp_dat$SNP,
    outcomes = 'ebi-a-GCST004131'
  )
  if(!is.null(chd_out_dat)){
    dat <- harmonise_data(
      exposure_dat = IBD_exp_dat,
      outcome_dat = chd_out_dat
    )
    dat$mr_keep="TRUE"
    dat$mr_keep=as.logical(dat$mr_keep)
    res_single <- mr_singlesnp(dat, single_method="mr_meta_fixed")
    dat$WaldRatioPvalue=res_single$p[1]
    dat$WaldRatioBeta=res_single$b[1]
    dat$WaldRatioSE=res_single$se[1]
    
    dat$Protein=uniq_pqtl$Phenotye[i]
    cat(yellow("Number",i,"MR processing","\n"))
    
    return.string = dat
  }else{cat(green(uniq_pqtl$Phenotye[i],"no MR ===================","\n"))}
  
}

MR_all=rbind(MR_CD,MR_UC,MR_IBD)
MR_all$FDR=p.adjust(MR_all$WaldRatioPvalue)
write.table(MR_all,file = "Single_snp.MR.txt",sep = "\t",row.names = F,quote = F)

# ======================================  RNA-seq  ======================================
library(lme4)
coupling=read.table("Coupling.genetics.txt",sep = "\t",stringsAsFactors = F,header = F)
genetics=read.table("FUT2.zoom/FUT2.CCL25.genotype.txt",header = T,sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F)
genetics=as.data.frame(t(genetics))
genetics=genetics[rownames(genetics) %in% coupling$V1,]
coupling=coupling[coupling$V1 %in% rownames(genetics),]
genetics=merge(genetics,coupling,by.x="row.names",by.y="V1")

expression=read.table("Biosy_rna/TMM_expression.table.Log2Transformed.ProbesCentered.SamplesZTransformed.18PCAsOverSamplesRemoved.txt",
                      sep="\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)
# CCL25 is ENSG00000131142.9
# FUT2 is ENSG00000176920.10
expression_sub=as.data.frame(t(expression))
expression_sub=expression_sub[,colnames(expression_sub) %in% c("ENSG00000131142","ENSG00000176920"),drop=F]
write.table("Biosy_rna/old_release.txt",sep = "\t",row.names = T,quote = F)

# cis of CCL25
patient=read.table("Biosy_rna/Patients.biopsy.txt",sep = "\t",stringsAsFactors = F,header = T)
expression_sub=merge(expression_sub,patient[,c("Location","ID","UMCG_ID")],by.x="row.names",by.y="ID",all=F)

genetics=merge(genetics,expression_sub,by.x="V2",by.y="UMCG_ID",all=F)
genetics$rs2032887=as.character(genetics$rs2032887)
genetics$rs2032887[genetics$rs2032887=="A/A"]=0
genetics$rs2032887[genetics$rs2032887=="G/A"]=1
genetics$rs2032887[genetics$rs2032887=="G/G"]=2

genetics$rs3745387=as.character(genetics$rs3745387)
genetics$rs3745387[genetics$rs3745387=="G/G"]=0
genetics$rs3745387[genetics$rs3745387=="A/G"]=1
genetics$rs3745387[genetics$rs3745387=="A/A"]=2

genetics$rs602662=as.character(genetics$rs602662)
genetics$rs602662[genetics$rs602662=="A/A"]=0
genetics$rs602662[genetics$rs602662=="G/A"]=1
genetics$rs602662[genetics$rs602662=="G/G"]=2

genetics_ileum=genetics[genetics$Location=="ileum",]
lmeModel = lmer( ENSG00000131142~ as.numeric(rs3745387)+
                  (1|V2), data=genetics_ileum)
coefs <- data.frame(coef(summary(lmeModel)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))

lmeModel = lmer( ENSG00000131142~ as.numeric(rs2032887)+
                   (1|V2), data=genetics_ileum)
coefs <- data.frame(coef(summary(lmeModel)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))

lmeModel = lmer( ENSG00000176920~ as.numeric(rs602662)+
                   (1|V2), data=genetics_ileum)
coefs <- data.frame(coef(summary(lmeModel)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))

lmeModel = lmer( ENSG00000131142~ as.numeric(rs602662)+
                   (1|V2), data=genetics_ileum)
coefs <- data.frame(coef(summary(lmeModel)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))

cor.test(genetics_ileum$ENSG00000131142,genetics_ileum$ENSG00000176920)

ggplot(genetics_ileum, aes(rs3745387, ENSG00000131142,color=rs3745387,fill=rs3745387)) +
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+theme_bw()+scale_fill_jama()+guides(fill=F,color=F)
ggsave("Plot/CCL25.cis.eQTL.pdf",width = 5,height = 5)
ggplot(genetics_ileum, aes(rs2032887, ENSG00000131142,color=rs2032887,fill=rs2032887)) +
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+theme_bw()+scale_fill_jama()+guides(fill=F,color=F)
ggsave("Plot/CCL25.cis.eQTL.rs2032887.pdf",width = 5,height = 5)
ggplot(genetics_ileum, aes(rs602662, ENSG00000176920,color=rs602662,fill=rs602662)) +
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+theme_bw()+scale_fill_jama()+guides(fill=F,color=F)
ggsave("Plot/FUT2.cis.eQTL.pdf",width = 5,height = 5)
ggplot(genetics_ileum, aes(rs602662, ENSG00000131142,color=rs602662,fill=rs602662)) +
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+theme_bw()+scale_fill_jama()+guides(fill=F,color=F)
ggsave("Plot/FUT2.trans.eQTL.pdf",width = 5,height = 5)

# ======================================  all eQTL  ======================================
library(lme4)
library(foreach)
expression=read.table("Biosy_rna/TMM_expression.table.Log2Transformed.ProbesCentered.SamplesZTransformed.18PCAsOverSamplesRemoved.txt",
                      sep="\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)
expression_sub=as.data.frame(t(expression))
patient=read.table("Biosy_rna/Patients.biopsy.txt",sep = "\t",stringsAsFactors = F,header = T)
coupling=read.table("Coupling.genetics.txt",sep = "\t",stringsAsFactors = F,header = F)
pQTLs=read.table("Unique.pQTLs.txt",header = T,stringsAsFactors = F,sep = "\t")
expression_sub=expression_sub[,colnames(expression_sub) %in% pQTLs$ID,drop=F]
genetics=read.table("FUT2.zoom/pQTL.genotypes.txt",sep = "\t",row.names = 1,stringsAsFactors = F,check.names = F,header = T)
genetics=as.data.frame(t(genetics))
#expression_sub=merge(expression_sub,patient[,c("Location","ID","UMCG_ID")],by.x="row.names",by.y="ID",all=F)

pQTLs=pQTLs[pQTLs$Gene!="CST5",]

eQTL = foreach(i=1:nrow(pQTLs),.combine = rbind) %do%  {

  gene=pQTLs$ID[i]
  tmp.gene=expression_sub[,gene,drop=F]
  tmp.gene=merge(tmp.gene,patient[,c("Location","ID","UMCG_ID")],by.x="row.names",by.y="ID",all=F)
  
  snp=pQTLs$SNP[i]
  tmp.snp=genetics[,snp,drop=F]
  tmp.snp=merge(tmp.snp,coupling,by.x="row.names",by.y="V1",all = F)
  
  tmp.data=merge(tmp.snp,tmp.gene,by.x="V2",by.y="UMCG_ID",all=F)
  tmp.ileum=tmp.data[tmp.data$Location=="ileum",]
  tmp.colon=tmp.data[tmp.data$Location!="ileum",]
  
  lmeModel_ileum = lmer( tmp.ileum[,5]~ tmp.ileum[,3]+
                     (1|V2), data=tmp.ileum)
  lmeModel_colon = lmer( tmp.colon[,5]~ tmp.colon[,3]+
                           (1|V2), data=tmp.colon)
  aa=data.frame(coef(summary(lmeModel_ileum)))
  aa$Pvalue= 2 * (1 - pnorm(abs(aa$t.value)))
  bb=data.frame(coef(summary(lmeModel_colon)))
  bb$Pvalue= 2 * (1 - pnorm(abs(bb$t.value)))
  
  return.string=data.frame(Probe=gene,SNP=snp,
                           ileum_beta=aa$Estimate[2],
                           ileum_se=aa$Std..Error[2],
                           ileum_p=aa$Pvalue[2],
                           colon_beta=bb$Estimate[2],
                           colon_se=bb$Std..Error[2],
                           colon_p=bb$Pvalue[2])
  
}
eQTL=merge(eQTL,pQTLs[,c("ID","Gene")],by.x = "Probe",by.y="ID")
write.table(eQTL,"tables/eQTLs.mucosal.txt",sep = "\t",row.names = F,quote = F)


# ======================================  FUT2 story  ======================================
ld_file_cis=read.table("FUT2.zoom/cis.ld.info.txt",sep = "\t",header = T,stringsAsFactors = F,fill = T)
ccl25_pqtl=read.table("FUT2.zoom/CCL25.zoom.txt",header = T,stringsAsFactors = F,sep = " ")
ccl25_pqtl$P.value=-log10(ccl25_pqtl$P.value)
ccl25_pqtl=merge(ccl25_pqtl,ld_file_cis,by.x="MarkerName",by.y="SNP_B",all=F)
ccl25_pqtl$shape=NA
ccl25_pqtl$shape[ccl25_pqtl$MarkerName=="rs2032887"]=23
ccl25_pqtl$shape[ccl25_pqtl$MarkerName!="rs2032887"]=21
ccl25_pqtl$shape[ccl25_pqtl$MarkerName=="rs3745387"]=23
ccl25_pqtl$shape=as.factor(ccl25_pqtl$shape)
ggplot(ccl25_pqtl, aes(Pos, P.value,fill=R2)) +
  geom_point(aes(shape=shape), color = "black", size = 2)+theme_bw()+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)
ggsave("Plot/CCL25.cis.pQTL.pdf",width = 10,height = 2)

ld_file_trans=read.table("FUT2.zoom/trans.ld.info.txt",sep = "\t",header = T,stringsAsFactors = F,fill = T)
fut2_eqtl=read.table("FUT2.zoom/FUT2.eQTL.zoom.txt",header = F,stringsAsFactors = F,sep = " ")
colnames(fut2_eqtl)=c("SNP","Pvalue","Pos")
fut2_eqtl$Pvalue=-log10(fut2_eqtl$Pvalue)
fut2_eqtl=merge(fut2_eqtl,ld_file_trans,by.x="SNP",by.y="SNP_B",all=F)
fut2_eqtl$shape=NA
fut2_eqtl$shape[fut2_eqtl$SNP=="rs602662"]=23
fut2_eqtl$shape[fut2_eqtl$SNP!="rs602662"]=21
fut2_eqtl$shape=as.factor(fut2_eqtl$shape)
ggplot(fut2_eqtl, aes(Pos, Pvalue,fill=R2)) +
  geom_point(aes(shape=shape), color = "black", size = 2)+theme_bw()+xlim(49142791,49267738)+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(23,21))+guides(shape=F)
ggsave("Plot/FUT2.eQTL.pdf",width = 10,height = 2)

fut2_gwas=read.table("FUT2.zoom/FUT2.gwas.zoom.txt",header = F,stringsAsFactors = F,sep = "\t")
fut2_gwas=merge(fut2_gwas,ld_file_trans,by.x="V16",by.y="SNP_B",all=F)
fut2_gwas$shape=NA
fut2_gwas$shape[fut2_gwas$V16=="rs602662"]=23
fut2_gwas$shape[fut2_gwas$V16!="rs602662"]=21
fut2_gwas$shape=as.factor(fut2_gwas$shape)
ggplot(fut2_gwas, aes(V2, V15,fill=R2)) +
  geom_point(aes(shape=shape), color = "black", size = 2)+theme_bw()+xlim(49142791,49267738)+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)
ggsave("Plot/FUT2.gwas.pdf",width = 10,height = 2)

fut2_CCL25=read.table("FUT2.zoom/CCL25.trans.FUT2.txt",header = T,stringsAsFactors = F,sep = "\t")
fut2_CCL25$P=-log10(fut2_CCL25$P)
fut2_CCL25=merge(fut2_CCL25,ld_file_trans,by.x="SNP",by.y="SNP_B",all=F)
fut2_CCL25$shape=NA
fut2_CCL25$shape[fut2_CCL25$SNP=="rs602662"]=23
fut2_CCL25$shape[fut2_CCL25$SNP!="rs602662"]=21
fut2_CCL25$shape=as.factor(fut2_CCL25$shape)
ggplot(fut2_CCL25, aes(BP, P,fill=R2)) +
  geom_point(aes(shape=shape), color = "black", size = 2)+theme_bw()+xlim(49142791,49267738)+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)
ggsave("Plot/FUT2.CCL25.trans.pdf",width = 10,height = 2)

# genetic score
coupling=read.table("Coupling.genetics.txt",sep = "\t",stringsAsFactors = F,header = F)
genetics=read.table("FUT2.zoom/FUT2.CCL25.genotype.txt",header = T,sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F)
genetics=as.data.frame(t(genetics))
genetics=genetics[rownames(genetics) %in% coupling$V1,]
coupling=coupling[coupling$V1 %in% rownames(genetics),]
genetics=merge(genetics,coupling,by.x="row.names",by.y="V1")

corrected_cd=read.table("CorrectedData/CD_numeric.txt",row.names = 1,stringsAsFactors = F,header = T,sep = "\t")
corrected_uc=read.table("CorrectedData/UC_numeric.txt",row.names = 1,stringsAsFactors = F,header = T,sep = "\t")
corrected_cd=as.data.frame(t(corrected_cd))
corrected_uc=as.data.frame(t(corrected_uc))
corrected_data=rbind(corrected_cd,corrected_uc)
genetics[genetics=="0/0"]=NA
genetics=na.omit(genetics)
genetics$Score=NA

for(i in 1:nrow(genetics)){
  if(genetics$rs2032887[i]=="A/A" & genetics$rs3745387[i]=="G/G" & genetics$rs602662[i]=="A/A"){
    genetics$Score[i]=0
  }
  else if(genetics$rs2032887[i]!="A/A" & genetics$rs3745387[i]=="G/G" & genetics$rs602662[i]=="A/A"){
    genetics$Score[i]=1
  }
  else if(genetics$rs2032887[i]=="A/A" & genetics$rs3745387[i]!="G/G" & genetics$rs602662[i]=="A/A"){
    genetics$Score[i]=1
  }
  else if(genetics$rs2032887[i]!="A/A" & genetics$rs3745387[i]!="G/G" & genetics$rs602662[i]=="A/A"){
    genetics$Score[i]=3
  }
  else if(genetics$rs2032887[i]!="A/A" & genetics$rs3745387[i]!="G/G" & genetics$rs602662[i]!="A/A"){
    genetics$Score[i]=4
  }
  else if(genetics$rs2032887[i]=="A/A" & genetics$rs3745387[i]=="G/G" & genetics$rs602662[i]!="A/A"){
    genetics$Score[i]=2
  }
}

genetics=merge(genetics,corrected_data[,c("CCL25"),drop=F],by.x="V2",by.y="row.names")
genetics_plot=na.omit(genetics)

wilcox.test(genetics$CCL25[genetics$Score==1],genetics$CCL25[genetics$Score==0])
wilcox.test(genetics$CCL25[genetics$Score==2],genetics$CCL25[genetics$Score==0])
wilcox.test(genetics$CCL25[genetics$Score==3],genetics$CCL25[genetics$Score==0])
wilcox.test(genetics$CCL25[genetics$Score==4],genetics$CCL25[genetics$Score==0])
wilcox.test(genetics$CCL25[genetics$Score==1],genetics$CCL25[genetics$Score==3])
wilcox.test(genetics$CCL25[genetics$Score==4],genetics$CCL25[genetics$Score==3])

genetics_plot$Score=as.factor(genetics_plot$Score)
ggplot(genetics_plot, aes(x=Score, y=CCL25,fill=Score)) + 
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+scale_fill_jama()+theme_bw()+guides(color=F)
ggsave("Plot/Genetics.score.CCL25.pdf",width = 15,height = 5)

genetics_plot$rs3745387=factor(genetics_plot$rs3745387,levels = c("G/G","A/G","A/A"))
ggplot(genetics_plot, aes(x=rs3745387, y=CCL25,fill=rs3745387)) + 
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+scale_fill_jama()+theme_bw()+guides(color=F)
ggsave("Plot/cis.CCL25.rs3745387.pdf",width = 7,height = 5)

genetics_plot$rs2032887=factor(genetics_plot$rs2032887,levels = c("A/A","G/A","G/G"))
ggplot(genetics_plot, aes(x=rs2032887, y=CCL25,fill=rs2032887)) + 
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+scale_fill_jama()+theme_bw()+guides(color=F)
ggsave("Plot/cis.CCL25.rs2032887.pdf",width = 7,height = 5)

genetics_plot$rs602662=factor(genetics_plot$rs602662,levels = c("A/A","G/A","G/G"))
ggplot(genetics_plot, aes(x=rs602662, y=CCL25,fill=rs602662)) + 
  geom_point(position = "jitter",size=3,shape = 21, color = "black") +
  geom_boxplot(alpha = 0,color="black")+scale_fill_jama()+theme_bw()+guides(color=F)
ggsave("Plot/trans.CCL25.rs602662.pdf",width = 7,height = 5)


# ====================================== cell type ============================================
# Note: can't perform any interactions analysis, because RNA-seq, and proteomics, and microbiota are not the same time
expression=read.table("Biosy_rna/NewRelease.GeneCount.txt",
                      sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
count=expression
rownames(count)=lapply(strsplit(as.character(rownames(count)), "\\."), "[", 1)
anno=read.table("Biosy_rna/annotation.file.txt",sep = "\t",header = T,stringsAsFactors = F)
anno=anno[!duplicated(anno$Gene),]
anno=anno[!duplicated(anno$id),]
anno=anno[anno$Gene %in% rownames(count),]
anno=anno[order(anno$Gene),]
count=count[rownames(count) %in% anno$Gene,]
count=count[order(rownames(count)),]
rownames(count)=anno$id

deconvolution=xCellAnalysis(count)
deconvolution=as.data.frame(deconvolution)
deconvolution=as.data.frame(t(deconvolution))
deconvolution=deconvolution[,c("cDC",	"Macrophages M1",	"NK cells",	"pDC",	"Macrophages M2",
                               "CD4+ naive T-cells",	"CD4+ Tcm",	"CD8+ naive T-cells",	"CD8+ Tcm",	
                               "Tgd cells",	"Th2 cells",	"Tregs",	"Th1 cells",
                               "NKT",	"CD8+ Tem",	"CD4+ Tem"),drop=F]
write.table(as.data.frame(t(deconvolution)),"Deconvoluton.cells.txt",sep = "\t",row.names = T,quote = F)

deconvolution=read.table("Biosy_rna/Deconvoluton.cells.Log2Transformed.ProbesCentered.SamplesZTransformed.20PCAsOverSamplesRemoved.txt",
                         sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
deconvolution=as.data.frame(t(deconvolution))
metadata=read.table("Biosy_rna/RNA-seq_samples(n=691).txt",sep = "\t",header = T,stringsAsFactors = F)
celltype=merge(deconvolution,metadata,by.x="row.names",by.y="RNA.seqID",all=F)
celltype=merge(celltype,genetics,by.x="Research.ID",by.y="V2",all=F)
celltype$rs2032887=as.character(celltype$rs2032887)
celltype$rs2032887[celltype$rs2032887=="A/A"]=0
celltype$rs2032887[celltype$rs2032887=="G/A"]=1
celltype$rs2032887[celltype$rs2032887=="G/G"]=2
celltype$rs2032887=as.numeric(celltype$rs2032887)
celltype$rs3745387=as.character(celltype$rs3745387)
celltype$rs3745387[celltype$rs3745387=="G/G"]=0
celltype$rs3745387[celltype$rs3745387=="A/G"]=1
celltype$rs3745387[celltype$rs3745387=="A/A"]=2
celltype$rs3745387=as.numeric(celltype$rs3745387)

celltype=celltype[celltype$Location=="ileum" & celltype$Inflammation=="No",]
cell_cor=matrix(nrow = 16,ncol = 4)
cell_cor=as.data.frame(cell_cor)
colnames(cell_cor)=c("snp","cell","pvalue","cor")
for(i in 3:18){
  cell=colnames(celltype)[i]
  cell.tmp=celltype[,c(cell,"rs3745387","Research.ID")]
  cell.tmp=na.omit(cell.tmp)
  lmeModel = lm(cell.tmp[,1] ~ cell.tmp[,2], data=cell.tmp)
  coefs <- data.frame(coef(summary(lmeModel)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  cell_cor$snp[i-2]="rs3745387"
  cell_cor$cell[i-2]=cell
  cell_cor$pvalue[i-2]=coefs$p.z[2]
  cell_cor$cor[i-2]=coefs$Estimate[2]
}

cell_cor$FDR=p.adjust(cell_cor$pvalue)

cell_rs2032887=cell_cor
cell_rs3745387=cell_cor

celltype$rs2032887=as.factor(celltype$rs2032887)
ggplot(celltype, aes(x=rs2032887, y=`CD8+ Tcm`,color=rs2032887)) + 
  geom_point(alpha = 0.5, position = "jitter",size=3) +
  geom_boxplot(alpha = 0)+scale_color_jama()+theme_bw()

