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
library("scales")
source("Microbiome.function.R")

# ======================================  data check ======================================
# Proteomics, data import and rename 
pro_all=read.table("data.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
coupling=read.table("Coupling.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
pheno_all=read.table("OLINK Proteomics phenotype data-v1.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
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

# remove QC failed samples perform PCA, based on manhattan or euclidian distance matrix
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
pcoas$qc_warning[pcoas$qc_warning=="TRUE"]="QC.warned"
pcoas$qc_warning[pcoas$qc_warning=="FALSE"]="QC.passed"
ggplot (pcoas, aes(PCoA1,PCoA2,color=qc_warning)) + 
  geom_point() + theme_bw() + 
  xlab(label = paste("PCoA1",PCoAList$Variance1,sep = " ")) + ylab(label = paste("PCoA2",PCoAList$Variance2,sep = " ")) + 
  theme_bw()+scale_color_jama()
ggsave("./Plot/QC.pdf")
ggplot (pcoas, aes(PCoA1,PCoA2,color=`IBD Medication`)) + 
  geom_point() + theme_bw() + 
  xlab(label = paste("PCoA1",PCoAList$Variance1,sep = " ")) + ylab(label = paste("PCoA2",PCoAList$Variance2,sep = " ")) + 
  theme_bw()


# ======================================  factor assess ======================================
# try clustering samples based on proteins
pro_all_distance=pro_all_distance[rownames(pro_all_distance)!="UMCGIBD01506",]
pheno_all_distance=pheno_all
pheno_all_distance=pheno_all_distance[pheno_all_distance$`UMCG1000IBD-ID` %in% rownames(pro_all_distance),]

pro_cd=pro_all_distance[rownames(pro_all_distance) %in% pheno_all_distance$`UMCG1000IBD-ID`[pheno_all_distance$`IBD type`=="CD"],]
pro_uc=pro_all_distance[rownames(pro_all_distance) %in% pheno_all_distance$`UMCG1000IBD-ID`[pheno_all_distance$`IBD type`=="UC"],]
pro_covariate=pheno_all_distance[,c("Age","Gender","BMI","Current smoking","Aminosalicylates","Thiopurines","Steroids","TNF-antagonists",
                                    "Calcineurin inhibitors","Methotrexate","Mycophenolate mofetil","Ustekinumab",
                                    "Vedolizumab","History of colectomy","History of ileocecal resection")]
rownames(pro_covariate)=pheno_all_distance$`UMCG1000IBD-ID`
pro_covariate$Gender[pro_covariate$Gender=="Female"]=1
pro_covariate$Gender[pro_covariate$Gender=="Male"]=2

heatmap_plot=pro_all_distance
heatmap_plot=cbind(heatmap_plot,pheno_all_distance[,c("Gender","IBD type","IBD Medication","PSC"),drop=F])
heatmap_plot=heatmap_plot[heatmap_plot$`IBD type`=="CD" | heatmap_plot$`IBD type`=="UC",]
ggheatmap(percentize(heatmap_plot[,1:92]),colors = colorRampPalette(brewer.pal(3, "RdBu"))(256),dendrogram = "both",
          show_dendrogram = c(TRUE, F),
          row_side_colors = heatmap_plot[, c("Gender","IBD type","IBD Medication","PSC")],
          hclust_method=NA,
          dist_method="euclidean")  

# use both univariable and multivariable analysis to assess factors influence on each protein
# univariable
annov_cd = foreach(i=1:ncol(pro_cd),.combine = rbind) %do%  {
  pro=colnames(pro_cd)[i]
  tmp.pro=pro_cd[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  
  tmp=matrix(nrow = ncol(tmp.cov),ncol = 3)
  tmp=as.data.frame(tmp)
  colnames(tmp)=c("R2","UniVariate.Pvalue","Factor")
  for(cov in 1:ncol(tmp.cov)){
    tmp$Factor[cov]=colnames(tmp.cov)[cov]
    tmp.cor=cor.test(tmp.pro[,1],as.numeric(tmp.cov[,cov]))
    tmp$R2[cov]=tmp.cor$estimate^2
    tmp$UniVariate.Pvalue[cov]=tmp.cor$p.value
  }
  
  tmp.glm=glm(tmp.pro[,1]~tmp.cov$Age+tmp.cov$Gender+tmp.cov$BMI+tmp.cov$Aminosalicylates+tmp.cov$Thiopurines+
                tmp.cov$Steroids+tmp.cov$`TNF-antagonists`+tmp.cov$`Calcineurin inhibitors`+tmp.cov$Methotrexate+
                tmp.cov$`Mycophenolate mofetil`+tmp.cov$Ustekinumab+tmp.cov$Vedolizumab+tmp.cov$`Current smoking`+
                tmp.cov$`History of colectomy`+tmp.cov$`History of ileocecal resection`,family = gaussian)
  tmp.av=anova(tmp.glm)
  tmp.resid=tmp.av$`Resid. Dev`[nrow(tmp.av)]
  tmp.null=tmp.av$`Resid. Dev`[1]
  tmp.explained=(tmp.null - tmp.resid)/tmp.null
  tmp.explained=data.frame(R2=tmp.explained,UniVariate.Pvalue=NA,Factor="All")
  
  tmp=rbind(tmp,tmp.explained)
  tmp$Protein=pro

  return.string = tmp
}
annov_uc = foreach(i=1:ncol(pro_uc),.combine = rbind) %do%  {
  pro=colnames(pro_uc)[i]
  tmp.pro=pro_uc[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  
  tmp=matrix(nrow = ncol(tmp.cov),ncol = 3)
  tmp=as.data.frame(tmp)
  colnames(tmp)=c("R2","UniVariate.Pvalue","Factor")
  for(cov in 1:ncol(tmp.cov)){
    tmp$Factor[cov]=colnames(tmp.cov)[cov]
    tmp.lm=glm(tmp.pro[,1]~as.numeric(tmp.cov[,cov]))
    tmp.cor=cor.test(tmp.pro[,1],as.numeric(tmp.cov[,cov]))
    tmp$R2[cov]=tmp.cor$estimate^2
    tmp$UniVariate.Pvalue[cov]=tmp.cor$p.value
  }
  
  tmp.glm=glm(tmp.pro[,1]~tmp.cov$Age+tmp.cov$Gender+tmp.cov$BMI+tmp.cov$Aminosalicylates+tmp.cov$Thiopurines+
                tmp.cov$Steroids+tmp.cov$`TNF-antagonists`+tmp.cov$`Calcineurin inhibitors`+tmp.cov$Methotrexate+
                tmp.cov$`Mycophenolate mofetil`+tmp.cov$Ustekinumab+tmp.cov$Vedolizumab+tmp.cov$`Current smoking`+
                tmp.cov$`History of colectomy`+tmp.cov$`History of ileocecal resection`,family = gaussian)
  tmp.av=anova(tmp.glm)
  tmp.resid=tmp.av$`Resid. Dev`[nrow(tmp.av)]
  tmp.null=tmp.av$`Resid. Dev`[1]
  tmp.explained=(tmp.null - tmp.resid)/tmp.null
  tmp.explained=data.frame(R2=tmp.explained,UniVariate.Pvalue=NA,Factor="All")
  
  tmp=rbind(tmp,tmp.explained)
  tmp$Protein=pro
  
  return.string = tmp
}

annov_cd_uni=annov_cd[annov_cd$Factor!="All",]
annov_cd_all=annov_cd[annov_cd$Factor=="All",]
annov_uc_uni=annov_uc[annov_uc$Factor!="All",]
annov_uc_all=annov_uc[annov_uc$Factor=="All",]
annov_cd_uni$FDR=p.adjust(annov_cd_uni$UniVariate.Pvalue)
annov_uc_uni$FDR=p.adjust(annov_uc_uni$UniVariate.Pvalue)
annov_cd_uni$R2[annov_cd_uni$FDR>0.1]=0
annov_uc_uni$R2[annov_uc_uni$FDR>0.1]=0

annov_cd_sig=annov_cd_uni[annov_cd_uni$FDR<0.1,]
annov_uc_sig=annov_uc_uni[annov_uc_uni$FDR<0.1,]

# select significant factors and then put them in multivariable model
glm_cd = foreach(i=1:ncol(pro_cd),.combine = rbind) %do%  {
  pro=colnames(pro_cd)[i]
  tmp.pro=pro_cd[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.cov=tmp.cov[,colnames(tmp.cov) %in% annov_cd_sig$Factor[annov_cd_sig$Protein==pro],drop=F]
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  
  if(ncol(tmp.cov)==0){
    
    return.string=data.frame(R2=0,MultiVariate.Pvalue=NA,Factor="All",Protein=pro)
    
  }else{
    
    tmp.glm=glm(tmp.pro[,1]~.,data=tmp.cov,family = gaussian)
    tmp.av=anova(tmp.glm)
    tmp.resid=tmp.av$`Resid. Dev`[nrow(tmp.av)]
    tmp.null=tmp.av$`Resid. Dev`[1]
    tmp.explained=(tmp.null - tmp.resid)/tmp.null
    
    return.string=data.frame(R2=tmp.explained,MultiVariate.Pvalue=NA,Factor="All",Protein=pro)
    
  }

}
glm_uc = foreach(i=1:ncol(pro_uc),.combine = rbind) %do%  {
  pro=colnames(pro_uc)[i]
  tmp.pro=pro_uc[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.cov=tmp.cov[,colnames(tmp.cov) %in% annov_uc_sig$Factor[annov_uc_sig$Protein==pro],drop=F]
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  
  if(ncol(tmp.cov)==0){
    
    return.string=data.frame(R2=0,MultiVariate.Pvalue=NA,Factor="All",Protein=pro)
    
  }else{
    
    tmp.glm=glm(tmp.pro[,1]~.,data=tmp.cov,family = gaussian)
    tmp.av=anova(tmp.glm)
    tmp.resid=tmp.av$`Resid. Dev`[nrow(tmp.av)]
    tmp.null=tmp.av$`Resid. Dev`[1]
    tmp.explained=(tmp.null - tmp.resid)/tmp.null
    
    return.string=data.frame(R2=tmp.explained,MultiVariate.Pvalue=NA,Factor="All",Protein=pro)
    
  }
  
}

sample_order=annov_cd_uni[annov_cd_uni$Factor=="Age",]
sample_order=sample_order[order(sample_order$R2,decreasing = T),]
sample_order=as.character(sample_order$Protein)
category_order=colnames(pro_covariate)
set.seed(12)
palette <- distinctColorPalette(16)
annov_cd_uni$Factor<-factor(annov_cd_uni$Factor,levels = rev(category_order))
annov_uc_uni$Factor<-factor(annov_cd_uni$Factor,levels = rev(category_order))
ggplot() +
  geom_bar(data=annov_cd_uni,aes(x=Protein,y=R2,fill=Factor),position = "stack",stat="identity") +
  scale_x_discrete(limits = sample_order)+
  scale_fill_manual(breaks = category_order,values = palette)+
  geom_point(data=glm_cd, aes(x=Protein,y=R2,group=1),color="grey50")+
  geom_bar(data=annov_uc_uni,aes(x=Protein,y=-R2,fill=Factor), stat="identity") +
  geom_point(data=glm_uc, aes(x=Protein,y=-R2,group=1),color="grey50")+theme_classic()+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 6),
        axis.line.y =  element_line(colour = 'white'))+
  geom_hline(yintercept = 0,  size=1,colour = 'white')
ggsave("Plot/Univariation.pdf",height = 5,width = 10)

matrix_cd=matrix(ncol = 0,nrow = 92)
matrix_cd=as.data.frame(matrix_cd)
rownames(matrix_cd)=colnames(pro_all_distance)
for(i in 1:ncol(pro_covariate)){
  cov=colnames(pro_covariate)[i]
  tmp=annov_cd_uni[annov_cd_uni$Factor==cov,c("Protein","FDR")]
  tmp[is.na(tmp)]=9
  tmp$FDR[tmp$FDR<0.1]=NA
  tmp$FDR[tmp$FDR>0.1]=0
  tmp$FDR[is.na(tmp$FDR)]=1
  rownames(tmp)=tmp$Protein
  tmp$Protein=NULL
  colnames(tmp)=cov
  matrix_cd=cbind(matrix_cd,tmp)
}
pdf("./Plot/Upset.plot.CD.pdf",width = 15,height = 6)
upset(matrix_cd, sets = colnames(matrix_cd), sets.bar.color = "#56B4E9",point.size=3,
      order.by = "freq",  keep.order = TRUE)
dev.off()
matrix_uc=matrix(ncol = 0,nrow = 92)
matrix_uc=as.data.frame(matrix_uc)
rownames(matrix_uc)=colnames(pro_all_distance)
for(i in 1:ncol(pro_covariate)){
  cov=colnames(pro_covariate)[i]
  tmp=annov_uc_uni[annov_uc_uni$Factor==cov,c("Protein","FDR")]
  tmp[is.na(tmp)]=9
  tmp$FDR[tmp$FDR<0.1]=NA
  tmp$FDR[tmp$FDR>0.1]=0
  tmp$FDR[is.na(tmp$FDR)]=1
  rownames(tmp)=tmp$Protein
  tmp$Protein=NULL
  colnames(tmp)=cov
  matrix_uc=cbind(matrix_uc,tmp)
}
pdf("./Plot/Upset.plot.UC.pdf",width = 15,height = 6)
upset(matrix_uc, sets = colnames(matrix_uc), sets.bar.color = "#56B4E9",point.size=3,
      order.by = "freq",  keep.order = TRUE)
dev.off()

combine_cd=merge(pro_covariate,pro_cd,by="row.names")
combine_cd$Type="CD"
combine_uc=merge(pro_covariate,pro_uc,by="row.names")
combine_uc$Type="UC"
combine=rbind(combine_cd,combine_uc)
ggplot(combine, aes(Age, CXCL10)) +
  geom_point(shape = 21, fill = "lightgray",
             color = "black", size = 3)+geom_smooth(method = lm,color="red",fill="red")+facet_grid(~Type)+
  theme_bw()
ggsave("Plot/Age.CXCL10.pdf")
ggplot(combine, aes(Gender, CCL11,color=Gender)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2))+
  facet_grid(~Type)+scale_color_jama()+
  theme_bw()
ggsave("Plot/Sex.CCL11.pdf")
combine$Thiopurines=as.factor(combine$Thiopurines)
ggplot(combine, aes(Thiopurines, Flt3L,color=Thiopurines)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2))+
  facet_grid(~Type)+scale_color_jama()+
  theme_bw()
ggsave("Plot/Thiopurines.Flt3L.pdf")

# ======================================  clinical outcome/genetics correlation ======================================
# correct for all significant covariates
pro_covariate$Gender=as.numeric(pro_covariate$Gender)
corrected_cd = foreach(i=1:ncol(pro_cd),.combine = cbind) %do%  {
  pro=colnames(pro_cd)[i]
  cat(green("Correction ===>",pro,"\n"))
  tmp.pro=pro_cd[,i,drop=F]
  tmp.cov=pro_covariate[rownames(pro_covariate) %in% rownames(tmp.pro),]
  tmp.sig.cov=tmp.cov[,colnames(tmp.cov) %in% annov_cd_sig$Factor[annov_cd_sig$Protein==pro],drop=F]

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
  tmp.sig.cov=tmp.cov[,colnames(tmp.cov) %in% annov_uc_sig$Factor[annov_uc_sig$Protein==pro],drop=F]
  
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

# HBI, SCCAI, and CRP
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
hbi_uc = foreach(i=1:ncol(corrected_uc),.combine = rbind) %do%  {
  
  pro=colnames(corrected_uc)[i]
  tmp.pro=corrected_uc[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","HBI"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method = "spearman")
  
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,R=tmp$estimate)
}
hbi_uc$FDR=p.adjust(hbi_uc$Pvalue)

tmp=pro_all_distance[,colnames(pro_all_distance) %in% as.character(hbi_cd$Protein[hbi_cd$FDR<0.05]),drop=F]
tmp=merge(tmp,pheno_all_distance[pheno_all_distance$`IBD type`=="CD",c("UMCG1000IBD-ID","HBI")],by.x = "row.names",by.y = "UMCG1000IBD-ID")
rownames(tmp)=tmp$Row.names
tmp$Row.names=NULL
tmp_plot=data.frame(Value=NA,Protein=NA,HBI=NA)
for(i in 1:11){
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
ggsave("Plot/HBI.elevenGenes.pdf",width = 6,height = 6)

SCCAI_cd = foreach(i=1:ncol(corrected_cd),.combine = rbind) %do%  {
  
  pro=colnames(corrected_cd)[i]
  tmp.pro=corrected_cd[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","SCCAI"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method = "spearman")
  
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,R=tmp$estimate)
}
SCCAI_cd$FDR=p.adjust(SCCAI_cd$Pvalue)
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

tmp=pro_all_distance[,colnames(pro_all_distance) %in% as.character(SCCAI_uc$Protein[SCCAI_uc$FDR<0.05]),drop=F]
tmp=merge(tmp,pheno_all_distance[pheno_all_distance$`IBD type`=="UC",c("UMCG1000IBD-ID","SCCAI")],by.x = "row.names",by.y = "UMCG1000IBD-ID")
rownames(tmp)=tmp$Row.names
tmp$Row.names=NULL
tmp_plot=data.frame(Value=NA,Protein=NA,SCCAI=NA)
for(i in 1:8){
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
ggsave("Plot/SCCAI.eightGenes.pdf",width = 6,height = 6)

CRP_cd = foreach(i=1:ncol(corrected_cd),.combine = rbind) %do%  {
  
  pro=colnames(corrected_cd)[i]
  tmp.pro=corrected_cd[,i,drop=F]
  tmp.pro=merge(tmp.pro,pheno_all_distance[,c("UMCG1000IBD-ID","CRP"),drop=F],by.x="row.names",by.y="UMCG1000IBD-ID",all=F)
  tmp.pro=na.omit(tmp.pro)
  tmp.pro=tmp.pro[!tmp.pro$CRP %in% c("<0,3","<5","Aangevraagd"),]
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
  tmp.pro=tmp.pro[!tmp.pro$CRP %in% c("<0,3","<5","Aangevraagd"),]
  tmp.pro[,2]=as.numeric(as.character(tmp.pro[,2]))
  tmp.pro[,3]=as.numeric(as.character(tmp.pro[,3]))
  tmp=cor.test(tmp.pro[,2],tmp.pro[,3],method = "spearman")
  
  return.string=data.frame(Protein=pro,Pvalue=tmp$p.value,R=tmp$estimate)
}
CRP_uc$FDR=p.adjust(CRP_uc$Pvalue)

CRP_cd$Type="CD"
CRP_uc$Type="UC"
CRP_all=rbind(CRP_cd,CRP_uc)
feature=union(CRP_cd$Protein[CRP_cd$FDR<0.05],CRP_uc$Protein[CRP_uc$FDR<0.05])
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
ggsave("Plot/CRP.correlation.pdf")

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
tmp_plot=tmp_plot[!tmp_plot$CRP %in% c("<0,3","<5","Aangevraagd"),]
tmp_plot=na.omit(tmp_plot)
tmp_plot$CRP=as.numeric(tmp_plot$CRP)
tmp_plot$Protein=factor(tmp_plot$Protein,levels = c("IL6","CSF-1","SCF","DNER"))
ggplot(tmp_plot, aes(CRP, Value,color=Protein,fill=Protein)) +
  geom_point(shape = 21, 
             color = "black", size = 3)+
  geom_smooth(method = lm)+
  facet_wrap( .~ Protein, scales="free_y")+scale_color_jama()+
  theme_bw()
ggsave("./Plot/CRP.Four.genes.pdf")
