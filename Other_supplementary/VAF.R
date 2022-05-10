library(ggplot2)
library(nlme)
library(lme4)
SNP_file<-read.table('/Users/yw2/Documents/FARM/small_bowel/filter/snp.tsv',header=F, stringsAsFactors = F)
# patient_vcf<-read.table('/Users/yw2/Documents/FARM/yw2/data/small_bowel/filter/PD28690/PD28690.snp.tsv',comment.char = "",header=T, stringsAsFactors = F)
VAF_stat <- data.frame()
# SNP_file=SNP_file[c(30,13),]
which(SNP_file$V1 ==path)
#---------median VAF---------------
for (path in SNP_file$V1[1:40]){
  # path = SNP_file$V1[16]
  patient_vcf<-read.table(paste0('/Users/yw2/Documents/FARM/small_bowel/filter/',path),comment.char = "",header=T, stringsAsFactors = F)
  sample_name<-unlist(strsplit(colnames(patient_vcf[grep('MTR',colnames(patient_vcf))]),split='_MTR'))
  
  all_VAF<-patient_vcf[,colnames(patient_vcf[grep('MTR',colnames(patient_vcf))])]/patient_vcf[,colnames(patient_vcf[grep('DEP',colnames(patient_vcf))])]
  

  df=data.frame()
  for (locus in sample_name[2:length(sample_name)]){
    sample_VAF<-all_VAF[,colnames(all_VAF)[grep(locus,colnames(all_VAF))]]
    sample_VAF<-na.omit(sample_VAF[sample_VAF>0])
    median(sample_VAF)
    VAF_stat <- rbind(VAF_stat,c(locus,median(sample_VAF)),stringsAsFactors = F)
    df_tmp=data.frame(VAF=sample_VAF,locus=locus)
    df <- rbind(df,df_tmp,stringsAsFactors = F)
  }
  print(ggplot(df, aes(x = VAF)) + geom_density(aes(group=factor(locus)), fill="grey", size=1, alpha=.4)+ggtitle(substring(path,1,7))+ylab('Density')+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black')))
}
  
# VAF_stat<-na.omit(VAF_stat)
colnames(VAF_stat)=c('sample','median_VAF')
#--------------coverage-------------
VAF_stat$coverage<-NA
for (i in 1:dim(VAF_stat)[1]){
  if (VAF_stat$sample[i] %in% data$sample){
    VAF_stat$coverage[i] = data$coverage[which(data$sample == VAF_stat$sample[i])]
  }
}
VAF_stat$coverage[VAF_stat$sample=='PD43949b_lo0021']=33.33

tmp<-read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/data/median_VAF_with_indel_20211116.txt',header=T,stringsAsFactors = F)
for (i in 1:length(VAF_stat$coverage)){
  if (VAF_stat$sample[i]%in%tmp$sample){
    VAF_stat$coverage[i]=tmp$coverage[which(tmp$sample==VAF_stat$sample[i])]
  }
}
write.table(VAF_stat,'/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/data/median_VAF_with_indel.txt', quote=F, col.names = T, row.names = F)
# VAF_stat<-read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/data/median_VAF_with_indel_20211116.txt',header=T,stringsAsFactors = F)


#------------sensitivity------------
set.seed(99)
min_reads=4 #CaVEMan minimum number of reads
VAF_stat$sensitivity=NA
VAF_stat[,2]=as.numeric(VAF_stat[,2])
for (i in 1:length(VAF_stat[,1])){
  VAF_stat$sensitivity[i]=mean(unlist(lapply(rpois(n=100000,lambda=VAF_stat$coverage[i]),function(x) pbinom(q=min_reads-0.1,size=x,p=VAF_stat[i,2],lower.tail = F))))
}

min_reads=5 #Pindel minimum number of reads
VAF_stat$indel_sensitivity=NA
for (i in 1:length(VAF_stat[,1])){
  VAF_stat$indel_sensitivity[i]=mean(unlist(lapply(rpois(n=100000,lambda=VAF_stat$coverage[i]),function(x) pbinom(q=min_reads-0.1,size=x,p=VAF_stat[i,2],lower.tail = F))))
}
#-------------SBS count--------------

# VAF_stat$sample == rownames(mut_persample)
VAF_stat$sbs_count=NA


options(stringsAsFactors=F)
library(ape)
library(ggtree)
rownames(VAF_stat)=VAF_stat$sample

for (patient in patients[1:37]){
  tree=read.tree(paste0("/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/data/tree/",patient,"_snp_tree_with_branch_length.tree"))
  tree_df=fortify(tree)
  # tree_df[tree_df$isTip,c("label","x")]
  tmp=tree_df[tree_df$isTip,c("label","x")]
  VAF_stat[tmp$label,'sbs_count'] <- tmp$x 
}

# VAF_stat<-na.omit(VAF_stat)

# mut_persample<-read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/data/mut_persample.txt')
# for (i in 1:length(VAF_stat[,1])){
#   name=as.character(VAF_stat[i,1])
#   VAF_stat$sbs_count[i] = sum(mut_persample[name,])
# }
VAF_stat$sbs_count_adj<- VAF_stat$sbs_count/ VAF_stat$sensitivity

# sum(mut_persample['PD46573b_lo0003',])
# sum(input_for_hdp['PD46573_11',])
#-------------indel count--------------

# VAF_stat$sample == rownames(mut_persample)
VAF_stat$indel_count=NA

for (patient in patients[1:37]){
  tree=read.tree(paste0("/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/data/tree/",patient,"_indel_tree_with_branch_length.tree"))
  tree_df=fortify(tree)
  # tree_df[tree_df$isTip,c("label","x")]
  tmp=tree_df[tree_df$isTip,c("label","x")]
  VAF_stat[tmp$label,'indel_count'] <- tmp$x 
}

VAF_stat<-na.omit(VAF_stat)

# indel_persample<-as.data.frame(t(read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/data/indel/indel_vcf/output/ID/indel_vcf.ID83.all')))
# colnames(indel_persample)=indel_persample[1,]
# indel_persample=indel_persample[-1,]
# indel_sample_name <- unlist(strsplit(indel_persample$MutationType,split = '_good_indel'))
# rownames(indel_persample)=indel_sample_name
# indel_persample=indel_persample[,-1]
# for (i in 1:length(VAF_stat$indel_count)){
#   name=as.character(VAF_stat$sample[i])
#   VAF_stat$indel_count[i] = sum(as.numeric(as.character(indel_persample[name,])))
# }
# VAF_stat$sample[! VAF_stat$sample%in% indel_sample_name]
# length(indel_sample_name)
# length(VAF_stat$sample)

VAF_stat$indel_count_adj<- VAF_stat$indel_count/VAF_stat$indel_sensitivity
VAF_stat$indel_count[which(VAF_stat$sample=='PD42833b_lo0023')]=0
VAF_stat$indel_count_adj[which(VAF_stat$sample=='PD42833b_lo0023')]=0

#---------------age-------------------
VAF_stat$age <- NA
metadata<-read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/data/metadata.txt')
PD = substr(VAF_stat$sample,1,7)
for (i in 1:length(PD)){
  if (PD[i] %in% metadata[,1]){
    VAF_stat$age[i]= metadata[which(metadata[,1]==PD[i]),4]
  }
  else if (PD[i] == 'PD43850'){
    VAF_stat$age[i]= 51
  }
  else if (PD[i] == 'PD43851'){
    VAF_stat$age[i]= 47
  }
}

#------------------region-------------------
VAF_stat$region <- NA
for (i in 1:length(PD)){
  if (PD[i] %in% metadata[,1]){
    VAF_stat$region[i]= as.character(metadata[which(metadata[,1]==PD[i]),2])
  }
}

for (i in 1:length(VAF_stat$region)){
  if (VAF_stat$sample[i]%in%tmp$sample){
    VAF_stat$region[i]=tmp$region[which(tmp$sample==VAF_stat$sample[i])]
  }
}
#-------------------diseases--------------
VAF_stat$condition='Normal'
VAF_stat$patient = substr(VAF_stat$sample,1,7)
VAF_stat$condition[VAF_stat$patient %in% c('PD46562','PD46563','PD46565','PD46566','PD46568','PD46573')] <- 'Coeliac'
VAF_stat$condition[VAF_stat$patient %in% c('PD28690','PD43853')] <- 'Chemo'

#-------------------ref--------------
VAF_stat$ref='Crypt'
VAF_stat$ref[VAF_stat$sample %in% c('PD34200a_i33','PD34200a_i22','PD43403e_lo0017','PD41853c_lo0018','PD41853c_lo0018','PD41851d_lo0059','PD41852c_lo0022','PD28690bp_SB1_C9','PD42834b_lo0035','PD28690bp_SB1_D9','PD46573b_lo0003','PD46562b_lo0001','PD42835b_lo0041')]="Control"
VAF_stat$ref[VAF_stat$sample %in% c('PD43851j_P52_DDM_B4','PD43851j_P52_DDM_B5','PD43851j_P52_DDM_C4','PD43851j_P52_DDM_E2')]="Brunners_gland"
# VAF_stat$ref[VAF_stat$sample %in% c('PD42835b_lo0032')]="Ambiguous"



# colnames(VAF_stat)=c('sample','median_VAF','sensitivity','sbs_count','condition','age','region','coverage','sbs_count_adj')
VAF_stat$project = '2146'
VAF_stat$project[which(VAF_stat$patient %in% c('PD28690','PD43851','PD43850'))]='1697'
VAF_stat$project[which(VAF_stat$patient %in% c('PD34200','PD37266'))]='1494'
VAF_stat$project[which(VAF_stat$patient %in% c('PD37449'))]='1728'
VAF_stat$project[which(VAF_stat$patient %in% c('PD41851','PD41852','PD41853','PD42834','PD42833','PD42835','PD43853'))]='2064'
VAF_stat$project[which(VAF_stat$sample %in% c('PD43402d_lo0001','PD43402d_lo0002','PD43402d_lo0003','PD43402d_lo0004','PD43402d_lo0005','PD43402d_lo0006','PD43402d_lo0007'))]='2064'
VAF_stat$project[which(VAF_stat$sample %in% c('PD41851d_lo0092'))]='2146'



write.table(VAF_stat,'/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/data/median_VAF_with_indel_20220103.txt', quote=F, col.names = T, row.names = F)
# VAF_stat=read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/data/median_VAF_with_indel.txt',header=T,stringsAsFactors = F)
VAF_stat<-VAF_stat[which(VAF_stat$coverage>15),]
VAF_stat<-VAF_stat[which(VAF_stat$ref=='Crypt'),]


sum(VAF_stat$sbs_count)
sum(VAF_stat$indel_count)
mean(VAF_stat$coverage)
length(unique(substr(VAF_stat$sample[which(VAF_stat$region=='Duodenum')],1,7)))
length(unique(substr(VAF_stat$sample[which(VAF_stat$region=='Ileum')],1,7)))
length(unique(substr(VAF_stat$sample[which(VAF_stat$region=='Jejunum')],1,7)))
table(VAF_stat$region)
# colnames(patient_vcf[grep('MTR',colnames(patient_vcf))])
# patient_vcf[,colnames(patient_vcf[grep('DEP',colnames(patient_vcf))])]


# VAF_stat$sbs_count[which(VAF_stat$sbs_count<1000)]<-NA


#---------------modelling mutation burden--------------

exposure_matrix_persample_tmp = exposure_matrix_samples[VAF_stat$sample,]
VAF_stat$SBS18 <- VAF_stat$sbs_count_adj*exposure_matrix_persample_tmp$SBS18
VAF_stat$SBS1 <- VAF_stat$sbs_count_adj*exposure_matrix_persample_tmp$SBS1
VAF_stat$SBS5 <- VAF_stat$sbs_count_adj*exposure_matrix_persample_tmp$SBS5
VAF_stat$SBS2 <- VAF_stat$sbs_count_adj*exposure_matrix_persample_tmp$SBS2
VAF_stat$SBS13 <- VAF_stat$sbs_count_adj*exposure_matrix_persample_tmp$SBS13
VAF_stat$SBS88 <- VAF_stat$sbs_count_adj*exposure_matrix_persample_tmp$SBS88
VAF_stat$SBS35 <- VAF_stat$sbs_count_adj*exposure_matrix_persample_tmp$SBS35
VAF_stat$SBS41 <- VAF_stat$sbs_count_adj*exposure_matrix_persample_tmp$SBS41


sum(VAF_stat$SBS18)/sum(VAF_stat$sbs_count_adj)
# VAF_stat$sbs_count[VAF_stat$sbs_count_adj<1000]<-NA

# library(MuMIn)
# data(Arabidopsis)
# attach(Arabidopsis)

# VAF_stat<-na.omit(VAF_stat)
lmm1 <- lme(sbs_count_adj ~ condition + age+ region, random= ~1|patient,method = "REML",data=na.omit(VAF_stat))
summary(lmm1)
tmp<-summary(lmm1)[["tTable"]]
p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
tmp[,'p-value']<-p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
VAF_stat$fit_sbs_count_adj<-predict(lmm1) 

39.4496+1.96*5.8153
39.4496-1.96*5.8153
table(VAF_stat$region)
length(unique(VAF_stat$patient))
length(unique(VAF_stat$patient[which(VAF_stat$region=='Jejunum')]))
length(unique(VAF_stat$patient[which(VAF_stat$region=='Ileum')]))
length(unique(VAF_stat$patient[which(VAF_stat$region=='Duodenum')]))
mean(VAF_stat$coverage)

cor(VAF_stat$sbs_count_adj,VAF_stat$age,method = 'pearson')
cor(VAF_stat$SBS18,VAF_stat$age,method = 'pearson')
cor(VAF_stat$SBS1,VAF_stat$age,method = 'pearson')
cor(VAF_stat$SBS5,VAF_stat$age,method = 'pearson')
cor(VAF_stat$SBS2,VAF_stat$age,method = 'pearson')
cor(VAF_stat$SBS13,VAF_stat$age,method = 'pearson')
cor(VAF_stat$SBS88,VAF_stat$age,method = 'pearson')
tmp<-exposure_matrix_persample[which(rownames(exposure_matrix_persample) %in% VAF_stat$sample),]
tmp<-exposure_matrix_per_patient
sum(tmp$SBS2+tmp$SBS13>0.05)/nrow(tmp)
# r.squaredGLMM(lmm1)

library(ggplot2)

tmp=data.frame(condition=c(16,3))
tmp_name = c('Coeliac','Normal')
rownames(tmp)=tmp_name


pdf('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/docs/manuscript/Fig/SBS_burden.pdf')
ggplot(data = VAF_stat, mapping = aes(x = age, y = sbs_count_adj, colour = region,fill=region,shape=condition))+geom_point() +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  scale_shape_manual(values = tmp$condition)+
  # geom_line(aes(y=fit_sbs_count_adj,x=age), size=0.8) +
  labs(y='SBS Mutational Burden',x="Age (yrs)")
dev.off()
#--------indel bruden-----------

# lmm1 <- lme(indel_count_adj ~ condition + age+ region, random= ~1|patient,method = "REML",data=na.omit(VAF_stat[VAF_stat$patient != 'PD46573',]))
lmm1 <- lme(indel_count_adj ~ condition + age+ region, random= ~1|patient,method = "REML",data=na.omit(VAF_stat))
summary(lmm1)
tmp<-summary(lmm1)[["tTable"]]
p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
tmp[,'p-value']<-p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
VAF_stat$fit_sbs_count_adj<-predict(lmm1) 

2.43201+1.96*0.64995 
2.43201-1.96*0.64995 

tmp=data.frame(condition=c(16,3))
tmp_name = c('Coeliac','Normal')
rownames(tmp)=tmp_name


pdf('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/docs/manuscript/Fig/indel_burden.pdf')
ggplot(data = VAF_stat, mapping = aes(x = age, y = indel_count_adj, colour = region,fill=region,shape=condition))+geom_point() +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  # geom_line(aes(y=fit_sbs_count_adj,x=age), size=0.8) +
  scale_shape_manual(values = tmp$condition)+
  labs(y='Indel Mutational Burden',x="Age (yrs)")
dev.off()


#-------------SBS1-----------
lmm1 <- lme(SBS1 ~ condition + age+ region, random= ~1|patient,method = "REML",data=na.omit(VAF_stat))
summary(lmm1)
tmp<-summary(lmm1)[["tTable"]]
p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
tmp[,'p-value']<-p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
VAF_stat$fit_SBS1<-predict(lmm1)


ggplot(data = VAF_stat, mapping = aes(x = age, y = SBS1))+geom_point() +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  # geom_line(aes(y=fit_sbs_count_adj,x=age), size=0.8) +
  labs(y='SBS1',x="Age (yrs)")

#-------------SBS5-----------
lmm1 <- lme(SBS5 ~ condition + age+ region, random= ~1|patient,method = "REML",data=na.omit(VAF_stat))
summary(lmm1)
tmp<-summary(lmm1)[["tTable"]]
p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
tmp[,'p-value']<-p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
VAF_stat$fit_SBS5<-predict(lmm1)


ggplot(data = VAF_stat, mapping = aes(x = age, y = SBS5))+geom_point(colour = 'black',fill='black') +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  # geom_line(aes(y=fit_sbs_count_adj,x=age), size=0.8) +
  labs(y='SBS5',x="Age (yrs)")

#-------------SBS18-----------
lmm1 <- lme(SBS18 ~ condition + age+ region, random= ~1|patient,method = "REML",data=na.omit(VAF_stat))
summary(lmm1)
tmp<-summary(lmm1)[["tTable"]]
p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
tmp[,'p-value']<-p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
VAF_stat$fit_SBS18<-predict(lmm1)


ggplot(data = VAF_stat, mapping = aes(x = age, y = SBS18))+geom_point(colour = 'black',fill='black') +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  # geom_line(aes(y=fit_sbs_count_adj,x=age), size=0.8) +
  labs(y='SBS18',x="Age (yrs)")

#-------------SBS2-----------
lmm1 <- lme(SBS2 ~ condition + age+ region, random= ~1|patient,method = "REML",data=na.omit(VAF_stat))
summary(lmm1)
tmp<-summary(lmm1)[["tTable"]]
p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
tmp[,'p-value']<-p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
VAF_stat$fit_SBS2<-predict(lmm1)


# ggplot(data = VAF_stat, mapping = aes(x = age, y = SBS2, colour = region,fill=region))+geom_point() +
ggplot(data = VAF_stat, mapping = aes(x = age, y = SBS2))+geom_point(colour = 'black',fill='black') +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  # geom_line(aes(y=fit_sbs_count_adj,x=age), size=0.8) +
  labs(y='SBS2',x="Age (yrs)")

#-------------SBS13-----------
lmm1 <- lme(SBS13 ~ condition + age+ region, random= ~1|patient,method = "REML",data=na.omit(VAF_stat))
summary(lmm1)
tmp<-summary(lmm1)[["tTable"]]
p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
tmp[,'p-value']<-p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
VAF_stat$fit_SBS13<-predict(lmm1)


# ggplot(data = VAF_stat, mapping = aes(x = age, y = SBS13, colour = region,fill=region))+geom_point() +
ggplot(data = VAF_stat, mapping = aes(x = age, y = SBS13))+geom_point(colour = 'black',fill='black') +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  # geom_line(aes(y=fit_sbs_count_adj,x=age), size=0.8) +
  labs(y='SBS13',x="Age (yrs)")

#-------------SBS88-----------
lmm1 <- lme(SBS88 ~ condition + age+ region, random= ~1|patient,method = "REML",data=na.omit(VAF_stat))
summary(lmm1)
tmp<-summary(lmm1)[["tTable"]]
p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
tmp[,'p-value']<-p.adjust(summary(lmm1)[["tTable"]][,'p-value'],method='fdr')
VAF_stat$fit_SBS88<-predict(lmm1)


# ggplot(data = VAF_stat, mapping = aes(x = age, y = SBS13, colour = region,fill=region))+geom_point() +
ggplot(data = VAF_stat, mapping = aes(x = age, y = SBS88))+geom_point(colour = 'black',fill='black') +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  # geom_line(aes(y=fit_sbs_count_adj,x=age), size=0.8) +
  labs(y='SBS88',x="Age (yrs)")




library(ggplot2)
library(reshape2)
densFindPeak <- function(x){
  td <- density(x)
  maxDens <- which.max(td$y)
  list(x=td$x[maxDens],y=td$y[maxDens])
}
densFindPeak(VAF_stat$SBS2)
densFindPeak(VAF_stat$SBS13)

densFindPeak(VAF_stat$SBS13+VAF_stat$SBS2)
td <- density(VAF_stat$SBS2)
p<-ggplot(VAF_stat, aes(x = (SBS2))) + geom_density(aes(group=factor("SBS2")), fill="grey", size=1, alpha=.4)+ggtitle("APOBEC burden")+ylab('Density')+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  labs(x="SBS2 mutation count")

#-------------------------PD43402-------------------
snp_to_branch<-read.table('/Users/yw2/Documents/FARM/yw2/data/small_bowel/filter/PD43402/PD43402_snp_assigned_to_branches.txt',header=T)
lo001_snp<-snp_to_branch$Pos[snp_to_branch$SampleID == 'PD43402_6']
lo002_snp<-snp_to_branch$Pos[snp_to_branch$SampleID == 'PD43402_7']
shared_snp<-snp_to_branch$Pos[snp_to_branch$SampleID == 'PD43402_19']


path="PD43402/PD43402.snp.tsv"
patient_vcf<-read.table(paste0('/Users/yw2/Documents/FARM/yw2/data/small_bowel/filter/',path),comment.char = "",header=T, stringsAsFactors = F)
sample_name<-unlist(strsplit(colnames(patient_vcf[grep('MTR',colnames(patient_vcf))]),split='_MTR'))

all_VAF<-patient_vcf[,colnames(patient_vcf[grep('MTR',colnames(patient_vcf))])]/patient_vcf[,colnames(patient_vcf[grep('DEP',colnames(patient_vcf))])]
all_VAF<-cbind(patient_vcf[,3:6],all_VAF[,c('PD43402d_lo0001_MTR',"PD43402d_lo0002_MTR")])
all_VAF<-all_VAF[all_VAF$Pos %in% shared_snp,]
plot(all_VAF$PD43402d_lo0001_MTR,all_VAF$PD43402d_lo0002_MTR,xlab='PD43402d_lo0001',ylab='PD43402d_lo0002')


