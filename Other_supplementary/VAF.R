library(ggplot2)
library(nlme)
library(lme4)
VAF_stat <- data.frame()
#---------VAF plot---------------
for (path in vcf_file$V1){
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
  
colnames(VAF_stat)=c('sample','median_VAF')
