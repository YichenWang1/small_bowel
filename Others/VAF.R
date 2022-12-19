options(stringsAsFactors=F)
library(ggplot2)
data=read.csv('./data/Extended_Data_Table3_crypt_summary.csv',header=T,stringsAsFactors = F)

#---------extract VAF---------------
for (patient in unique(data$patient)){
  locus_list = data$sample[data$patient==patient]
  main_dir=paste0('/Users/yw2/Documents/FARM/yw2/small_bowel/filter/',patient)
  patient_vcf<-read.table(paste0(main_dir,"/",patient,".snp.tsv"),comment.char = "", header=T,stringsAsFactors = F)
  rownames(patient_vcf)=paste(patient_vcf$Chrom,patient_vcf$Pos,patient_vcf$Ref,patient_vcf$Alt, sep = '_')
  df= data.frame(matrix(nrow=0,ncol=4))
  for (locus in locus_list){
    if(locus %in% data$sample[data$ref!="Crypt"]){
      next
    }
    good_snp <-read.table(paste0("./data/vcf/",locus,"_final_snp.vcf"),header=F, stringsAsFactors = F)
    good_snp <- good_snp[!duplicated(good_snp),]
    rownames(good_snp)=paste(good_snp$V1,good_snp$V2,good_snp$V4,good_snp$V5, sep = '_')
    sample_vcf = patient_vcf[rownames(good_snp),]
    
    all_alt<-sample_vcf[,colnames(sample_vcf[grep('MTR',colnames(sample_vcf))])]
    all_dep<-sample_vcf[,colnames(sample_vcf[grep('DEP',colnames(sample_vcf))])]
    sample_alt<-all_alt[,colnames(all_alt)[grep(locus,colnames(all_alt))]]
    sample_dep<-all_dep[,colnames(all_dep)[grep(locus,colnames(all_dep))]]
    df_sample=data.frame(ALT=sample_alt,DEP=sample_dep)
    df_sample <- na.omit(df_sample[df_sample$ALT>=8,])
    df_sample$VAF = df_sample$ALT/df_sample$DEP
    df_sample$locus = locus
    df = rbind(df,df_sample)
    }
  df <-df[df$locus %in% data$sample[which(data$ref == 'Crypt')],]

  pdf(paste0(patient,'_VAF_plot.pdf'),height = 2, width = 2)
  print(ggplot(df, aes(x = VAF)) + geom_density(aes(group=factor(locus)), fill="grey", size=0.5, alpha=.4)+ggtitle(patient)+ylab('Density')+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black')))
  dev.off()
}
