exposure_matrix = read.table('./data/signatures/exposure_matrix_branches.txt')
exposure_matrix = exposure_matrix[!(rownames(exposure_matrix) %in% c('PD43851_1','PD43851_4','PD43851_14')),]
input_for_hdp<-read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/small_bowel/data/input_for_hdp.txt',check.names = F)
input_for_hdp=input_for_hdp[apply(input_for_hdp,1,sum)>50,]
input_for_hdp=input_for_hdp[!(rownames(input_for_hdp) %in% c('PD43851_1','PD43851_4','PD43851_14')),]

motif_enrichment = exposure_matrix * apply(input_for_hdp,1,sum)
rtCa_all = read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/small_bowel/data/motif/P-MACD/work/res_rtCa/small_bowel_snp_MAF_sorted_rtCa_sum_all_fisher_Pcorr.txt', comment.char = '!',sep='\t',header=T,check.names = F)
rownames(rtCa_all) = rtCa_all$Sample

motif_enrichment$all_rtCa_toGT = rtCa_all[rownames(motif_enrichment),]$`rtCa_to_G+rtCa_to_T`
motif_enrichment$all_rtca = rtCa_all[rownames(motif_enrichment),]$`rtca+tgay`
motif_enrichment$all_CtoGT = rtCa_all[rownames(motif_enrichment),]$`[(C_to_G)+(C_to_T)]-[(rtCa_to_G)+(rtCa_to_T)]`+rtCa_all[rownames(motif_enrichment),]$`rtCa_to_G+rtCa_to_T`
motif_enrichment$all_c = rtCa_all[rownames(motif_enrichment),]$`c-rtca`+rtCa_all[rownames(motif_enrichment),]$`rtca+tgay`
motif_enrichment$rtCa_enrich = rtCa_all[rownames(motif_enrichment),]$APOBECrtCa_enrich

ytCa_all = read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/small_bowel/data/motif/P-MACD/work/res_ytCa/small_bowel_snp_MAF_sorted_ytCa_sum_all_fisher_Pcorr.txt', comment.char = '!',sep='\t',header=T,check.names = F)
rownames(ytCa_all) = ytCa_all$Sample
motif_enrichment$all_ytCa_toGT = ytCa_all[rownames(motif_enrichment),]$`ytCa_to_G+ytCa_to_T`
motif_enrichment$all_ytca = ytCa_all[rownames(motif_enrichment),]$`ytca+tgar`
motif_enrichment$ytCa_enrich = ytCa_all[rownames(motif_enrichment),]$APOBECytCa_enrich

tCa_all = read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/small_bowel/data/motif/P-MACD/work/res_tCa/small_bowel_snp_MAF_sorted_tCa_sum_all_fisher_Pcorr.txt', comment.char = '!',sep='\t',header=T,check.names = F)
rownames(tCa_all) = tCa_all$Sample
motif_enrichment$all_tCa_toGT = tCa_all[rownames(motif_enrichment),]$`tCa_to_G+tCa_to_T`
motif_enrichment$all_tca = tCa_all[rownames(motif_enrichment),]$`tca+tga`
motif_enrichment$tca_enrich = tCa_all[rownames(motif_enrichment),]$APOBECtCa_enrich


tCw_all = read.table('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/small_bowel/data/motif/P-MACD/work/res_tCw/small_bowel_snp_MAF_sorted_sum_all_fisher_Pcorr.txt', comment.char = '!',sep='\t',header=T,check.names = F)
rownames(tCw_all) =tCw_all$Sample
motif_enrichment$all_tCw_toGT = tCw_all[rownames(motif_enrichment),]$`tCw_to_G+tCw_to_T+revcomp`
motif_enrichment$all_tcw = tCw_all[rownames(motif_enrichment),]$`tcw+revcomp`
motif_enrichment$tcw_enrich = tCw_all[rownames(motif_enrichment),]$APOBEC_enrich

#--------------summary-----------------
APOBEC_positive_index = which((motif_enrichment$SBS2+motif_enrichment$SBS2) >0)

all_rtCa_toGT = sum(motif_enrichment$all_rtCa_toGT[APOBEC_positive_index])
all_rtca = sum(motif_enrichment$all_rtca[APOBEC_positive_index])
all_CtoGT = sum(motif_enrichment$all_CtoGT[APOBEC_positive_index])
all_c = sum(motif_enrichment$all_c[APOBEC_positive_index])
enrich_rtCa = (all_rtCa_toGT/all_rtca)/(all_CtoGT/all_c)
enrich_rtCa #0.9172801
fisher.test(matrix(data=c(all_rtCa_toGT,all_rtca,all_CtoGT-all_rtCa_toGT,all_c-all_rtca),nrow=2),alternative = 'greater')


all_ytCa_toGT = sum(motif_enrichment$all_ytCa_toGT[APOBEC_positive_index])
all_ytca = sum(motif_enrichment$all_ytca[APOBEC_positive_index])
enrich_ytCa = (all_ytCa_toGT/all_ytca)/(all_CtoGT/all_c)
enrich_ytCa #1.335129
fisher.test(matrix(data=c(all_ytCa_toGT,all_ytca,all_CtoGT-all_ytCa_toGT,all_c-all_ytca),nrow=2),alternative = 'greater')


all_tCa_toGT = sum(motif_enrichment$all_tCa_toGT[APOBEC_positive_index])
all_tca = sum(motif_enrichment$all_tca[APOBEC_positive_index])
enrich_tCa = (all_tCa_toGT/all_tca)/(all_CtoGT/all_c)
enrich_tCa
fisher.test(matrix(data=c(all_tCa_toGT,all_tca,all_CtoGT-all_tCa_toGT,all_c-all_tca),nrow=2),alternative = 'greater')


all_tCw_toGT = sum(motif_enrichment$all_tCw_toGT[APOBEC_positive_index])
all_tcw = sum(motif_enrichment$all_tcw[APOBEC_positive_index])
enrich_tCw = (all_tCw_toGT/all_tcw)/(all_CtoGT/all_c)
enrich_tCw
fisher.test(matrix(data=c(all_tCw_toGT,all_tcw,all_CtoGT-all_tCw_toGT,all_c-all_tcw),nrow=2),alternative = 'greater')


df_enrichment = data.frame(matrix(ncol=2,nrow = 4))
colnames(df_enrichment) = c('Fold_enrichment', 'Motif')
df_enrichment$Motif = c('TCW','TCA','YTCA','RTCA')
df_enrichment$Motif = factor(df_enrichment$Motif,levels = c('TCW','TCA','YTCA','RTCA'))
df_enrichment$Fold_enrichment = c(enrich_tCw,enrich_tCa,enrich_ytCa,enrich_rtCa)
library(ggplot2)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(4, "Set3"))

pdf('./fig/Extended_Fig7_extended_context/fold_enrichment.pdf',width = 5,height = 3)
ggplot(data=df_enrichment, aes(x=Motif, y=Fold_enrichment, fill=Motif )) +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  geom_bar(stat="identity")+ coord_flip() +geom_text(aes(label=round(Fold_enrichment,2)), hjust=-0.3, size=3.5)+scale_fill_brewer()+ylim(c(0,1.5))+labs(y='Fold enrichment',x="Motif")#+scale_fill_manual(values=getPalette(4))
dev.off()
plot(motif_enrichment$tcw_enrich,motif_enrichment$tca_enrich)

write.table(motif_enrichment,'./data/motif/motif_enrichment.csv', quote=F)