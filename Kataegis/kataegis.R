if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("MutationalPatterns")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library("GenomicRanges")
library("Rsamtools")
library("MASS")
library('vcfR')
library('ggplot2')
setwd('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/small_bowel/data')
options(stringsAsFactors = F)


# modified from package MutationalPatterns
plot_rainfall <- function (vcf, chromosomes, title = "", colors, cex = 2.5, cex_text = 3, 
                           ylim = 1e+08) 
{
  if (length(colors) != 6) 
    stop("colors vector length not 6")
  names(colors) = c("C>A","C>G","C>T","T>A","T>C","T>G")
  chr_length = seqlengths(vcf)
  chr_length = chr_length[names(chr_length) %in% chromosomes]
  chr_cum = c(0, cumsum(as.numeric(chr_length)))
  names(chr_cum) = names(chr_length)
  labels = gsub("chr", "", names(chr_length))
  m = c()
  for (i in 2:length(chr_cum)) m = c(m, (chr_cum[i - 1] + chr_cum[i])/2)
  type = loc = dist = chrom = c()
  for (i in 1:length(chromosomes)) {
    chr_subset = vcf[seqnames(vcf) == chromosomes[i]]
    n = length(chr_subset@seqnames)
    if (n <= 1) {
      next
    }
    type = c(type, mut_type(chr_subset)[-1])
    loc = c(loc, (start(chr_subset) + chr_cum[i])[-1])
    dist = c(dist, diff(start(chr_subset)))
    chrom = c(chrom, rep(chromosomes[i], n - 1))
  }
  data = data.frame(type = as.factor(type), location = loc, distance = dist, 
                    chromosome = chrom)
  colors = colors[levels(data$type)]
  location = NULL
  if (nrow(data)==0){
    return()
  }
  plot = ggplot(data, aes(x = location, y = distance)) + geom_point(aes(colour = factor(type)), 
                                                                    cex = cex) + geom_vline(xintercept = as.vector(chr_cum), 
                                                                                            linetype = "dotted") + annotate("text", x = m, y = ylim, 
                                                                                                                            label = labels, cex = cex_text) + xlab("Genomic Location") + 
    ylab("Genomic Distance") + scale_y_log10() + scale_colour_manual(values = colors) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(chr_cum))) +
    ggtitle(title) + theme_bw() + theme(legend.position = "bottom",
                                        legend.title = element_blank(), legend.key = element_blank(),
                                        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
                                        axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    guides(colour = guide_legend(nrow = 1))
  return(plot)
}


# function to identify kataegis by a negative binomial test
find_kataegis <- function (vcf,sample_name,outpath="./rainfall_plot/") 
{
  chromosomes <- seqnames(get(ref_genome))[1:23]
  chr_length = seqlengths(vcf)
  chr_length = chr_length[names(chr_length) %in% chromosomes]
  chr_cum = c(0, cumsum(as.numeric(chr_length)))
  names(chr_cum) = names(chr_length)
  prob  = as.numeric(length(vcf)/chr_cum[length(chr_cum)])
  
  type = loc = dist = chrom = c()
  data = data.frame(matrix(ncol =3, nrow = 0))
  for (i in 1:length(chromosomes)) {
    chr_subset = vcf[seqnames(vcf) == chromosomes[i]]
    n = length(chr_subset@seqnames)
    if (n <= 1) {
      next
    }
    type = mut_type(chr_subset)
    loc = start(chr_subset)
    chrom =  rep(chromosomes[i], n)
    dist=c(0,diff(loc))
    data_chr = data.frame(type = as.factor(type), location = loc, dist = dist,
                          chromosome = chrom)
    data_chr = data_chr[data_chr$dist>10 | data_chr$dist==0,]
    idx=1
    while (idx<length(data_chr$location)){
      start = idx
      end = start
      size = 0
      the_next = data_chr$location[idx+1]
      
      while (the_next<(data_chr$location[idx]+1e4) && idx<length(data_chr$location)){
        size = size + 1 # No.of successes
        end = idx + 1
        idx = idx + 1
        the_next = data_chr$location[idx+1]
      }
      if(end <= start){
        idx = idx + 1
        next
      }
      else{
        k = (0:(data_chr$location[end]-data_chr$location[start]-size)) # No.of failures
        pvalue = sum(dnbinom(k, size, prob))
        
        data_add = data_chr[start:end,]
        data_add$pvalue = pvalue
        data = rbind(data,data_add)
        idx = idx + 1
      }
    }
  }
  data$qvalue_bonferroni = p.adjust(data$pvalue,method = "bonferroni")
  data$qvalue_fdr = p.adjust(data$pvalue,method = "fdr")
  
  if (nrow(data)>0){
    trinuc_ref = as.vector(scanFa(genomeFile, GRanges(unlist(strsplit(data$chromosome,split='chr'))[c(F,T)], IRanges(data$location-1, data$location+1))))
    ntcomp = c(T="A",G="C",C="G",A="T")
    trinuc_ref_py = trinuc_ref
    for (j in 1:length(trinuc_ref)) {
      if (substr(trinuc_ref[j],2,2) %in% c("A","G")) { # Purine base
        trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(trinuc_ref[j],split="")[[1]])],collapse="")
      }
    }
    dinuc_ref_py <- substr(trinuc_ref_py,1,2)
    
    data$trinuc_ref_py = trinuc_ref_py
    data$APOBEC_motif = ((dinuc_ref_py =='TC') & (data$type =='C>T' | data$type =='C>G' | data$type =='C>A'))
    data$sample = sample_name
  }
  write.table(data, paste0(outpath,sample_name,"_kataegis.txt"),quote=F, row.names = F, sep = '\t')
}
  
## --------------- Identify kataegis and make rainfall plots---------------
# Load vcf files
vcf_files <- list.files(path="./snp",
                        pattern = "*_final_snp.vcf$", full.names = TRUE)
genomeFile = "/Users/yw2/OneDrive - University of Cambridge/PhD/data/public/genome.fa"

for (i in 1:length(vcf_files)){
  sample_name <- grep('PD',unlist(strsplit(grep('PD',unlist(strsplit(vcf_files[i],split='/')), value=T),split='_final_snp')),value=T) #c("PD43953b_lo0009")
  
  grl = try(read_vcfs_as_granges(vcf_files[i], sample_name, ref_genome))
  
  if (class(grl) != "CompressedGRangesList"){
    next
  }

  # Define autosomal chromosomes
  chromosomes <- seqnames(get(ref_genome))[1:23]
  find_kataegis(grl[[1]],sample_name)
  # Make a rainfall plot
  pdf(paste0("./rainfall_plot/",sample_name,"_rainfall.pdf"),width=10, height=2.5)
  print(plot_rainfall(grl[[1]],title = names(grl),chromosomes = chromosomes,  colors = c("dodgerblue","black","red","grey70","olivedrab3","plum2"), cex = 1.5, ylim = 1e+09))
  dev.off()
}

# Concatenate all individual *kataegis.txt into a single all.txt for multiple test correction
tmp<-read.table("./rainfall_plot/all.txt",header=T,sep='\t')
tmp$qvalue_bonferroni = p.adjust(tmp$pvalue,method = "bonferroni")
tmp$qvalue_fdr = p.adjust(tmp$pvalue,method = "fdr")


# Final list of kataegis
filtered <- tmp[tmp$qvalue_bonferroni < 0.0001,]
length(unique(filtered$pvalue))
table(filtered$APOBEC_motif) #FALSE=50   TRUE=402

for (sample in unique(filtered$sample)){
  sample_file = filtered[which(filtered$sample==sample),]
  write.table(sample_file,paste0("./rainfall_plot/final/",sample,"_kataegis.txt"),quote=F,sep='\t')
}
write.table(filtered,"./rainfall_plot/final/kataegis.txt",quote=F,sep='\t',row.names = F)


#------------APOBEC positive crypts with kataegis--------------
exposure_matrix_crypts = read.table("./data/sigs/exposure_matrix_crypts.txt",header=T,stringsAsFactors = F)
APOBEC_postive <- rownames(exposure_matrix_crypts)[exposure_matrix_crypts$SBS2+exposure_matrix_crypts$SBS13>0.05]
sum(APOBEC_postive %in% filtered$sample)
sum(APOBEC_postive %in% filtered$sample)/length(APOBEC_postive)
chisq.test(rbind(c(sum(APOBEC_postive %in% filtered$sample),length(unique(APOBEC_postive))-sum(APOBEC_postive %in% filtered$sample)), c(sum(!unique(filtered$sample) %in% APOBEC_postive),length(data$sample)-length(APOBEC_postive)-20)))


sum(filtered$sample %in% APOBEC_postive)
sum(!filtered$sample %in% APOBEC_postive)

df=data.frame(matrix(ncol = 2, nrow = 2))
colnames(df) <- c('SBS2&13', 'kataegis')
df$SBS2_13 = c('SBS2&13 positive','SBS2&13 negative')
df$kataegis = c(sum(filtered$sample %in% APOBEC_postive),sum(!filtered$sample %in% APOBEC_postive))

sum(filtered$sample %in% APOBEC_postive)/sum(df$kataegis)

df$SBS2_13 <-factor(df$SBS2_13,levels = c('SBS2&13 positive','SBS2&13 negative'))

# visualisation
ggplot(df) +
  geom_bar(aes(x = SBS2_13, y = kataegis,fill = c('kataegis','kataegis')),
           stat = "identity") +scale_fill_brewer(palette="Paired")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none",legend.title=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  labs(x = NULL, y = "Number of Mutations")

sum(data$sbs_count[data$sample %in% APOBEC_postive])

chisq.test(rbind(c(sum(filtered$sample %in% APOBEC_postive),sum(!filtered$sample %in% APOBEC_postive)), c(sum(data$sbs_count[data$sample %in% APOBEC_postive])-355,sum(data$sbs_count[!data$sample %in% APOBEC_postive])-97)))

#------------kataegis in APOBEC motif--------------
df=data.frame(matrix(ncol = 2, nrow = 2))
colnames(df) <- c('region', 'kataegis')
df$region = c('APOBEC motif','Others')
df$kataegis = c(402,50) #table(filtered$APOBEC_motif)
402/452
df$region<-factor(df$region,levels = c('APOBEC motif','Others'))

# visualisation
ggplot(df) +
  geom_bar(aes(x = region, y = kataegis,fill = c('kataegis','kataegis')),
           stat = "identity") +scale_fill_brewer(palette="Paired")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none",legend.title=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  labs(x = NULL, y = "Number of Mutations")


sbs<-read.table('./sbs/sbs_on_branch.txt',check.names = F)
sbs=sbs[!(rownames(sbs) %in% c('PD43851_1','PD43851_4','PD43851_14')),]

sum(sbs[,grep('C>[ATG],T-',colnames(sbs))])-402
sum(sbs[,-grep('C>[ATG],T-',colnames(sbs))])-50
chisq.test(rbind(c(402, 50), c(sum(input_for_hdp[,grep('C>[ATG],T-',colnames(input_for_hdp))])-402, sum(input_for_hdp[,-grep('C>[ATG],T-',colnames(input_for_hdp))])-50)))

