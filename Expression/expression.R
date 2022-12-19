library(SeuratDisk)
library(Seurat)
require(foreign)
require(ggplot2)
require(MASS)

#------APOBEC Expression - Gut Cell Survey - ----------
# Dataset downloaded from https://cellgeni.cog.sanger.ac.uk/gutcellatlas/epi_raw_counts02.h5ad
# Convert("epi_raw_counts02.h5ad", "h5seurat",overwrite =T)

# Load data
epi <- LoadH5Seurat("epi_raw_counts02.h5seurat")

LargeInt.subset.all <- subset(x = epi, subset = Region == "LargeInt")
SmallInt.subset.all <- subset(x = epi, subset = Region == "SmallInt")
rm(epi)
gc()

# LargeInt.subset.all <- LoadH5Seurat("/Users/yw2/Documents/FARM/public/scRNAseq/gut/LargeInt_raw_counts02.h5seurat")
# SmallInt.subset.all <- LoadH5Seurat("/Users/yw2/Documents/FARM/public/scRNAseq/gut/SmallInt_raw_counts02.h5seurat")
hist(SmallInt.subset.all@meta.data[["nCount_RNA"]], xlim = c(0,80000),xlab = "nCount_RNA",main='Small Intestine Epithelium',breaks=40)
median(SmallInt.subset.all@meta.data[["nCount_RNA"]])
hist(SmallInt.subset.all@meta.data[["nFeature_RNA"]],xlab = "nFeature_RNA",main='Small Intestine Epithelium')
median(SmallInt.subset.all@meta.data[["nFeature_RNA"]])

hist(LargeInt.subset.all@meta.data[["nCount_RNA"]], xlim = c(0,80000),xlab = "nCount_RNA",main='Large Intestine Epithelium',breaks=40)
median(LargeInt.subset.all@meta.data[["nCount_RNA"]])
hist(LargeInt.subset.all@meta.data[["nFeature_RNA"]],xlab = "nFeature_RNA",main='Large Intestine Epithelium')
median(LargeInt.subset.all@meta.data[["nFeature_RNA"]])

# Normalize data
SmallInt.subset.all<- NormalizeData(object = SmallInt.subset.all, normalization.method = "RC", 
                                    scale.factor = 10000)
LargeInt.subset.all<- NormalizeData(object = LargeInt.subset.all, normalization.method = "RC", 
                                    scale.factor = 10000)
# Stem cell subsets
LargeInt.subset <- subset(x = LargeInt.subset.all, subset = annotation == "Stem cells")
SmallInt.subset <- subset(x = SmallInt.subset.all, subset = annotation == "Stem cells")

hist(SmallInt.subset@meta.data[["nCount_RNA"]], xlim = c(0,80000),xlab = "nCount_RNA",main='Small Intestine Epithelium',breaks=20)
median(SmallInt.subset@meta.data[["nCount_RNA"]])
hist(SmallInt.subset@meta.data[["nFeature_RNA"]],xlab = "nFeature_RNA",main='Small Intestine Epithelium')
median(SmallInt.subset@meta.data[["nFeature_RNA"]])

hist(LargeInt.subset@meta.data[["nCount_RNA"]], xlim = c(0,80000),xlab = "nCount_RNA",main='Large Intestine Epithelium',breaks=20)
median(LargeInt.subset@meta.data[["nCount_RNA"]])
hist(LargeInt.subset@meta.data[["nFeature_RNA"]],xlab = "nFeature_RNA",main='Large Intestine Epithelium')
median(LargeInt.subset@meta.data[["nFeature_RNA"]])


# APOBEC expression in epithelium. Change the gene names as needed.
APOBEC_LargeInt = LargeInt.subset.all@assays$RNA[c('APOBEC1','APOBEC3A','APOBEC3B'),]
APOBEC_SmallInt = SmallInt.subset.all@assays$RNA[c('APOBEC1','APOBEC3A','APOBEC3B'),]

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(LargeInt.subset.all)[2]),rep('SmallInt',dim(SmallInt.subset.all)[2]))
df$nFeature=c(LargeInt.subset.all@meta.data$nFeature_RNA,SmallInt.subset.all@meta.data$nFeature_RNA)
df$nCount=c(LargeInt.subset.all@meta.data$nCount_RNA,SmallInt.subset.all@meta.data$nCount_RNA)
df$region=as.factor(df$region)

sum(df$APOBEC1[which(df$region=='SmallInt')])
sum(df$region=='SmallInt')
sum(df$APOBEC1[which(df$region=='SmallInt')])/sum(df$region=='SmallInt')

sum(df$APOBEC1[which(df$region=='LargeInt')])
sum(df$region=='LargeInt')
sum(df$APOBEC1[which(df$region=='LargeInt')])/sum(df$region=='LargeInt')

(sum(df$APOBEC1[which(df$region=='SmallInt')])/sum(df$region=='SmallInt')) / (sum(df$APOBEC1[which(df$region=='LargeInt')])/sum(df$region=='LargeInt')) 

# Whether APOBEC expression level is different in small vs. large intestine epithelium. Change the gene names as needed
APOBEC_LargeInt = LargeInt.subset.all@assays$RNA@counts[c('APOBEC1','APOBEC3A','APOBEC3B'),]
APOBEC_SmallInt = SmallInt.subset.all@assays$RNA@counts[c('APOBEC1','APOBEC3A','APOBEC3B'),]

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(LargeInt.subset.all)[2]),rep('SmallInt',dim(SmallInt.subset.all)[2]))
df$nFeature=c(LargeInt.subset.all@meta.data$nFeature_RNA,SmallInt.subset.all@meta.data$nFeature_RNA)
df$nCount=c(LargeInt.subset.all@meta.data$nCount_RNA,SmallInt.subset.all@meta.data$nCount_RNA)
df$region=as.factor(df$region)

summary(m1 <- glm.nb(APOBEC1 ~ nFeature + nCount + region +APOBEC3A +APOBEC3B, data = df), maxit = 500)[["coefficients"]]
m2 <- update(m1, . ~ . - region)
anova(m1, m2)


# APOBEC expression in stem cells. Change the gene names as needed.
APOBEC_LargeInt = LargeInt.subset@assays$RNA[c('APOBEC1','APOBEC3A','APOBEC3B'),]
APOBEC_SmallInt = SmallInt.subset@assays$RNA[c('APOBEC1','APOBEC3A','APOBEC3B'),]

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(LargeInt.subset)[2]),rep('SmallInt',dim(SmallInt.subset)[2]))
df$nFeature=c(LargeInt.subset@meta.data$nFeature_RNA,SmallInt.subset@meta.data$nFeature_RNA)
df$nCount=c(LargeInt.subset@meta.data$nCount_RNA,SmallInt.subset@meta.data$nCount_RNA)
df$region=as.factor(df$region)


sum(df$APOBEC1[which(df$region=='SmallInt')])
sum(df$region=='SmallInt')
sum(df$APOBEC1[which(df$region=='SmallInt')])/sum(df$region=='SmallInt')

sum(df$APOBEC1[which(df$region=='LargeInt')])
sum(df$region=='LargeInt')
sum(df$APOBEC1[which(df$region=='LargeInt')])/sum(df$region=='LargeInt')

(sum(df$APOBEC1[which(df$region=='SmallInt')])/sum(df$region=='SmallInt')) / (sum(df$APOBEC1[which(df$region=='LargeInt')])/sum(df$region=='LargeInt')) 

# Whether APOBEC expression level is different in small vs. large intestine epithelial stem cells. Change the gene names as needed
APOBEC_LargeInt = LargeInt.subset@assays$RNA@counts[c('APOBEC1','APOBEC3A','APOBEC3B'),]
APOBEC_SmallInt = SmallInt.subset@assays$RNA@counts[c('APOBEC1','APOBEC3A','APOBEC3B'),]

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(LargeInt.subset)[2]),rep('SmallInt',dim(SmallInt.subset)[2]))
df$nFeature=c(LargeInt.subset@meta.data$nFeature_RNA,SmallInt.subset@meta.data$nFeature_RNA)
df$nCount=c(LargeInt.subset@meta.data$nCount_RNA,SmallInt.subset@meta.data$nCount_RNA)
df$region=as.factor(df$region)

summary(m1 <- glm.nb(APOBEC1 ~ nFeature + nCount + region +APOBEC3B +APOBEC3A, data = df))[["coefficients"]]
m2 <- update(m1, . ~ . - region)
anova(m1, m2)


#------APOBEC Expression - Human Protien Atlas Collection -----------
# Dataset downloaded from https://www.proteinatlas.org/download/rna_single_cell_read_count.tsv.zip
# Load data
HPA = read.table('/Users/yw2/Documents/FARM/public/scRNAseq/HPA/subset.tsv',sep='\t',header=T)
HPA$nCount = apply(HPA[,4:dim(HPA)[2]],1,sum)
rownames(HPA) = paste0(HPA$Tissue,'_',HPA$Cell)
HPA_meta = data.frame(HPA[1:3])
HPA_seurat <- CreateSeuratObject(counts = t(HPA[,4:(dim(HPA)[2]-1)]), meta.data = HPA_meta)

HPA_SmallInt.all <- subset(x = HPA_seurat, subset = Tissue == "Small intestine")
HPA_LargeInt.all <- subset(x = HPA_seurat, subset = Tissue == "Colon")

# Normalize data
HPA_SmallInt.all<- NormalizeData(object = HPA_SmallInt.all, normalization.method = "RC", 
                                 scale.factor = 10000)
HPA_LargeInt.all<- NormalizeData(object = HPA_LargeInt.all, normalization.method = "RC", 
                                 scale.factor = 10000)

# Undifferentiated cells
HPA_LargeInt.subset <- subset(x = HPA_LargeInt.all, subset = Cluster %in% c(2,8,13))
HPA_SmallInt.subset <- subset(x = HPA_SmallInt.all, subset = Cluster %in% c(3,7,8))


# APOBEC expression in epithelium. Change the gene names as needed.
APOBEC_LargeInt = HPA_LargeInt.all@assays$RNA[c('ENSG00000111701','ENSG00000128383','ENSG00000179750'),]
APOBEC_SmallInt = HPA_SmallInt.all@assays$RNA[c('ENSG00000111701','ENSG00000128383','ENSG00000179750'),]
rownames(APOBEC_LargeInt) = c('APOBEC1','APOBEC3A','APOBEC3B')
rownames(APOBEC_SmallInt) = c('APOBEC1','APOBEC3A','APOBEC3B')

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(HPA_LargeInt.all)[2]),rep('SmallInt',dim(HPA_SmallInt.all)[2]))
df$nFeature=c(HPA_LargeInt.all@meta.data$nFeature_RNA,HPA_SmallInt.all@meta.data$nFeature_RNA)
df$nCount=c(HPA_LargeInt.all@meta.data$nCount_RNA,HPA_SmallInt.all@meta.data$nCount_RNA)
df$region=as.factor(df$region)

sum(df$APOBEC3B[which(df$region=='SmallInt')])
sum(df$region=='SmallInt')
sum(df$APOBEC3B[which(df$region=='SmallInt')])/sum(df$region=='SmallInt')

sum(df$APOBEC3B[which(df$region=='LargeInt')])
sum(df$region=='LargeInt')
sum(df$APOBEC3B[which(df$region=='LargeInt')])/sum(df$region=='LargeInt')

(sum(df$APOBEC3B[which(df$region=='SmallInt')])/sum(df$region=='SmallInt')) / (sum(df$APOBEC3B[which(df$region=='LargeInt')])/sum(df$region=='LargeInt')) 



# Whether APOBEC expression level is different in small vs. large intestine epithelial cells. Change the gene names as needed
APOBEC_LargeInt = HPA_LargeInt.all@assays$RNA@counts[c('ENSG00000111701','ENSG00000128383','ENSG00000179750'),]
APOBEC_SmallInt = HPA_SmallInt.all@assays$RNA@counts[c('ENSG00000111701','ENSG00000128383','ENSG00000179750'),]
rownames(APOBEC_LargeInt) = c('APOBEC1','APOBEC3A','APOBEC3B')
rownames(APOBEC_SmallInt) = c('APOBEC1','APOBEC3A','APOBEC3B')

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(HPA_LargeInt.all)[2]),rep('SmallInt',dim(HPA_SmallInt.all)[2]))
df$nFeature=c(HPA_LargeInt.all@meta.data$nFeature_RNA,HPA_SmallInt.all@meta.data$nFeature_RNA)
df$nCount=c(HPA_LargeInt.all@meta.data$nCount_RNA,HPA_SmallInt.all@meta.data$nCount_RNA)
df$region=as.factor(df$region)

summary(m1 <- glm.nb(APOBEC1 ~ nFeature + nCount + region + APOBEC3A + APOBEC3B, data = df))[["coefficients"]]
m2 <- update(m1, . ~ . - region)
anova(m1, m2)

# APOBEC expression in undifferentiated cells in epithelial undifferentiated cells. Change the gene names as needed.
APOBEC_LargeInt = HPA_LargeInt.subset@assays$RNA[c('ENSG00000111701','ENSG00000128383','ENSG00000179750'),]
APOBEC_SmallInt = HPA_SmallInt.subset@assays$RNA[c('ENSG00000111701','ENSG00000128383','ENSG00000179750'),]
rownames(APOBEC_LargeInt) = c('APOBEC1','APOBEC3A','APOBEC3B')
rownames(APOBEC_SmallInt) = c('APOBEC1','APOBEC3A','APOBEC3B')

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(HPA_LargeInt.subset)[2]),rep('SmallInt',dim(HPA_SmallInt.subset)[2]))
df$nFeature=c(HPA_LargeInt.subset@meta.data$nFeature_RNA,HPA_SmallInt.subset@meta.data$nFeature_RNA)
df$nCount=c(HPA_LargeInt.subset@meta.data$nCount_RNA,HPA_SmallInt.subset@meta.data$nCount_RNA)
df$region=as.factor(df$region)

sum(df$APOBEC1[which(df$region=='SmallInt')])
sum(df$region=='SmallInt')
sum(df$APOBEC1[which(df$region=='SmallInt')])/sum(df$region=='SmallInt')

sum(df$APOBEC1[which(df$region=='LargeInt')])
sum(df$region=='LargeInt')
sum(df$APOBEC1[which(df$region=='LargeInt')])/sum(df$region=='LargeInt')

(sum(df$APOBEC1[which(df$region=='SmallInt')])/sum(df$region=='SmallInt')) / (sum(df$APOBEC1[which(df$region=='LargeInt')])/sum(df$region=='LargeInt')) 

# Whether APOBEC expression level is different in small vs. large intestine epithelial undifferentiated cells. Change the gene names as needed
APOBEC_LargeInt = HPA_LargeInt.subset@assays$RNA@counts[c('ENSG00000111701','ENSG00000128383','ENSG00000179750'),]
APOBEC_SmallInt = HPA_SmallInt.subset@assays$RNA@counts[c('ENSG00000111701','ENSG00000128383','ENSG00000179750'),]
rownames(APOBEC_LargeInt) = c('APOBEC1','APOBEC3A','APOBEC3B')
rownames(APOBEC_SmallInt) = c('APOBEC1','APOBEC3A','APOBEC3B')

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(HPA_LargeInt.subset)[2]),rep('SmallInt',dim(HPA_SmallInt.subset)[2]))
df$nFeature=c(HPA_LargeInt.subset@meta.data$nFeature_RNA,HPA_SmallInt.subset@meta.data$nFeature_RNA)
df$nCount=c(HPA_LargeInt.subset@meta.data$nCount_RNA,HPA_SmallInt.subset@meta.data$nCount_RNA)
df$region=as.factor(df$region)

summary(m1 <- glm.nb(APOBEC1 ~ nFeature + nCount + region +APOBEC3A +APOBEC3B, data = df))[["coefficients"]]
m2 <- update(m1, . ~ . - region)
anova(m1, m2)

#------APOBEC Expression - Tabula Sapiens -----------
# Dataset downloaded from https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5# Load data
# Load data
all<-readRDS("/Users/yw2/Documents/FARM/public/scRNAseq/various/epithelial.rds")
Key(all)
all@assays$RNA@key <- "rna_"
all@meta.data[["cell_type"]]

Large_Intestine_all <- subset(x = all, subset = tissue_in_publication == "Large_Intestine")
hist(Large_Intestine_all@meta.data[["nCount_RNA"]],xlab = "nCount_RNA",main='Large Intestine Epithelium',breaks = 800)
median(Large_Intestine_all@meta.data[["nCount_RNA"]])
hist(Large_Intestine_all@meta.data[["nFeature_RNA"]],xlab = "nFeature_RNA",main='Large Intestine Epithelium',breaks=6)
median(Large_Intestine_all@meta.data[["nFeature_RNA"]])

Small_Intestine_all <- subset(x = all, subset = tissue_in_publication == "Small_Intestine")
hist(Small_Intestine_all@meta.data[["nCount_RNA"]],xlab = "nCount_RNA",main='Small Intestine Epithelium',breaks = 800)
median(Small_Intestine_all@meta.data[["nCount_RNA"]])
hist(Small_Intestine_all@meta.data[["nFeature_RNA"]],xlab = "nFeature_RNA",main='Small Intestine Epithelium',breaks=6)
median(Small_Intestine_all@meta.data[["nFeature_RNA"]])

# Normalize data
Small_Intestine_all<- NormalizeData(object = Small_Intestine_all, normalization.method = "RC", 
                                    scale.factor = 10000)
Large_Intestine_all <- NormalizeData(object = Large_Intestine_all, normalization.method = "RC", 
                                     scale.factor = 10000)

# APOBEC expression in epithelium. Change the gene names as needed.
APOBEC_LargeInt = Large_Intestine_all@assays$RNA[Large_Intestine_all@assays[["RNA"]]@meta.features[["feature_name"]] %in% c('APOBEC1','APOBEC3A','APOBEC3B'),]
APOBEC_SmallInt = Small_Intestine_all@assays$RNA[Small_Intestine_all@assays[["RNA"]]@meta.features[["feature_name"]] %in% c('APOBEC1','APOBEC3A','APOBEC3B'),]
rownames(APOBEC_LargeInt) = c('APOBEC1','APOBEC3A','APOBEC3B')
rownames(APOBEC_SmallInt) = c('APOBEC1','APOBEC3A','APOBEC3B')

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(APOBEC_LargeInt)[2]),rep('SmallInt',dim(APOBEC_SmallInt )[2]))
df$nFeature=c(Large_Intestine_all@meta.data$nFeature_RNA,Small_Intestine_all@meta.data$nFeature_RNA)
df$nCount=c(Large_Intestine_all@meta.data$nCount_RNA,Small_Intestine_all@meta.data$nCount_RNA)
df$region=as.factor(df$region)
df$cell_type=c(as.character(Large_Intestine_all@meta.data$cell_type),as.character(Small_Intestine_all@meta.data$cell_type))
df$tissue_type=c(as.character(Large_Intestine_all@meta.data$cell_type),as.character(Small_Intestine_all@meta.data$cell_type))
df = df[grep('10X',rownames(df)),]

sum(df$APOBEC1[which(df$region=='SmallInt')])
sum(df$region=='SmallInt')
mean(df$APOBEC1[which(df$region=='SmallInt')])

sum(df$APOBEC1[which(df$region=='LargeInt')])
sum(df$region=='LargeInt')
mean(df$APOBEC1[which(df$region=='LargeInt')])

mean(df$APOBEC1[which(df$region=='SmallInt')]) / mean(df$APOBEC1[which(df$region=='LargeInt')])


# Whether APOBEC expression level is different in small vs. large intestine epithelium. Change the gene names as needed
APOBEC_LargeInt = Large_Intestine_all@assays$RNA@counts[Large_Intestine_all@assays[["RNA"]]@meta.features[["feature_name"]] %in% c('APOBEC1','APOBEC3A','APOBEC3B'),]
APOBEC_SmallInt = Small_Intestine_all@assays$RNA@counts[Small_Intestine_all@assays[["RNA"]]@meta.features[["feature_name"]] %in% c('APOBEC1','APOBEC3A','APOBEC3B'),]
rownames(APOBEC_LargeInt) = c('APOBEC1','APOBEC3A','APOBEC3B')
rownames(APOBEC_SmallInt) = c('APOBEC1','APOBEC3A','APOBEC3B')

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(APOBEC_LargeInt)[2]),rep('SmallInt',dim(APOBEC_SmallInt )[2]))
df$nFeature=c(Large_Intestine_all@meta.data$nFeature_RNA,Small_Intestine_all@meta.data$nFeature_RNA)
df$nCount=c(Large_Intestine_all@meta.data$nCount_RNA,Small_Intestine_all@meta.data$nCount_RNA)
df$region=as.factor(df$region)
df$cell_type=c(as.character(Large_Intestine_all@meta.data$cell_type),as.character(Small_Intestine_all@meta.data$cell_type))
df$tissue_type=c(as.character(Large_Intestine_all@meta.data$cell_type),as.character(Small_Intestine_all@meta.data$cell_type))
df = df[grep('10X',rownames(df)),]

summary(m1 <- glm.nb(APOBEC1 ~ nFeature + nCount + region +APOBEC3B +APOBEC3A, data = df, maxit = 500))[["coefficients"]]
m2 <- update(m1, . ~ . - region)
anova(m1, m2)


# intestinal crypt stem cell of large intestine/transit amplifying cell of colon/transit amplifying cell of small intestine/intestinal crypt stem cell
Large_Intestine_subset <- subset(x = Large_Intestine_all, subset = cell_type %in% c("intestinal crypt stem cell of large intestine","transit amplifying cell of colon"))
median(Large_Intestine_subset@meta.data[["nCount_RNA"]])
median(Large_Intestine_subset@meta.data[["nFeature_RNA"]])
dim(Large_Intestine_subset) 

Small_Intestine_subset <- subset(x = Small_Intestine_all, subset = cell_type %in% c("intestinal crypt stem cell","transit amplifying cell of small intestine"))
median(Small_Intestine_subset@meta.data[["nCount_RNA"]])
median(Small_Intestine_subset@meta.data[["nFeature_RNA"]])
dim(Small_Intestine_subset) 

# APOBEC expression in stem cells. Change the gene names as needed.
APOBEC_LargeInt = Large_Intestine_subset@assays$RNA[Large_Intestine_subset@assays[["RNA"]]@meta.features[["feature_name"]] %in% c('APOBEC1','APOBEC3A','APOBEC3B'),]
APOBEC_SmallInt = Small_Intestine_subset@assays$RNA[Small_Intestine_subset@assays[["RNA"]]@meta.features[["feature_name"]] %in% c('APOBEC1','APOBEC3A','APOBEC3B'),]
rownames(APOBEC_LargeInt) = c('APOBEC1','APOBEC3A','APOBEC3B')
rownames(APOBEC_SmallInt) = c('APOBEC1','APOBEC3A','APOBEC3B')

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(APOBEC_LargeInt)[2]),rep('SmallInt',dim(APOBEC_SmallInt )[2]))
df$nFeature=c(Large_Intestine_subset@meta.data$nFeature_RNA,Small_Intestine_subset@meta.data$nFeature_RNA)
df$nCount=c(Large_Intestine_subset@meta.data$nCount_RNA,Small_Intestine_subset@meta.data$nCount_RNA)
df$region=as.factor(df$region)
df = df[grep('10X',rownames(df)),]

sum(df$APOBEC1[which(df$region=='SmallInt')])
sum(df$region=='SmallInt')
sum(df$APOBEC1[which(df$region=='SmallInt')])/sum(df$region=='SmallInt')

sum(df$APOBEC1[which(df$region=='LargeInt')])
sum(df$region=='LargeInt')
sum(df$APOBEC1[which(df$region=='LargeInt')])/sum(df$region=='LargeInt')

(sum(df$APOBEC1[which(df$region=='SmallInt')])/sum(df$region=='SmallInt')) / (sum(df$APOBEC1[which(df$region=='LargeInt')])/sum(df$region=='LargeInt')) 

# Whether APOBEC expression level is different in small vs. large intestine epithelial stem/transit amplifying cells. Change the gene names as needed
APOBEC_LargeInt = Large_Intestine_subset@assays$RNA@counts[Large_Intestine_subset@assays[["RNA"]]@meta.features[["feature_name"]] %in% c('APOBEC1','APOBEC3A','APOBEC3B'),]
APOBEC_SmallInt = Small_Intestine_subset@assays$RNA@counts[Small_Intestine_subset@assays[["RNA"]]@meta.features[["feature_name"]] %in% c('APOBEC1','APOBEC3A','APOBEC3B'),]
rownames(APOBEC_LargeInt) = c('APOBEC1','APOBEC3A','APOBEC3B')
rownames(APOBEC_SmallInt) = c('APOBEC1','APOBEC3A','APOBEC3B')

df = data.frame(APOBEC1 = c(APOBEC_LargeInt['APOBEC1',],APOBEC_SmallInt['APOBEC1',]), APOBEC3A = c(APOBEC_LargeInt['APOBEC3A',],APOBEC_SmallInt['APOBEC3A',]), APOBEC3B = c(APOBEC_LargeInt['APOBEC3B',],APOBEC_SmallInt['APOBEC3B',]))
df$region=c(rep('LargeInt',dim(APOBEC_LargeInt)[2]),rep('SmallInt',dim(APOBEC_SmallInt )[2]))
df$nFeature=c(Large_Intestine_subset@meta.data$nFeature_RNA,Small_Intestine_subset@meta.data$nFeature_RNA)
df$nCount=c(Large_Intestine_subset@meta.data$nCount_RNA,Small_Intestine_subset@meta.data$nCount_RNA)
df$region=as.factor(df$region)
df = df[grep('10X',rownames(df)),]

summary(m1 <- glm.nb(APOBEC3B ~ nFeature + nCount +APOBEC1 +APOBEC3A+ region , data = df, maxit=1000))[["coefficients"]]
m2 <- update(m1, . ~ . - region)
anova(m1, m2)
