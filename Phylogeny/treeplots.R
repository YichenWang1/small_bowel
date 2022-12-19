# From Tim Coorens

library(ape)
library(ggplot2)
library(ggtree)
 
source("/lustre/scratch117/casm/team268/tc16/Scripts/TreeAssign/treemut.R")
 
options(stringsAsFactors = F)
patient=commandArgs(T)[1]
opt=commandArgs(T)[2]
NR_flt=read.table(paste0(patient,"/",opt,"_NR_filtered_all.txt"))
NV_flt=read.table(paste0(patient,"/",opt,"_NV_filtered_all.txt"))
tree=read.tree(paste0(patient,"/",opt,"_for_MPBoot.fa.treefile"))
tree$edge.length=rep(1,nrow(tree$edge))
tree=drop.tip(tree,"Ancestral")
NR_flt = as.matrix(NR_flt[,tree$tip.label])
NV_flt = as.matrix(NV_flt[,tree$tip.label])
 
XY_chromosomal = grepl("X|Y",rownames(NR_flt))
autosomal = !XY_chromosomal
xy_depth=mean(rowMeans(NR_flt[XY_chromosomal,]))
#xy_depth=mean(NR_flt[XY_chromosomal,])
autosomal_depth=mean(rowMeans(NR_flt[autosomal,]))
gender='male'
if(xy_depth>0.8*autosomal_depth) gender='female'
 
NR_flt_nonzero=NR_flt
NR_flt_nonzero[NR_flt_nonzero==0]=1
 
genotype_bin=as.matrix(NV_flt/NR_flt_nonzero)
if(gender=="male"){
  genotype_bin[autosomal,][genotype_bin[autosomal,]<0.1]=0
  genotype_bin[autosomal,][genotype_bin[autosomal,]>=0.3]=1
  genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]<0.2]=0
  genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]>=0.6]=1
  genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
}
if(gender=="female"){
  genotype_bin[genotype_bin<0.1]=0
  genotype_bin[genotype_bin>=0.3]=1
  genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
}
present_vars=rowSums(genotype_bin>0)>0
NR_flt=NR_flt[present_vars,]
NV_flt=NV_flt[present_vars,]
df=reconstruct_genotype_summary(tree)
res=assign_to_tree(df=df,
                   mtr=NV_flt,
                   dep=NR_flt)
 
saveRDS(res,paste0(patient,"/",opt,"_assigned_to_tree.Rdata"))
 
tree2=tree
edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<0.01])
edge_length = rep(0,nrow(tree$edge))
names(edge_length)=1:nrow(tree$edge)
edge_length[names(edge_length_nonzero)]=edge_length_nonzero
tree2$edge.length=as.numeric(edge_length)
write.tree(tree2, paste0(patient,"/",patient,"_",opt,"_tree_with_branch_length.tree"))
 
tree2 = as.polytomy(tree2, feature='branch.length', fun=function(x) as.numeric(x)==0)
 
write.tree(tree2, paste0(patient,"/",patient,"_",opt,"_tree_with_branch_length_polytomised.tree"))
 
p=ggtree(tree2) + geom_tiplab(aes(x=branch),vjust=-0.3)+theme_tree2()+xlim(0,max(fortify(tree2)$x)*1.3)
pdf(paste0(patient,"/",patient,"_",opt,"_tree_with_branch_length.pdf"))
print(p)
dev.off()
 
tree_collapsed=tree2
tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
pdf(paste0(patient,"/",patient,"_",opt,"_tree_with_equal_branch_length.pdf"))
ggtree(tree_collapsed) + geom_tiplab(aes(x=branch),vjust=-0.3)+xlim(0,max(fortify(tree_collapsed)$x)*1.3)
dev.off()
 
Mutations_per_branch=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(NR_flt),split="_")),byrow = T))
colnames(Mutations_per_branch)=c("Chr","Pos","Ref","Alt")
Mutations_per_branch$Branch = tree$edge[res$summary$edge_ml,2]
Mutations_per_branch=Mutations_per_branch[res$summary$p_else_where<0.01,]
Mutations_per_branch$Patient = patient
Mutations_per_branch$SampleID = paste(patient,Mutations_per_branch$Branch,sep="_")
write.table(Mutations_per_branch,paste0(patient,"/",patient,"_",opt,"_assigned_to_branches.txt"),quote=F,row.names=F,sep="\t")

