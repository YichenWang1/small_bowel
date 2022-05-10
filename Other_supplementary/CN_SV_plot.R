#Plot copy number landscape with structural variants
#Tim Coorens, Sep 2019

options(stringsAsFactors = F)
data=read.table('./stat_summary.txt',header=T,stringsAsFactors = F)


for (sample_idx in 1:nrow(data)){
  sample=data$sample[sample_idx]
  project=data$project[sample_idx]

  brass=try(read.table(paste0("/nfs/cancer_ref01/nst_links/live/",project,"/",sample,"/",sample,".brass.annot.bedpe.gz"),sep="\t"))

  if(inherits(brass, "try-error")){
    brass=NA
  }
  else{
    brass=brass[grepl("score",brass$V23)|brass$V16>=15,] #This filters based on score, but also allows SVs to pass if having more than 15 supporting reads (arbitrary number)
    brass=brass[!grepl(",",brass$V11),] #Filter out anything called in the normal too
    brass=brass[brass$V13==-1|brass$V13>10^6,] #Filter out intrachromosomal SVs smaller than 1MB
  }
  

  # tum_cn=try(read.csv(paste0("/nfs/cancer_ref01/nst_links/live/",project,"/",sample,"/",sample,".ascat_ngs.summary.csv"),header=F))
  tum_cn=try(read.csv(paste0("/nfs/cancer_ref01/nst_links/live/",project,"/",sample,"/",sample,".battenberg.summary.csv"),header=F))
  
  
  
  if(inherits(tum_cn, "try-error")){
    next
  }
  CNVs=tum_cn[!(tum_cn$V5==tum_cn$V7&tum_cn$V6==tum_cn$V8),-1]
  
  colnames(CNVs)=c("Chr","Start","End","N1","N2","Maternal","Paternal")
  CNVs$Maternal=CNVs$Maternal-CNVs$Paternal
  
  chr_length=read.table("/lustre/scratch117/casm/team267/yw2/data/public/hg19_chr_length.txt",header=F) 
  pos=c(0,cumsum(as.numeric(chr_length$V2)))
  pos2=pos
  pos=pos[-25]
  
  pdf(paste0(sample,"_copynumber_plot.pdf"),width=15,height=5)
  plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
       ylab = "",xlim=c(-max(pos)*0.1,max(pos)*1.1),ylim=c(-1,6))
  
  rect(xleft=-0.01*max(pos),xright=max(pos)*1.01,ybottom=-0.2,ytop=4.4,lwd=1.2)
  paternal_cnv=c(rep(1,22),0)
  maternal_cnv=c(rep(1,22),1)
  for (i in 1:23){
    if(any(CNVs$Chr==i)){
      select=which(CNVs$Chr==i)
      for (s in 1:length(select)){
        n=select[s]
        rect(xleft=pos[i]+CNVs$Start[n],
             xright=pos[i]+CNVs$End[n],
             ybottom=CNVs$Paternal[n]-0.2,
             ytop=CNVs$Paternal[n],border=NA,col="steelblue")
        rect(xleft=pos[i]+CNVs$Start[n],
             xright=pos[i]+CNVs$End[n],
             ybottom=CNVs$Maternal[n],
             ytop=CNVs$Maternal[n]+0.2,border=NA,col="firebrick")
        if(CNVs$Start[n]>1&s==1){
          rect(xleft=pos[i],xright=pos[i]+CNVs$Start[n],ybottom=paternal_cnv[i]-0.2,ytop=paternal_cnv[i],border=NA,col="steelblue")
          rect(xleft=pos[i],xright=pos[i]+CNVs$Start[n],ybottom=maternal_cnv[i],ytop=maternal_cnv[i]+0.2,border=NA,col="firebrick")
        }
        if(s!=1){
          if(CNVs$Start[n]>CNVs$End[n-1]){
            rect(xleft=pos[i]+CNVs$End[n-1],xright=pos[i]+CNVs$Start[n],ybottom=paternal_cnv[i]-0.2,ytop=paternal_cnv[i],border=NA,col="steelblue")
            rect(xleft=pos[i]+CNVs$End[n-1],xright=pos[i]+CNVs$Start[n],ybottom=maternal_cnv[i],ytop=maternal_cnv[i]+0.2,border=NA,col="firebrick")
          }
        }
        if(CNVs$End[n]<chr_length$V2[i]&s==length(select)){
          rect(xleft=pos[i]+CNVs$End[n],xright=pos[i+1],ybottom=paternal_cnv[i]-0.2,ytop=paternal_cnv[i],border=NA,col="steelblue")
          rect(xleft=pos[i]+CNVs$End[n],xright=pos[i+1],ybottom=maternal_cnv[i],ytop=maternal_cnv[i]+0.2,border=NA,col="firebrick")
        }
      }
      
    }else{
      rect(xleft=pos[i],xright=pos[i+1],ybottom=paternal_cnv[i]-0.2,ytop=paternal_cnv[i],border=NA,col="steelblue")
      rect(xleft=pos[i],xright=pos[i+1],ybottom=maternal_cnv[i],ytop=maternal_cnv[i]+0.2,border=NA,col="firebrick")
    }
  }
  segments(y0=0:4,x0=-0.01*max(pos),x1=max(pos)*1.01,col="grey80")
  segments(y0=-0.2,y1=4.4,x0=pos,col="grey80")
  midpoints=(pos[-1]+pos[-length(pos)])/2
  text(x = midpoints,y=4.2,labels=c(1:22,"X"))
  #text(x=max(pos)/2,y=4.7,label=sample,cex=1.5)
  text(y=0:3,labels = 0:3,cex=1.2,x=-0.02*max(pos))
  
  colours=c("darkslateblue","darkorange4","chartreuse2","black")
  names(colours)=c("deletion","tandem-duplication","inversion","translocation")

  
  if (!is.na(brass) && nrow(brass)>0){
    brass$V1[brass$V1=="X"]=23
    brass$V4[brass$V4=="X"]=23
    brass$V1[brass$V1=="Y"]=24
    brass$V4[brass$V4=="Y"]=24
    for (i in 1:nrow(brass)){
      start=1
      #start=min(CNVs$Start[CNVs$Chr==brass$V1[i]])
      #if(start==Inf) start=1
      length=pos2[as.numeric(brass$V1[i])+1]-pos2[as.numeric(brass$V1[i])]
      ratio=1-start/length
      x_left=pos2[as.numeric(brass$V1[i])]+(brass$V2[i]-start)/ratio
      
      #start=min(CNVs$Start[CNVs$Chr==brass$V4[i]])
      #if(start==Inf) start=1
      length=pos2[as.numeric(brass$V4[i])+1]-pos2[as.numeric(brass$V4[i])]
      ratio=1-start/length
      x_right=pos2[as.numeric(brass$V4[i])]+(brass$V5[i]-start)/ratio
      
      segments(x0=x_left,y0=0,y1=2.5,col=colours[brass$V12[i]])
      segments(x0=x_right,y0=0,y1=2.5,col=colours[brass$V12[i]])
      x_mid = (x_left + x_right)/2  # Center of arc
      xr = x_right - x_mid
      print(sample)
      x_points = seq(x_left, x_right, length.out = 200)
      y_sign=1
      if(brass$V12[i]=="deletion") y_sign=-1
      y_points = 2.5 + y_sign*0.5 * sqrt(pmax(0,1-((x_points-x_mid)/xr)^2))
      lines(x_points,y_points,col=colours[brass$V12[i]])
    }
  }
  dev.off()
}