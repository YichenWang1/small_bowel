# From Henry Lee-Six
library(RColorBrewer)
library(plot3D)

# setwd("/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/small_bowel/data/stem_cell/observed_VAF")
data=read.csv('/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/code/small_bowel/data/Extended_Data_Table3_crypt_summary.csv',header=T,stringsAsFactors = F)
PD_ID = unique(data$patient)
PD_ID = unique(data$patient[data$condition=='Normal'])
# remove the oen pateint undergone long-term chemotherapy and children
PD_ID = PD_ID[-which(PD_ID=='PD43853')]
PD_ID = PD_ID[-which(substr(PD_ID,1,5)=='PD439')]

meanyears_to_MRCA = data.frame(patient=PD_ID,median_meanyears_to_MRCA=NA,mean_meanyears_to_MRCA=NA,upper_95_CI=NA,lower_95_CI=NA)
for (patient in PD_ID[1:length(PD_ID)]){
  
# read in the observed data.
obs_patient <- read.csv(paste0("/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/small_bowel/data/stem_cell/observed_VAF/4mtr/",patient,"_VAF.txt"), sep=" ", header=T, stringsAsFactors = F)
n_crypts <- length(unique(obs_patient$locus))
crypt_lst = unique(obs_patient$locus)
obsstats=c()
for (i in 1:n_crypts){
  df_sample = obs_patient[obs_patient$locus==crypt_lst[i],]
  obsstats=c(obsstats,as.vector(table(cut(df_sample$VAF, breaks=seq(from=0, to=0.7, by=0.05)))))
}


# read in the output of the simulations
sumstats <- read.csv(paste0("/Users/yw2/Documents/FARM/yw2/small_bowel/stem_cell/2022.08_abc/",patient,"/summary_statistics_100-150.txt"), sep=" ", header=F, stringsAsFactors = F)
dim(sumstats)

# I have a  "years_to_MRCA_parameter" for every crypt, and the stem cell number, gendays, and mutation rate parameters for every crypt
breaks=4+c(1:n_crypts)*15
newnams <- c('stem','gendays','mutsperdoubling')
for (i in 1:length(breaks)) {
  newnams <- c(newnams,paste0(crypt_lst[i],"_", c("years_to_MRCA","0_0.05","0.05_0.1","0.1_0.15","0.15_0.2","0.2_0.25","0.25_0.3","0.3_0.35","0.35_0.4","0.4_0.45","0.45_0.5","0.5_0.55","0.55_0.6","0.6_0.65","0.65_0.7")))
}
newnams

colnames(sumstats) <- newnams

# look at  the prior
# # plot the prior.
pdf("Prior.pdf")
cutpointstem <- seq(from=0, to=30, length.out = 7)
cutpointgendays <- seq(from=0, to=620, length.out = 21)

ex_c <- cut(sumstats$stem, breaks=cutpointstem)
ey_c <- cut(sumstats$gendays, breaks=cutpointgendays)
ez <- table(ey_c, ex_c)

par(mar=c(5,4,5,5))
image2D(z=ez, border="black", xaxt="n", yaxt="n", xlab="Generation time (months)", ylab="Stem",
        main="Joint prior distribution")
axis(side=1, at=seq(-0.02,1.02, length.out = length(cutpointgendays)), labels = c(0:20), las=2)
axis(side=2, at=seq(-0.1,1.1, length.out = length(cutpointstem)),labels = cutpointstem, las=2)

par(mfrow=c(2,1))
hist(sumstats$stem, col="grey", 50, main="Prior for stem cells")
hist(sumstats$gendays, col="grey", 50, main="Prior for generation time")
dev.off()

# flatten the prior
i <- 1
j <- 1

nsumstats <- data.frame()
for (i in 1:length(cutpointstem)) {
  print(i)
  for (j in 1:length(cutpointgendays)) {
    tbox <- sumstats[sumstats$stem>cutpointstem[i] & sumstats$stem<=cutpointstem[i+1] 
             & sumstats$gendays>cutpointgendays[j] & sumstats$gendays<=cutpointgendays[j+1],]
    if (nrow(tbox)>1000) {
      tbox <- tbox[sample(1:nrow(tbox), 1000, replace=F),]
    }
    nsumstats <- rbind(nsumstats, tbox)
    j <- j+1
  }
  i <- i+1
}
dim(nsumstats)
dim(sumstats)

x_c <- cut(nsumstats$stem, breaks=cutpointstem)
y_c <- cut(nsumstats$gendays, breaks=cutpointgendays)
z <- table(y_c, x_c)

pdf("Prior_after_trimming.pdf")
par(mar=c(5,4,5,5))
image2D(z=z, border="black", xaxt="n", yaxt="n", xlab="Generation time (months)", ylab="Stem",
        main="Joint prior distribution")
axis(side=1, at=seq(-0.02,1.02, length.out = length(cutpointgendays)), labels = c(0:20), las=2)
axis(side=2, at=seq(-0.1,1.1, length.out = length(cutpointstem)),labels = cutpointstem, las=2)
dev.off()

sumstats <- nsumstats

# compare the observed to the simulated.
reord1 <- c(colnames(sumstats)[1:3], grep("years_to_MRCA", colnames(sumstats), value=T))
reord2 <- colnames(sumstats)[!colnames(sumstats) %in% reord1]
reord <- c(reord1, reord2)
ss <- sumstats[,reord]
names(obsstats) <- reord2

# calculate the distances between the two.
for (name in names(obsstats)) {
  ss[,name] <- abs(ss[,name] - as.numeric(obsstats[name]))
}
saved <- ss

# normalise each bin.
scale.fn <- function(x) {x / sqrt(sum(x^2))}
for (name in names(obsstats)) {
  ss[,name] <- scale.fn(ss[,name])
}
head(ss)
ss[is.na(ss)] <- 0

pdf("distances_for_each_vaf_bin_after_normalising.pdf")
par(mfrow=c(5,5))
par(mar=c(1,1,1,1))
for (name in names(obsstats)) {
  hist(ss[,name], col="grey", 50, main=name, cex.main=0.5)
}
dev.off()

ss$dist4 <- rowSums(ss[,grep("PD", names(obsstats), value=T)])

hist(ss$dist4, col="grey", 50)


# reorder them based on the distance.
head(ss)

ss$meanyears_to_MRCA <- apply(ss[,grep("years_to_MRCA", colnames(ss))], 1, function(row)
  mean(row, na.rm=T))

ss4 <- ss[with(ss, order(dist4, decreasing=F)),]

# pick the best 1%
astats4 <- ss4[1:round(nrow(ss)*0.01),]

dim(astats4)

# colour them by the best.
# I want to colour the bottom 99% one colour, the top 1% another colour, 
# and the top 0.1% another colour
display.brewer.pal(10, "Spectral")
mypal <- RColorBrewer::brewer.pal(10, 'Spectral')[c(1,2,3,9)]


ss4 <- ss[with(ss, order(dist4), decreasing=T),]
ss4$discol <- mypal[4]
ss4$discol[1:round(nrow(ss)*0.01)] <- mypal[2]
ss4$discol[1:round(nrow(ss)*0.001)] <- mypal[1]


pdf(paste0("/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/small_bowel/data/stem_cell/output/",patient,".pdf"),height = 8,width = 8)
par(mfrow=c(3,1))

plot(stem~gendays, data=ss4[nrow(ss4):1,], col=ss4$discol[nrow(ss4):1], pch=16, 
     main="Simulations coloured by \n similarity to observed with 4 reads", 
     xlab="Stem cell replacement rate (days)", ylab="Stem cell effective population size")
legend("topright", legend=c("best 0.1%", "0.1-1%", "1-100%"), col=mypal[c(1,2,4)], pch=16, cex=1,
       bg="white")


hist(ss$meanyears_to_MRCA, col="grey", 50,  xlab="Years to MRCA",
     main=paste0("Mean time to MRCA, prior distribution \n CI95: ", paste(round(sort(ss$meanyears_to_MRCA)[c(round(length(ss$meanyears_to_MRCA)*0.025), round(length(ss$meanyears_to_MRCA)*0.975))], digits=1), collapse="-"), " years"))

mean(ss$meanyears_to_MRCA)
median(ss$meanyears_to_MRCA)

hist(astats4$meanyears_to_MRCA, col="grey", 50, xlab="Years to MRCA",
     main=paste0("Mean time to MRCA, best 1% estimation \n CI95: ", paste(round(sort(astats4$meanyears_to_MRCA)[c(round(length(astats4$meanyears_to_MRCA)*0.025), round(length(astats4$meanyears_to_MRCA)*0.975))], digits=1), collapse="-"), " years"))

mean(astats4$meanyears_to_MRCA) 
median(astats4$meanyears_to_MRCA) 
dev.off()



meanyears_to_MRCA$median_meanyears_to_MRCA[meanyears_to_MRCA$patient==patient] = median(astats4$meanyears_to_MRCA)
meanyears_to_MRCA$mean_meanyears_to_MRCA[meanyears_to_MRCA$patient==patient] = mean(astats4$meanyears_to_MRCA)
meanyears_to_MRCA$upper_95_CI[meanyears_to_MRCA$patient==patient] = sort(astats4$meanyears_to_MRCA)[length(astats4$meanyears_to_MRCA)*0.975]
meanyears_to_MRCA$lower_95_CI[meanyears_to_MRCA$patient==patient] = sort(astats4$meanyears_to_MRCA)[length(astats4$meanyears_to_MRCA)*0.025]
}

write.csv(meanyears_to_MRCA,'/Users/yw2/OneDrive\ -\ University\ of\ Cambridge/PhD/small_bowel/data/stem_cell/output/meanyears_to_MRCA_full.csv',quote = F,row.names = F)

