#-------------------------------------------------
# Binomial mixture model using EM algorithm
# Can incorporate truncated binomial
# Tim Coorens - April 2018
#-------------------------------------------------

#-------------------------------------------------
# Functions
#-------------------------------------------------

## Define the truncated binomial distribution
dbinomtrunc = function(x, size, prob, minx=4) {
  dbinom(x, size, prob) / pbinom(minx-0.1, size, prob, lower.tail=F)
}

estep = function(x,size,p.vector,prop.vector,ncomp, mode){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    if(mode=="Truncated") p.mat_estep[,i]=prop.vector[i]*dbinomtrunc(x,size,prob=p.vector[i])
    if(mode=="Full") p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood
  
  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }
  
  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}

## Maximisation step
mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)
  
  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)   
}

## EM algorithm
em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust,binom_mode){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities 
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust,mode=binom_mode)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust,mode=binom_mode)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]
    
    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}
    
    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]
    
  }
  if(!flag) warning("Didn’t converge\n")
  
  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}

binom_mix = function(x,size,nrange=1:5,criterion="BIC",maxit=5000,tol=1e-6, mode="Full"){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC) 
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol,binom_mode=mode)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}
#-------------------------------------------------
if(0){

# Some dummy data
n1 = 100; p1 = 0.15 # Two clones of different sizes (n) and underlying VAF (p)
n2 = 50; p2 = 0.35
depth=30

#Generate vectors of read depth (trials) and number of reads supporting variants (successes)

NR=rpois(n=n1+n2,lambda=depth)
NV=rep(0,n1+n2)
for (n in 1:n1) NV[n]=rbinom(n=1,prob=p1,size=NR[n])
for (n in (n1+1):(n1+n2)) NV[n]=rbinom(1,prob=p2,size=NR[n])
NR=NR[NV>3]
NV=NV[NV>3]

Mode="Full" #or set to "Full"
if (Mode=='Truncated'){
  NR=NR[NV>3]
  NV=NV[NV>3]
}

res = binom_mix(NV,NR,mode=Mode)

######
# Output 'res' has the following information:
# LL - log-likelihood vector of the iteration
# prop - the optimal mixing proportion of clones
# p - the probabilities (VAFs) associated with each clone
# BIC - the Bayesian information criterion associated with best fit
# AIC - the Akaike information criterion associated with best fit
# n - the optimal number of clusters/clones
# Which_cluster - same length as NV/NR input, which mutation belongs in which clone/cluster
# BIC_vec - the vector of BIC for entire range


p=hist(NV/NR,breaks=20,xlim=c(0,1),col='gray',freq=F,xlab="Variant Allele Frequency",
       main=paste0("Test (n=",length(NV),")"))
cols=c("red","blue","green","magenta","cyan")

y_coord=max(p$density)-0.5
y_intv=y_coord/5
text(y=y_coord,x=0.9,label='Data')
segments(lwd=2,lty='dashed',col='black',y0=y_coord+0.25,x0=0.85,x1=0.95)

for (i in 1:res$n){
  depth=rpois(n=5000,lambda=median(NR))
  sim_NV=unlist(lapply(depth,rbinom,n=1,prob=res$p[i]))
  sim_VAF=sim_NV/depth
  if (Mode=="Truncated") sim_VAF=sim_VAF[sim_NV>3]
  dens=density(sim_VAF)
  lines(x=dens$x,y=res$prop[i]*dens$y,lwd=2,lty='dashed',col=cols[i])
  y_coord=y_coord-y_intv/2
  text(y=y_coord,x=0.9,label=paste0("p1: ",round(res$p[i],digits=2)))
  segments(lwd=2,lty='dashed',col=cols[i],y0=y_coord+y_intv/4,x0=0.85,x1=0.95)
}

}
