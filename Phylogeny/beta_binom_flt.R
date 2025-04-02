# From Tim Coorens

require(VGAM)
estimateRho_gridml = function(NV_vec,NR_vec) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}

beta.binom.filter = function(NR,NV,cutoff=0.1, binom.pval=F,pval.cutoff=0.05){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  # Optionally calculates pvalue of likelihood beta-binomial with estimated rho
  # fits better than binomial. This was supposed to protect against low-depth variants,
  # but use with caution. Returns logical vector with good variants = TRUE
  
  rho_est = pval = rep(NA,nrow(NR))
  for (k in 1:nrow(NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(NV[k,]),
                                  NR_vec=as.numeric(NR[k,]))
    if(binom.pval){
      mu = sum(as.numeric(NV[k,]))/sum(as.numeric(NR[k,]))
      LL0 = sum(dbinom(as.numeric(NV[k,]),as.numeric(NR[k,]),prob=mu))
      LL1 = sum(dbetabinom(as.numeric(NV[k,]),as.numeric(NR[k,]),prob=mu,rho=rho_est[k]))
      pval[k] = (1-pchisq(2*(LL1-LL0),df=1)) / 2
    }
    if (k%%1000==0){
      print(k)
    }
  }
  if(binom.pval){
    qval=p.adjust(pval,method="BH")
    flt_rho=qval<=pval.cutoff&rho_est>cutoff
  }else{
    flt_rho=rho_est>=cutoff
  }
  return(rho_est)
}

