library(TwoSampleMR)
library(PMR)
library(mixIE)
library(MRcML)
source('./mr_egger_func.R')

two_runif <- function(m,min1,min2,max1,max2){
  y <- runif(m, 0, max1-min1+max2-min2)
  y[which(y<(max1-min1))] = min1+y[which(y<(max1-min1))]
  y[which(y>=(max1-min1))] = min2+y[which(y>=(max1-min1))]-(max1-min1)
  return(y)
}

generate_gwas_BP_1 <- function(seed,n,m,theta,K){
  set.seed(seed)
  gx = lapply(1:m,function(x){rbinom(n,2,0.3)})
  gx = do.call(cbind,gx)
  bxg = two_runif(m,min1=-0.2,max1=-0.1,min2=0.1,max2=0.2)
  byg = rep(0,m)
  bug = rep(0,m)
  if(K>=1){
    byg[1:K] = rnorm(K,0,0.1)
  }
  eu = rnorm(n,0,1)
  ex = rnorm(n,0,1)
  ey = rnorm(n,0,1)
  u = gx %*% bug + eu
  x = gx %*% bxg + u + ex
  y = theta*x + gx %*% byg + u + ey
  heri=var(gx %*% bxg)/var(x)
  
  ind1 = 1:(n/2) #first half sample
  ind2 = (n/2+1):n #second half sample
  gwasx = lapply(1:m,function(i){lm1=lm(x[ind1]~gx[ind1,i]);
  summary(lm1)$coefficients[2,1:2]})
  betax = unname(unlist(lapply(gwasx,function(x){x[1]})))
  sdx = unname(unlist(lapply(gwasx,function(x){x[2]})))
  gwasy = lapply(1:m,function(i){lm2=lm(y[ind2]~gx[ind2,i]);
  summary(lm2)$coefficients[2,1:2]})
  betay = unname(unlist(lapply(gwasy,function(x){x[1]})))
  sdy = unname(unlist(lapply(gwasy,function(x){x[2]})))
  are = sum(betax^2)/sum((betax-mean(betax))^2)
  
  return(list(byg=byg,bxg=bxg,betax=betax,sdx=sdx,betay=betay,sdy=sdy,heri=heri,gx=gx,are=are))
}

generate_gwas_BP_2 <- function(seed,n,m,theta,K){
  set.seed(seed)
  gx = lapply(1:m,function(x){rbinom(n,2,0.3)})
  gx = do.call(cbind,gx)
  bxg = two_runif(m,min1=-0.1,max1=-0.03,min2=0.1,max2=0.2)
  byg = rep(0,m)
  bug = rep(0,m)
  if(K>=1){
    byg[1:K] = rnorm(K,0,0.1)
  }
  eu = rnorm(n,0,1)
  ex = rnorm(n,0,1)
  ey = rnorm(n,0,1)
  u = gx %*% bug + eu
  x = gx %*% bxg + u + ex
  y = theta*x + gx %*% byg + u + ey
  heri=var(gx %*% bxg)/var(x)
  
  ind1 = 1:(n/2) #first half sample
  ind2 = (n/2+1):n #second half sample
  gwasx = lapply(1:m,function(i){lm1=lm(x[ind1]~gx[ind1,i]);
  summary(lm1)$coefficients[2,1:2]})
  betax = unname(unlist(lapply(gwasx,function(x){x[1]})))
  sdx = unname(unlist(lapply(gwasx,function(x){x[2]})))
  gwasy = lapply(1:m,function(i){lm2=lm(y[ind2]~gx[ind2,i]);
  summary(lm2)$coefficients[2,1:2]})
  betay = unname(unlist(lapply(gwasy,function(x){x[1]})))
  sdy = unname(unlist(lapply(gwasy,function(x){x[2]})))
  are = sum(betax^2)/sum((betax-mean(betax))^2)
  
  return(list(byg=byg,bxg=bxg,betax=betax,sdx=sdx,betay=betay,sdy=sdy,heri=heri,gx=gx,are=are))
}

generate_gwas_DP_1 <- function(seed,n,m,theta,K){
  set.seed(seed)
  gx = lapply(1:m,function(x){rbinom(n,2,0.3)})
  gx = do.call(cbind,gx)
  bxg = two_runif(m,min1=-0.2,max1=-0.1,min2=0.1,max2=0.2)
  byg = rep(0,m)
  bug = rep(0,m)
  if(K>=1){
    byg[1:K] = rnorm(K,0.1,0.1)
  }
  eu = rnorm(n,0,1)
  ex = rnorm(n,0,1)
  ey = rnorm(n,0,1)
  u = gx %*% bug + eu
  x = gx %*% bxg + u + ex
  y = theta*x + gx %*% byg + u + ey
  heri=var(gx %*% bxg)/var(x)
  
  ind1 = 1:(n/2) #first half sample
  ind2 = (n/2+1):n #second half sample
  gwasx = lapply(1:m,function(i){lm1=lm(x[ind1]~gx[ind1,i]);
  summary(lm1)$coefficients[2,1:2]})
  betax = unname(unlist(lapply(gwasx,function(x){x[1]})))
  sdx = unname(unlist(lapply(gwasx,function(x){x[2]})))
  gwasy = lapply(1:m,function(i){lm2=lm(y[ind2]~gx[ind2,i]);
  summary(lm2)$coefficients[2,1:2]})
  betay = unname(unlist(lapply(gwasy,function(x){x[1]})))
  sdy = unname(unlist(lapply(gwasy,function(x){x[2]})))
  are = sum(betax^2)/sum((betax-mean(betax))^2)
  
  return(list(byg=byg,bxg=bxg,betax=betax,sdx=sdx,betay=betay,sdy=sdy,heri=heri,gx=gx,are=are))
}

generate_gwas_DP_2 <- function(seed,n,m,theta,K){
  set.seed(seed)
  gx = lapply(1:m,function(x){rbinom(n,2,0.3)})
  gx = do.call(cbind,gx)
  bxg = two_runif(m,min1=-0.1,max1=-0.03,min2=0.1,max2=0.2)
  byg = rep(0,m)
  bug = rep(0,m)
  if(K>=1){
    byg[1:K] = rnorm(K,0.1,0.1)
  }
  eu = rnorm(n,0,1)
  ex = rnorm(n,0,1)
  ey = rnorm(n,0,1)
  u = gx %*% bug + eu
  x = gx %*% bxg + u + ex
  y = theta*x + gx %*% byg + u + ey
  heri=var(gx %*% bxg)/var(x)
  
  ind1 = 1:(n/2) #first half sample
  ind2 = (n/2+1):n #second half sample
  gwasx = lapply(1:m,function(i){lm1=lm(x[ind1]~gx[ind1,i]);
  summary(lm1)$coefficients[2,1:2]})
  betax = unname(unlist(lapply(gwasx,function(x){x[1]})))
  sdx = unname(unlist(lapply(gwasx,function(x){x[2]})))
  gwasy = lapply(1:m,function(i){lm2=lm(y[ind2]~gx[ind2,i]);
  summary(lm2)$coefficients[2,1:2]})
  betay = unname(unlist(lapply(gwasy,function(x){x[1]})))
  sdy = unname(unlist(lapply(gwasy,function(x){x[2]})))
  are = sum(betax^2)/sum((betax-mean(betax))^2)
  
  return(list(byg=byg,bxg=bxg,betax=betax,sdx=sdx,betay=betay,sdy=sdy,heri=heri,gx=gx,are=are))
}

run_simulation_egger <- function(seed,n,m,theta,K,type){
    generate_gwas = switch(type,'DP1'=generate_gwas_DP_1,'DP2'=generate_gwas_DP_2,
                           'BP1'=generate_gwas_BP_1,'BP2'=generate_gwas_BP_2)
    s = generate_gwas(seed,n,m,theta,K)
    b_exp=s$betax; se_exp = s$sdx; b_out = s$betay;se_out = s$sdy
    
    ivw = TwoSampleMR::mr_ivw(b_exp=s$betax,se_exp = s$sdx,b_out = s$betay,se_out = s$sdy)
    egger_oracle = mr_egger_ll(b_exp=s$betax,se_exp = s$sdx,b_out = s$betay,se_out = s$sdy,flip=0)
    egger_default = TwoSampleMR::mr_egger_regression(b_exp,b_out,se_exp,se_out)
    b_exp_flip = egger_default$dat$b_exp
    b_out_flip = egger_default$dat$b_out
  
    b_exp_random = b_exp_flip
    m_random_flip = sample(1:(m-1),size=1)
    random_flip = sample(1:m,size=m_random_flip,replace=FALSE)
    b_exp_random[random_flip] = -b_exp_random[random_flip]
    b_out_random = b_out_flip
    b_out_random[random_flip] = -b_out_random[random_flip]
    egger_random = mr_egger_ll(b_exp=b_exp_random,se_exp = s$sdx,b_out =b_out_random,se_out = s$sdy,flip=0)
  
    cMLdp_result = mr_cML_DP(b_exp = s$betax,b_out=s$betay,se_exp = s$sdx,se_out = s$sdy,n=n/2,num_pert=200)
    cMLMA_b = cMLdp_result$MA_BIC_theta; cMLMA_se=cMLdp_result$MA_BIC_se; cMLMA_pval=cMLdp_result$MA_BIC_p;
    cMLMADP_b = cMLdp_result$MA_BIC_DP_theta; cMLMADP_se=cMLdp_result$MA_BIC_DP_se; cMLMADP_pval=cMLdp_result$MA_BIC_DP_p;

    CEMdp_result = mixIE_MA_DP(b_exp = s$betax,b_out=s$betay,se_exp = s$sdx,se_out = s$sdy,flip=1,n=n/2)
    CEMMA_b = CEMdp_result$mixIE_MA_theta; CEMMA_se = CEMdp_result$mixIE_MA_se; CEMMA_pval = CEMdp_result$mixIE_MA_pval; CEMMA_pi = CEMdp_result$mixIE_MA_pi;
    CEMMADP_b = CEMdp_result$mixIE_MA_DP_theta; CEMMADP_se = CEMdp_result$mixIE_MA_DP_se; CEMMADP_pval = CEMdp_result$mixIE_MA_DP_pval;
     
    out = list()
    out$ivw_b = ivw$b
    out$ivw_pval = ivw$pval
    out$oracle_b = egger_oracle$b
    out$oracle_pval = egger_oracle$pval
    out$default_b = egger_default$b
    out$default_pval = egger_default$pval
    out$random_b = egger_random$b
    out$random_pval = egger_random$pval
    out$cMLMA_b = cMLMA_b
    out$cMLMA_pval = cMLMA_pval
    out$cMLMADP_b = cMLMADP_b
    out$cMLMADP_pval = cMLMADP_pval
    out$mixIEMA_b = CEMMA_b
    out$mixIEMA_pval = CEMMA_pval
    out$mixIEMADP_b = CEMMADP_b
    out$mixIEMADP_pval = CEMMADP_pval
    return(out)
}
    
  