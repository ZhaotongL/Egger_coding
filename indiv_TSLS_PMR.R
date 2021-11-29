library(TwoSampleMR)
library(nlme)
library(mixIE)
library(MRcML)
library(PMR)

two_runif <- function(m,min1,min2,max1,max2){
  y <- runif(m, 0, max1-min1+max2-min2)
  y[which(y<(max1-min1))] = min1+y[which(y<(max1-min1))]
  y[which(y>=(max1-min1))] = min2+y[which(y>=(max1-min1))]-(max1-min1)
  return(y)
}

generate_gwas_DP_1 <- function(seed,n,m,theta,K){
  set.seed(seed)
  gx = lapply(1:m,function(x){rbinom(n,2,0.3)})
  gx = do.call(cbind,gx)
  bxg = two_runif(m,min1=-0.2,max1=-0.1,min2=0.1,max2=0.2)
  byg = rep(0,m)
  bug = rep(0,m)
  if(K>=1){
    byg[1:K] = rnorm(K,mean=0.1,sd=0.1)
  }
  eu = rnorm(n,0,1)
  ex = rnorm(n,0,1)
  ey = rnorm(n,0,1)
  u = gx %*% bug + eu
  x = gx %*% bxg + u + ex
  y = theta*x + gx %*% byg + u + ey
  
  ind1 = 1:(n/2) #first half sample
  ind2 = (n/2+1):n #second half sample
  gwasx = lapply(1:m,function(i){lm1=lm(x[ind1]~gx[ind1,i]);
  summary(lm1)$coefficients[2,1:2]})
  betax = unname(unlist(lapply(gwasx,function(x){x[1]})))
  to_flip = which(betax<0)
  
  x1 = x[ind1]
  y2 = y[ind2]
  gx1 = gx[ind1,]
  gx2 = gx[ind2,]
  beta_hat = coef(lm(x1~gx1))
  xhat_2 = gx2 %*% beta_hat[-1] + beta_hat[1]
  
  Id <- factor(rep(1, length = length(y2)))
  X1snps_sum_vec = apply(gx2,1,sum)
  m1.block<-list(Id = pdIdent(~gx2-1))
  stage2_egger_lm = lm(y2~xhat_2+X1snps_sum_vec)
  oracle_egger_lm_b = summary(stage2_egger_lm)$coefficients[2,1]
  oracle_egger_lm_pval = summary(stage2_egger_lm)$coefficients[2,4]
  
  stage2_egger_lme = tryCatch(lme(y2 ~ xhat_2 + X1snps_sum_vec, random = m1.block,method='ML',control = lmeControl(opt = "optim")),
                              error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  oracle_egger_lme_b = summary(stage2_egger_lme)$tTable[2,1]
  oracle_egger_lme_pval = summary(stage2_egger_lme)$tTable[2,5]
  
  stage2_ivw_lme = tryCatch(lme(y2 ~ xhat_2 , random = m1.block,method='ML',control = lmeControl(opt = "optim")),
                            error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  ivw_lme_b = summary(stage2_ivw_lme)$tTable[2,1]
  ivw_lme_pval = summary(stage2_ivw_lme)$tTable[2,5]
  
  oracle_PMR_res = PMR_individual(yin=scale(x1),
                                  zin=scale(y2),
                                  x1in=scale(gx1),
                                  x2in=scale(gx2))
  oracle_PMR_b = oracle_PMR_res$causal_effect
  oracle_PMR_pval = oracle_PMR_res$causal_pvalue
  
  gx1_flip = gx[ind1,]
  gx1_flip[,to_flip] = 2 - gx1_flip[,to_flip]
  gx2_flip = gx[ind2,]
  gx2_flip[,to_flip] = 2 - gx2_flip[,to_flip]
  beta_hat_flip = coef(lm(x1~gx1_flip))
  xhat_2_flip = gx2_flip %*% beta_hat_flip[-1] + beta_hat_flip[1]
  
  X1snps_sum_vec_flip = apply(gx2_flip,1,sum)
  m1.block_flip<-list(Id = pdIdent(~gx2_flip-1))
  stage2_egger_default_lm = lm(y2~xhat_2_flip+X1snps_sum_vec_flip)
  default_egger_lm_b = summary(stage2_egger_default_lm)$coefficients[2,1]
  default_egger_lm_pval = summary(stage2_egger_default_lm)$coefficients[2,4]
  
  stage2_egger_default_lme = tryCatch(lme(y2 ~ xhat_2_flip + X1snps_sum_vec_flip, random = m1.block_flip,
                                          method='ML',control = lmeControl(opt = "optim")),
                                      error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  default_egger_lme_b = summary(stage2_egger_default_lme)$tTable[2,1]
  default_egger_lme_pval = summary(stage2_egger_default_lme)$tTable[2,5]
  
  default_PMR_res = PMR_individual(yin=scale(x1),
                                   zin=scale(y2),
                                   x1in=scale(gx1_flip),
                                   x2in=scale(gx2_flip))
  default_PMR_b = default_PMR_res$causal_effect
  default_PMR_pval = default_PMR_res$causal_pvalue
  
  out = list()
  out$oracle.lm_b = oracle_egger_lm_b
  out$oracle.lm_pval = oracle_egger_lm_pval
  out$oracle.lme_b = oracle_egger_lme_b
  out$oracle.lme_pval = oracle_egger_lme_pval
  out$default.lm_b = default_egger_lm_b
  out$default.lm_pval = default_egger_lm_pval
  out$default.lme_b = default_egger_lme_b
  out$default.lme_pval = default_egger_lme_pval
  out$ivw.lme_b = ivw_lme_b
  out$ivw.lme_pval = ivw_lme_pval
  out$oracle.PMR_b = oracle_PMR_b
  out$oracle.PMR_pval = oracle_PMR_pval
  out$default.PMR_b = default_PMR_b
  out$default.PMR_pval = default_PMR_pval
  
  return(out)
}

generate_gwas_DP_2 <- function(seed,n,m,theta,K){
  set.seed(seed)
  gx = lapply(1:m,function(x){rbinom(n,2,0.3)})
  gx = do.call(cbind,gx)
  bxg = two_runif(m,min1=-0.1,max1=-0.03,min2=0.1,max2=0.2)
  byg = rep(0,m)
  bug = rep(0,m)
  if(K>=1){
    byg[1:K] = rnorm(K,mean=0.1,sd=0.1)
  }
  eu = rnorm(n,0,1)
  ex = rnorm(n,0,1)
  ey = rnorm(n,0,1)
  u = gx %*% bug + eu
  x = gx %*% bxg + u + ex
  y = theta*x + gx %*% byg + u + ey
  
  ind1 = 1:(n/2) #first half sample
  ind2 = (n/2+1):n #second half sample
  gwasx = lapply(1:m,function(i){lm1=lm(x[ind1]~gx[ind1,i]);
  summary(lm1)$coefficients[2,1:2]})
  betax = unname(unlist(lapply(gwasx,function(x){x[1]})))
  to_flip = which(betax<0)
  
  x1 = x[ind1]
  y2 = y[ind2]
  gx1 = gx[ind1,]
  gx2 = gx[ind2,]
  beta_hat = coef(lm(x1~gx1))
  xhat_2 = gx2 %*% beta_hat[-1] + beta_hat[1]
  
  Id <- factor(rep(1, length = length(y2)))
  X1snps_sum_vec = apply(gx2,1,sum)
  m1.block<-list(Id = pdIdent(~gx2-1))
  stage2_egger_lm = lm(y2~xhat_2+X1snps_sum_vec)
  oracle_egger_lm_b = summary(stage2_egger_lm)$coefficients[2,1]
  oracle_egger_lm_pval = summary(stage2_egger_lm)$coefficients[2,4]
  
  stage2_egger_lme = tryCatch(lme(y2 ~ xhat_2 + X1snps_sum_vec, random = m1.block,method='ML',control = lmeControl(opt = "optim")),
                              error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  oracle_egger_lme_b = summary(stage2_egger_lme)$tTable[2,1]
  oracle_egger_lme_pval = summary(stage2_egger_lme)$tTable[2,5]
  
  stage2_ivw_lme = tryCatch(lme(y2 ~ xhat_2 , random = m1.block,method='ML',control = lmeControl(opt = "optim")),
                            error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  ivw_lme_b = summary(stage2_ivw_lme)$tTable[2,1]
  ivw_lme_pval = summary(stage2_ivw_lme)$tTable[2,5]
  
  oracle_PMR_res = PMR_individual(yin=scale(x1),
                                  zin=scale(y2),
                                  x1in=scale(gx1),
                                  x2in=scale(gx2))
  oracle_PMR_b = oracle_PMR_res$causal_effect
  oracle_PMR_pval = oracle_PMR_res$causal_pvalue
  
  gx1_flip = gx[ind1,]
  gx1_flip[,to_flip] = 2 - gx1_flip[,to_flip]
  gx2_flip = gx[ind2,]
  gx2_flip[,to_flip] = 2 - gx2_flip[,to_flip]
  beta_hat_flip = coef(lm(x1~gx1_flip))
  xhat_2_flip = gx2_flip %*% beta_hat_flip[-1] + beta_hat_flip[1]
  
  X1snps_sum_vec_flip = apply(gx2_flip,1,sum)
  m1.block_flip<-list(Id = pdIdent(~gx2_flip-1))
  stage2_egger_default_lm = lm(y2~xhat_2_flip+X1snps_sum_vec_flip)
  default_egger_lm_b = summary(stage2_egger_default_lm)$coefficients[2,1]
  default_egger_lm_pval = summary(stage2_egger_default_lm)$coefficients[2,4]
  
  stage2_egger_default_lme = tryCatch(lme(y2 ~ xhat_2_flip + X1snps_sum_vec_flip, random = m1.block_flip,
                                          method='ML',control = lmeControl(opt = "optim")),
                                      error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  default_egger_lme_b = summary(stage2_egger_default_lme)$tTable[2,1]
  default_egger_lme_pval = summary(stage2_egger_default_lme)$tTable[2,5]
  
  default_PMR_res = PMR_individual(yin=scale(x1),
                                   zin=scale(y2),
                                   x1in=scale(gx1_flip),
                                   x2in=scale(gx2_flip))
  default_PMR_b = default_PMR_res$causal_effect
  default_PMR_pval = default_PMR_res$causal_pvalue
  
  out = list()
  out$oracle.lm_b = oracle_egger_lm_b
  out$oracle.lm_pval = oracle_egger_lm_pval
  out$oracle.lme_b = oracle_egger_lme_b
  out$oracle.lme_pval = oracle_egger_lme_pval
  out$default.lm_b = default_egger_lm_b
  out$default.lm_pval = default_egger_lm_pval
  out$default.lme_b = default_egger_lme_b
  out$default.lme_pval = default_egger_lme_pval
  out$ivw.lme_b = ivw_lme_b
  out$ivw.lme_pval = ivw_lme_pval
  out$oracle.PMR_b = oracle_PMR_b
  out$oracle.PMR_pval = oracle_PMR_pval
  out$default.PMR_b = default_PMR_b
  out$default.PMR_pval = default_PMR_pval
  
  return(out)
}

run_simulation_egger <- function(seed,n,m,theta,K,type){
  generate_gwas = switch(type,'DP1'=generate_gwas_DP_1,'DP2'=generate_gwas_DP_2,
                         'BP1'=generate_gwas_BP_1,'BP2'=generate_gwas_BP_2)
  out = generate_gwas(seed,n,m,theta,K)
  return(out)
}