mr_egger_ll <-function (b_exp, b_out, se_exp, se_out,flip=0)
{
  stopifnot(length(b_exp) == length(b_out))
  stopifnot(length(se_exp) == length(se_out))
  stopifnot(length(b_exp) == length(se_out))
  nulllist <- list(b = NA, se = NA, pval = NA, nsnp = NA, b_i = NA,
                   se_i = NA, pval_i = NA,sig=NA)
  if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) &
          !is.na(se_out)) < 3) {
    return(nulllist)
  }
  if(flip==1){
    sign0 <- function(x) {
      x[x == 0] <- 1
      return(sign(x))
    }
    to_flip <- sign0(b_exp) == -1
    b_out = b_out * sign0(b_exp)
    b_exp = abs(b_exp)
  }
  mod <- stats::lm(b_out ~ b_exp, weights = 1/se_out^2)
  smod <- summary(mod)
  if (nrow(stats::coefficients(smod)) > 1) {
    b <- stats::coefficients(smod)[2, 1]
    se <- stats::coefficients(smod)[2, 2]/min(1, smod$sigma)
    pval <- 2 * stats::pt(abs(b/se), length(b_exp) - 2, lower.tail = FALSE)
    b_i <- stats::coefficients(smod)[1, 1]
    se_i <- stats::coefficients(smod)[1, 2]/min(1, smod$sigma)
    pval_i <- 2 * stats::pt(abs(b_i/se_i), length(b_exp) - 2, lower.tail = FALSE)
    lcom <- sum(log(L1_i(b_exp,b_out,se_exp,se_out,b,b_i,max(1,(sqrt((length(b_exp)-2)/length(b_exp))*smod$sigma)^2))))
  }
  else {
    warning("Collinearities in MR Egger, try LD pruning the exposure variables.")
    return(nulllist)
  }
  return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), sig=smod$sigma,
              b_i = b_i, se_i = se_i, pval_i = pval_i,mod=mod,smod=smod,lcom=lcom))
}

egger_flip <- function(b_exp, b_out, se_exp, se_out,flip_n){
  m = length(b_exp)
  theta_B = pval_B = se_B = m_B = lcom_B = r_B = rpval_B = rep(NA,flip_n)
  flip_B = matrix(FALSE,nrow=flip_n,ncol=m)
  #set.seed(seed)
  for(i in 1:flip_n){
    stf = sample((0):m,size=1)
    m_B[i] = stf
    to_flip = sample(1:m, size = stf, replace = FALSE)
    b_exp_B = b_exp
    b_out_B = b_out
    b_exp_B[to_flip] = -b_exp_B[to_flip]
    b_out_B[to_flip] = -b_out_B[to_flip]
    egger_B = mr_egger_ll(b_exp_B,b_out_B,se_exp,se_out,flip=0)
    theta_B[i] = egger_B$b
    pval_B[i] = egger_B$pval
    se_B[i] = egger_B$se
    r_B[i] = egger_B$b_i
    rpval_B[i] = egger_B$pval_i
    lcom_B[i] = egger_B$lcom
    
    flip_B[i,to_flip] = TRUE
  }
  out = list()
  dup = which(duplicated(theta_B))
  if(length(dup)>0){
    theta_B = theta_B[-dup]
    pval_B = pval_B[-dup]
    se_B = se_B[-dup]
    r_B = r_B[-dup]
    rpval_B = rpval_B[-dup]
    m_B = m_B[-dup]
    lcom_B = lcom_B[-dup]
    flip_B = flip_B[-dup,,drop=F]
  }
  lcom_order = order(lcom_B,decreasing = T)
  
  out$theta = theta_B[lcom_order]
  out$pval = pval_B[lcom_order]
  out$se = se_B[lcom_order]
  out$m = m_B[lcom_order]
  out$r = r_B[lcom_order]
  out$rpval = rpval_B[lcom_order]
  out$lcom = lcom_B[lcom_order]
  return(out)
}
