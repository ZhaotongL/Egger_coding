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