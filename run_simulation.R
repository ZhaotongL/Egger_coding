source('./main_simulation.R')
library(parallel)
library(tidyr)
ncore <- detectCores() -2


byx = c(0,0.2)
n = c(100000)
p_invalid = c(0,0.3,0.7,1)
m = c(30,100,500)
paras = crossing(n,m,byx,p_invalid)

input_type = 'DP2' # Simulation(b):directional pleiotropy

for(i in 1:nrow(paras)){
  K = ceiling(paras$m[i] * paras$p_invalid[i])
  out1 = mclapply(1:500,function(x){run_simulation_egger(seed=x, n=paras$n[i], m=paras$m[i], theta=paras$byx[i],  K=K ,type=input_type)},mc.cores=ncore)
  
  error_ind = which(unlist(lapply(out1,function(x){class(x)=='try-error'})))
  if(length(error_ind)>0){
    print(error_ind)
    out1 = out1[-error_ind]}
  var_name = names(out1[[1]])
  for(vn in var_name){
    t = paste0(vn,'<-unlist(lapply(out1,function(x){x$',vn,'}))')
    eval(parse(text=t))
    t1 = paste0("write.table(t(c(paras[i,],",vn,")),paste0(input_type,'",vn,".txt'),quote=F,row.names=F,append=T,col.names=F)")
    eval(parse(text=t1))
  }
  
  print(i)
}