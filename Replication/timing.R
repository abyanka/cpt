library(cpt)
library(cpt.paper)
library(Hotelling)
library(energy)
library(crossmatch)
set.seed(123456)

perm.N=500

time_methods = function(n=100,q=3,perm.N=500)
{
  print(n)
  print(q)
  methods = c("hotelling", "crossmatch", "energy", "logistic", "forest", "glmnet")
  times = rep(NA, length(methods))
  names(times)=methods
  T = as.factor(c(rep(0,n/2), rep(1,n/2)))
  Z = matrix(rnorm(n*q),n,q)
  groupindex = c(which(as.numeric(T) == 1), which(as.numeric(T) == 2))
  sizes = c(sum(as.numeric(T) == 1), sum(as.numeric(T) == 2))    
  #################################################  
  if (q < n)
  {
    print("hotelling")
    time0 = proc.time()
    a=hotelling.test(Z[as.numeric(T) == 1, ], Z[as.numeric(T) == 2, ])$pval
    times[1] = (proc.time()-time0)[3]
  }
  #################################################  
  print("crossmatch")
  time0 = proc.time()
  a=crosstest(Z, T)$approxpval
  times[2] = (proc.time()-time0)[3]
  #################################################  
  print("energy")
  time0 = proc.time()
  a=eqdist.etest(Z[groupindex, ], sizes = sizes, R = 999)$p.value
  times[3] = (proc.time()-time0)[3]
  #################################################
  if (q < n)
  {
    print("logistic")  
    time0 = proc.time()
    a=cpt(Z, T, class.methods = c("logistic"), perm.N = perm.N)$pval
    times[4] = (proc.time()-time0)[3]
  }
  #################################################  
  print("forest")
  time0 = proc.time()
  a=cpt(Z, T, class.methods = c("forest"), perm.N = perm.N)$pval
  times[5] = (proc.time()-time0)[3]
  #################################################  
  print("glmnet")
  time0 = proc.time()
  a=cpt(Z, T, class.methods = c("glmnet"), perm.N = perm.N)$pval
  times[6] = (proc.time()-time0)[3]
  #################################################  
  return(times)
}

ns = c(100,500)
qs = c(5,50,1000)
runtimes = rep(NA,length(ns)*length(qs)*6)
dim(runtimes) = c(length(ns), length(qs), 6)
dimnames(runtimes) = list(paste0("n", ns), paste0("q", qs), c("hotelling", "crossmatch", "energy", "logistic", "forest", "glmnet"))
for (ni in 1:length(ns)) for (qi in 1:length(qs)) 
{
  print(ni)
  print(qi)
  runtimes[ni, qi, ] = time_methods(n=ns[ni], q=qs[qi], perm.N=perm.N)
}

save(runtimes, file="output/runtimes.rda")
print(runtimes)
