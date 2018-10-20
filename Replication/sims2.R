library(cpt)
library(parallel)


args=(commandArgs(TRUE))
perm.N = eval( parse(text=args[1]) )
simN   = eval( parse(text=args[2]) )


n = 100
treat = as.factor(c(rep(0,n/2), rep(1,n/2)))
length.beta.vec = 8
beta.vec = seq(0,0.7,length=length.beta.vec)


    sim4d = rep(NA, length.beta.vec*7*simN)
    dim(sim4d) = c(length.beta.vec, 7, simN)
    dimnames(sim4d)[[1]] = paste("beta", beta.vec, sep="")
    dimnames(sim4d)[[2]] = c("glmnet",
                             "forest",
                             "ensemble",
                             "combined",
                             "energy",
                             "crossmatch",
                             "hotelling")
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    a=clusterCall(cl, function() library(cpt))  
    a=clusterCall(cl, function() library(cpt.paper))      
    a=clusterCall(cl, function() library(energy))  
    a=clusterCall(cl, function() library(Hotelling))      
    a=clusterCall(cl, function() library(MASS)) 
    clusterSetRNGStream(cl, iseed=123456)
    for (beta.i in 1:length.beta.vec)
    {
        beta = beta.vec[beta.i]*10/7
        sim4d[beta.i, , ] = parSapply(cl, 1:simN,
                          function(sim.i, n, T, perm.N, beta)
                          {
                             Z = matrix(rnorm(n*1000),n,1000)
                             Z[1:(n/2),1:2] = Z[1:(n/2),1:2] + beta
                             #################################################################
                             class.methods = c("glmnet", "forest") 
                             metric = "logscore"
                             ensemble.metric = "mean.log" 
                             comb.methods = c(class.methods, "ensemble")
                             #################################################################
                             pvals = rep(NA, length(class.methods) + 5)
                             names(pvals) = c(class.methods, "ensemble", "combined", "energy", 
                                                             "crossmatch", "hotelling")
                             rval = cpt(Z, T, class.methods = class.methods, metric = metric, 
                                        ensemble.metric = ensemble.metric, perm.N = perm.N, comb.methods = comb.methods)
                             pvals[class.methods] = rval$pvals[class.methods]
                             pvals["ensemble"] = rval$pvals["ensemble"]
                             pvals["combined"] = rval$pval
                             groupindex = c(which(as.numeric(T) == 1), which(as.numeric(T) == 2))
                             sizes = c(sum(as.numeric(T) == 1), sum(as.numeric(T) == 2))
                             pvals["energy"] = eqdist.etest(Z[groupindex, ], sizes = sizes, R = 999)$p.value
                             #pvals["crossmatch"] = crosstest(Z, T)$approxpval
                             return(pvals)
                          },
                n=n, T=treat, perm.N=perm.N, beta=beta)
    }       
    stopCluster(cl)
    save(sim4d, file="sim4d.rda")    



    sim3d = rep(NA, length.beta.vec*7*simN)
    dim(sim3d) = c(length.beta.vec, 7, simN)
    dimnames(sim3d)[[1]] = paste("beta", beta.vec, sep="")
    dimnames(sim3d)[[2]] = c("glmnet",
                             "forest",
                             "ensemble",
                             "combined",
                             "energy",
                             "crossmatch",
                             "hotelling")
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    a=clusterCall(cl, function() library(cpt))  
    a=clusterCall(cl, function() library(cpt.paper))      
    a=clusterCall(cl, function() library(energy))  
    a=clusterCall(cl, function() library(Hotelling))      
    a=clusterCall(cl, function() library(MASS)) 
    clusterSetRNGStream(cl, iseed=123456)
    for (beta.i in 1:length.beta.vec)
    {
        beta = beta.vec[beta.i]
        sim3d[beta.i, , ] = parSapply(cl, 1:simN,
                          function(sim.i, n, T, perm.N, beta)
                          {
                             Z = matrix(rnorm(n*50),n,50)
                             Z[1:(n/2),1:2] = Z[1:(n/2),1:2] + beta
                             #################################################################
                             class.methods = c("glmnet", "forest") 
                             metric = "logscore"
                             ensemble.metric = "mean.log" 
                             comb.methods = c(class.methods, "ensemble")
                             #################################################################
                             pvals = rep(NA, length(class.methods) + 5)
                             names(pvals) = c(class.methods, "ensemble", "combined", "energy", 
                                                             "crossmatch", "hotelling")
                             rval = cpt(Z, T, class.methods = class.methods, metric = metric, 
                                        ensemble.metric = ensemble.metric, perm.N = perm.N, comb.methods = comb.methods)
                             pvals[class.methods] = rval$pvals[class.methods]
                             pvals["ensemble"] = rval$pvals["ensemble"]
                             pvals["combined"] = rval$pval
                             groupindex = c(which(as.numeric(T) == 1), which(as.numeric(T) == 2))
                             sizes = c(sum(as.numeric(T) == 1), sum(as.numeric(T) == 2))
                             pvals["energy"] = eqdist.etest(Z[groupindex, ], sizes = sizes, R = 999)$p.value
                             pvals["crossmatch"] = crosstest(Z, T)$approxpval
                             pvals["hotelling"] = hotelling.test(Z[as.numeric(T) == 1, ], Z[as.numeric(T) == 2, ])$pval
                             return(pvals)
                          },
                n=n, T=treat, perm.N=perm.N, beta=beta)
    }       
    stopCluster(cl)
    save(sim3d, file="sim3d.rda")

    
    

    
    
 

    
