library(cpt)
library(parallel)


args=(commandArgs(TRUE))
perm.N = eval( parse(text=args[1]) )
simN   = eval( parse(text=args[2]) )



n = 100
treat = as.factor(c(rep(0,n/2), rep(1,n/2)))
length.beta.vec = 8
beta.vec = seq(0,0.7,length=length.beta.vec)
length.rho.vec = 9
rho.vec = seq(0,0.8,length=length.rho.vec)


    sim1 = rep(NA, length.beta.vec*8*simN)
    dim(sim1) = c(length.beta.vec, 8, simN)
    dimnames(sim1)[[1]] = paste("beta", beta.vec, sep="")
    dimnames(sim1)[[2]] = c("logistic",
                            "logistic2",
                            "forest",
                            "ensemble",
                            "combined",
                            "energy",
                            "crossmatch",
                            "hotelling")
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    a=clusterCall(cl, function() library(cpt.paper))  
    a=clusterCall(cl, function() library(MASS)) 
    clusterSetRNGStream(cl, iseed=123456)
    for (beta.i in 1:length.beta.vec)
    {
        beta = beta.vec[beta.i]
        sim1[beta.i, , ] = parSapply(cl, 1:simN,
                          function(sim.i, n, treat, perm.N, beta)
                          {
                             Z = rbind(mvrnorm(n/2,mu=c(0,    0,    0),Sigma=matrix(c(1,0,0,0,1,0,0,0,1),                        ncol=3,byrow=T)),
                                       mvrnorm(n/2,mu=c(beta, beta, 0),Sigma=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=T)))
                             return(cpt.sim.analysis(Z, treat, perm.N=perm.N, metric="rate", ensemble.metric="vote")) 
                          },
                n=n, treat=treat, perm.N=perm.N, beta=beta)
    }       
    stopCluster(cl)
    save(sim1, file="sim1.rda")

    
    sim1d = rep(NA, length.beta.vec*8*simN)
    dim(sim1d) = c(length.beta.vec, 8, simN)
    dimnames(sim1d)[[1]] = paste("beta", beta.vec, sep="")
    dimnames(sim1d)[[2]] = c("logistic",
                            "logistic2",
                            "forest",
                            "ensemble",
                            "combined",
                            "energy",
                            "crossmatch",
                            "hotelling")
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    a=clusterCall(cl, function() library(cpt.paper))  
    a=clusterCall(cl, function() library(MASS)) 
    clusterSetRNGStream(cl, iseed=123456)
    for (beta.i in 1:length.beta.vec)
    {
        beta = beta.vec[beta.i]
        sim1d[beta.i, , ] = parSapply(cl, 1:simN,
                          function(sim.i, n, treat, perm.N, beta)
                          {
                             Z = rbind(mvrnorm(n/2,mu=c(0,    0,    0),Sigma=matrix(c(1,0,0,0,1,0,0,0,1),                        ncol=3,byrow=T)),
                                       mvrnorm(n/2,mu=c(beta, beta, 0),Sigma=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=T)))
                             return(cpt.sim.analysis(Z, treat, perm.N=perm.N, metric="logscore", ensemble.metric="mean.log")) 
                          },
                n=n, treat=treat, perm.N=perm.N, beta=beta)
    }       
    stopCluster(cl)
    save(sim1d, file="sim1d.rda")


    
    
    
    
    
    sim2 = rep(NA, length.rho.vec*8*simN)
    dim(sim2) = c(length.rho.vec, 8, simN)
    dimnames(sim2)[[1]] = paste("rho", rho.vec, sep="")
    dimnames(sim2)[[2]] = c("logistic",
                            "logistic2",
                            "forest",
                            "ensemble",
                            "combined",
                            "energy",
                            "crossmatch",
                            "hotelling")
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    a=clusterCall(cl, function() library(cpt.paper))  
    a=clusterCall(cl, function() library(MASS)) 
    clusterSetRNGStream(cl, iseed=123456)
    for (rho.i in 1:length.rho.vec)
    {
        r = rho.vec[rho.i]
        sim2[rho.i, , ] = parSapply(cl, 1:simN,
                          function(sim.i, n, treat, perm.N, r)
                          {
                             Z = rbind(mvrnorm(n/2,mu=c(0, 0, 0),Sigma=matrix(c(1,0,0,0,1,0,0,0,1),                        ncol=3,byrow=T)),
                                       mvrnorm(n/2,mu=c(0, 0, 0),Sigma=matrix(c(1,r,r,r,1,r,r,r,1),ncol=3,byrow=T)))
                             return(cpt.sim.analysis(Z, treat, perm.N=perm.N, metric="rate", ensemble.metric="vote")) 
                          },
                    n=n, treat=treat, perm.N=perm.N, r=r)
    }       
    stopCluster(cl)
    save(sim2, file="sim2.rda")

    
    
    
    sim2d = rep(NA, length.rho.vec*8*simN)
    dim(sim2d) = c(length.rho.vec, 8, simN)
    dimnames(sim2d)[[1]] = paste("rho", rho.vec, sep="")
    dimnames(sim2d)[[2]] = c("logistic",
                            "logistic2",
                            "forest",
                            "ensemble",
                            "combined",
                            "energy",
                            "crossmatch",
                            "hotelling")

    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    a=clusterCall(cl, function() library(cpt.paper))  
    a=clusterCall(cl, function() library(MASS)) 
    clusterSetRNGStream(cl, iseed=123456)
    for (rho.i in 1:length.rho.vec)
    {
        r = rho.vec[rho.i]
        sim2d[rho.i, , ] = parSapply(cl, 1:simN,
                          function(sim.i, n, treat, perm.N, r)
                          {
                             Z = rbind(mvrnorm(n/2,mu=c(0, 0, 0),Sigma=matrix(c(1,0,0,0,1,0,0,0,1),                        ncol=3,byrow=T)),
                                       mvrnorm(n/2,mu=c(0, 0, 0),Sigma=matrix(c(1,r,r,r,1,r,r,r,1),ncol=3,byrow=T)))
                             return(cpt.sim.analysis(Z, treat, perm.N=perm.N, metric="logscore", ensemble.metric="mean.log")) 
                          },
                    n=n, treat=treat, perm.N=perm.N, r=r)
    }       
    stopCluster(cl)
    save(sim2d, file="sim2d.rda")

    
