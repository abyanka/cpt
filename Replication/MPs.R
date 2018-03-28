
rm(list=ls())
set.seed(1297493)

library(cpt)
library(cpt.paper)
library(crossmatch)
library(MASS)
library(energy)
library(dplyr)
library(parallel)

data(list="cpt.MPs")
Z.all = select(cpt.MPs$rawdata, -wobl_margin, -treated, -wobl_race_region,-pre_wobl)
nonmiss = complete.cases(Z.all)
Z.all = Z.all[nonmiss,]
tr = cpt.MPs$rawdata$treated[nonmiss]


##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
perm = eval( parse(text=args[1]) )
num.window  = eval( parse(text=args[2]) )
cat("\n Number of permutations is: ",perm,"\n \n")

# Parameters:
perm.N = perm


do_analysis = function (windowsize, Z, tr, wobl_margin, perm.N = 1000) 
{
    keep = abs(wobl_margin) <= windowsize
    T = as.factor(tr[keep])
    Z = as.matrix(Z[keep,])
    rval = list()
    rval[["cpt"]] = cpt(Z, T, class.methods = "forest", 
      metric = "logscore", perm.N = perm.N) 
    groupindex = c(which(as.numeric(T) == 1), which(as.numeric(T) == 2))
    sizes = c(sum(as.numeric(T) == 1), sum(as.numeric(T) == 2))
    rval[["energy"]] = try(eqdist.etest(Z[groupindex, ], sizes = sizes, R = 999))
    rval[["crossmatch"]] = try(crosstest(Z, T))
    rval2 = c(cpt.combined=rval$cpt$pval, rval$cpt$pvals, crossmatch=NA, energy=rval$energy$p.value)
    if (class(rval$crossmatch) != "try-error") rval2["crossmatch"] = rval$crossmatch$approxpval
    return(rval2)
}


### MPs
MPs.analysis = list()
MPs.analysis$windowsizes = seq(0.01,0.15, length.out=num.window)
MPs.analysis$observations = rep(0, length(MPs.analysis$windowsizes))
MPs.analysis$analyses = list()
for (j in 1:length(MPs.analysis$windowsizes))
{
  keep = abs(cpt.MPs$wobl_margin) <= MPs.analysis$windowsizes[j]
  MPs.analysis$observations[j] = sum(keep)
}

    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    a=clusterCall(cl, function() library(cpt)) 
    a=clusterCall(cl, function() library(cpt.paper))  
    a=clusterCall(cl, function() library(MASS)) 
    a=clusterCall(cl, function() library(crossmatch))  
    a=clusterCall(cl, function() library(energy)) 
    clusterSetRNGStream(cl, iseed=12345)
    MPs.analysis$analyses = parSapply(cl, MPs.analysis$windowsizes, do_analysis, Z=Z.all, tr=tr, wobl_margin=cpt.MPs$wobl_margin, perm.N=perm.N)
    stopCluster(cl)

save(MPs.analysis, file="MPs.rda")







