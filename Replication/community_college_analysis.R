
rm(list=ls())
set.seed(1297493)

library(cpt)
library(cpt.paper)

library(dplyr)
library(tidyr)
library(reshape)
library(ggplot2)
library(ggthemes)
library(rdd)
library(xtable)

data(list="cpt.college")

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
perm = eval( parse(text=args[1]) )
cat("\n Number of permutations is: ","\n \n")

# matching with replacement:
Z.rep = cpt.college$Z.matched.withreplacement
Z.norep = cpt.college$Z.matched.withoutreplacement
tr  = (cpt.college$T.matched.withreplacement=="1")*1
index = c(c(1:table(tr)[1]),c(1:table(tr)[1]))
  
perm.N=perm

# permutation matrix:
mat.perm = matrix(sample(c(0,1),size=length(unique(index))*perm.N,replace=TRUE),
                  nrow=length(unique(index)),ncol=perm.N)

fun_prob = function(prob, tstT) {
        indexmat = cbind(1:nrow(prob), tstT)
        return(mean(log(prob[indexmat]+.0001)))
      }

Z.rep=as.data.frame(Z.rep)
Z.norep=as.data.frame(Z.norep)

prob0.rep1  <- prob0.norep1 <- rep(NA,perm.N)
prob0.rep2  <- prob0.norep2 <- rep(NA,perm.N)

for (s in c(1:perm.N)){
  if (s%%20==0){cat("Iteration: ",s,"\n")}
  ## # permutation within pairs (blocks):
  tr0 = rep(NA,length(tr))
  for (i in unique(index)){
    tr0[index==i]=c(mat.perm[i,s],1-mat.perm[i,s])
  }
  prob0.rep1[s] = fun_prob(tstT = matrix(tr0,ncol=1), prob = matrix( glm(tr0~(.)^2,data=Z.rep,family=binomial(link="logit"))$fit, ncol = 1) )
  prob0.norep1[s] = fun_prob(tstT = matrix(tr0,ncol=1), prob = matrix( glm(tr0~(.)^2,data=Z.norep,family=binomial(link="logit"))$fit, ncol = 1) )
  
  ## # permutation across pairs:
  tr0 = sample(tr,size=length(tr),replace=FALSE)

  prob0.rep2[s] = fun_prob(tstT = matrix(tr0,ncol=1), prob = matrix( glm(tr0~(.)^2,data=Z.rep,family=binomial(link="logit"))$fit, ncol = 1) )
  prob0.norep2[s] = fun_prob(tstT = matrix(tr0,ncol=1), prob = matrix( glm(tr0~(.)^2,data=Z.norep,family=binomial(link="logit"))$fit, ncol = 1) )

}

# observed values:
prob.obs.rep = fun_prob(tstT = matrix(tr,ncol=1), matrix( glm(tr~(.)^2,data=Z.rep,family=binomial(link="logit"))$fit, ncol=1) )  
prob.obs.norep = fun_prob(tstT = matrix(tr,ncol=1), matrix( glm(tr~(.)^2,data=Z.norep,family=binomial(link="logit"))$fit, ncol=1) )

##############################################
### matching with replacement 
##############################################

dp.rep = data.frame(within.blocks=prob0.rep1,across.blocks = prob0.rep2)
dp.rep <- gather(dp.rep, key = Permutation.design, value = probability)

dp.rep$Permutation.design = as.factor(dp.rep$Permutation.design)
dp.rep$Permutation.design = dplyr::recode(dp.rep$Permutation.design, within.blocks = "Within blocks", across.blocks = "Across blocks"  )

fig.rep <- ggplot(dp.rep,aes(x=probability,fill=Permutation.design))+
  geom_density(alpha=0.75)+
  labs(
    title="Matching with replacement \n",
    y = " ",
    x = "Test statistic under the null hypothesis")+
  theme_bw()+ 
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  xlim( quantile(dp.rep$probability,0.999), quantile(dp.rep$probability,0.001) )+
  scale_fill_grey()+
  guides(fill=guide_legend(title="Randomization structure: "),shape=guide_legend(title="Test type: "))+ # adding legend title
  theme(legend.position="bottom")+  # legend position
  theme(axis.text.x = element_text(colour="grey20",size=14 ))+
  theme(axis.text.y = element_blank())+
  theme(axis.line.y = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(axis.title.y = element_text(colour="grey20",size=16))+
  theme(axis.title.x = element_text(colour="grey20",size=16))+
  geom_vline(xintercept=prob.obs.rep,lwd=1,col="black")+
  annotate("text", x = prob.obs.rep-0.02, y = 100, label = "Observed test statistic",col="black",size=5)

ggsave(plot = fig.rep, file="output/figures/college_rep_null_probability.pdf")


######################################################
### matching without replacement 
######################################################

dp.norep = data.frame(within.blocks=prob0.norep1,across.blocks = prob0.norep2)
dp.norep <- gather(dp.norep, key = Permutation.design, value = probability)

dp.norep$Permutation.design = as.factor(dp.norep$Permutation.design)
dp.norep$Permutation.design = dplyr::recode(dp.norep$Permutation.design, within.blocks = "Within blocks", across.blocks = "Across blocks"  )

fig.norep <- ggplot(dp.norep,aes(x=probability,fill=Permutation.design))+
  geom_density(alpha=0.75)+
  labs(
    title="Matching without replacement \n",
    y = " ",
    x = "Test statistic under the null hypothesis")+
  theme_bw()+ 
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_grey()+
  xlim(quantile(dp.norep$probability,0.999), quantile(dp.norep$probability,0.001) )+
  guides(fill=guide_legend(title="Randomization structure: "),shape=guide_legend(title="Test type: "))+ # adding legend title
  theme(legend.position="bottom")+  # legend position
  theme(axis.text.x = element_text(colour="grey20",size=14 ))+
  theme(axis.text.y = element_blank())+
  theme(axis.line.y = element_blank())+
  theme(axis.ticks.y = element_blank())+  
  theme(axis.title.y = element_text(colour="grey20",size=16))+
  theme(axis.title.x = element_text(colour="grey20",size=16))+
  geom_vline(xintercept=prob.obs.rep,lwd=1,col="black")+
  annotate("text", x = prob.obs.rep-0.02, y = 100, label = "Observed test statistic",col="black",size=5)

ggsave(plot = fig.norep, file="output/figures/college_norep_null_probability.pdf")

### Save results:

save.image(file="community_college_analysis.rda")


###############################################################
# Covariate balance diagnostics:
###############################################################

f.stat = function(x1,x0){
  #sd.bias = summary(x1-x0)[4]/sd(c(x1,x0),na.rm=T) # need to correct
  ks = ks.test(x1,x0)$p.value
  ttest = t.test(x1,x0)$p.value
  wilcox = wilcox.test(x1,x0)$p.value
  
  return(c(mean(x1,na.rm=T),mean(x0,na.rm=T),ttest,wilcox,ks))
}

plot.pval <- function(results, title=NULL, legend,legendx=0.7,legendy=3, textsize=0.9, parcex=0.8, at1=-0.35, at2=-0.15, at3=-0.9,xlim1=-0.85) {
    xlim = c(xlim1,1); pchset = c(21,24,22,23); pchcolset = c("gray76","darkgray","black","grey")
    par(cex=parcex, mai = c(0.5, 0.35, 1.1, 0.35))
    ny = nrow(results)
    if(!is.null(title))  plot(x=NULL,axes=F, xlim=xlim, ylim=c(1,ny),xlab="",ylab="", main=title)
    if(is.null(title))   plot(x=NULL,axes=F, xlim=xlim, ylim=c(1,ny),xlab="",ylab="")
    abline(v=c(0,0.05,0.1),lty=c(1,4,4), lwd=c(1,2,2))
    axis(side=1,at=c(0,0.05,0.1,1),tick=TRUE, las=2, cex.axis=0.7)
    axis(side=3,at=at1,labels="Mean\nTreated",tick=FALSE, padj=0.5,cex.axis=textsize)
    axis(side=3,at=at2,labels="Mean\nControl",tick=FALSE, padj=0.5,cex.axis=textsize)
    axis(side=3,at=0.5,labels="P-values",tick=FALSE, padj=0.5,cex.axis=textsize)
    for(i in 4:ncol(results)) points(results[,i],ny:1, pch = pchset[i-4+1], col = pchcolset[i-4+1], bg = pchcolset[i-4+1])
    for(i in 1:ny) {
      text(at3,ny-i+1,results[i,1],adj = 0,cex=textsize) # variable name
      text(at1,ny-i+1,results[i,2], cex=textsize)        # treatment mean
      text(at2,ny-i+1,results[i,3], cex=textsize)        # control mean
    }
    for(i in seq(2,by=2,length.out=floor((ny-1)/2))) abline(h = i+0.5, lty = 3)
    if(legend) legend(x=legendx, y=legendy, c("T-test","Wilcoxon","KS"), pch=pchset, pt.bg = pchcolset, cex=0.8)
}
  
### Prior to matching
x.unmatch = as.data.frame(cpt.college$Z.unmatched)
tr.unmatch = as.data.frame(cpt.college$T.unmatched)

tab=round(t(mapply(f.stat,as.list(x.unmatch[tr.unmatch=="1",]),
                   as.list(x.unmatch[tr.unmatch=="0",]))),dig=3)
tab=apply(tab,2,function(x){ return( ifelse(x>1,as.integer(round(x,dig=0)),x)) }  )
tab = data.frame(name=rownames(tab),tab)


pdf(file="output/figures/college_balance_no_match.pdf")
plot.pval(tab,legend=TRUE,title="Prior to matching \n \n")
dev.off()

### After to matching without replacement
x.match.no.rep = as.data.frame(cpt.college$Z.matched.withoutreplacement)
tr.match.no.rep = as.data.frame(cpt.college$T.matched.withoutreplacement)

tab=round(t(mapply(f.stat,as.list(x.match.no.rep[tr.match.no.rep=="1",]),
                   as.list(x.match.no.rep[tr.match.no.rep=="0",]))),dig=3)
tab=apply(tab,2,function(x){ return( ifelse(x>1,as.integer(round(x,dig=0)),x)) }  )
tab = data.frame(name=rownames(tab),tab)


pdf(file="output/figures/college_balance_match_norep.pdf")
plot.pval(tab,legend=TRUE,title="After matching without replacement \n \n")
dev.off()

### After to matching with replacement
x.match.rep = as.data.frame(cpt.college$Z.matched.withreplacement)
tr.match.rep = as.data.frame(cpt.college$T.matched.withreplacement)

tab=round(t(mapply(f.stat,as.list(x.match.rep[tr.match.rep=="1",]),
                   as.list(x.match.rep[tr.match.rep=="0",]))),dig=3)
tab=apply(tab,2,function(x){ return( ifelse(x>1,as.integer(round(x,dig=0)),x)) }  )
tab = data.frame(name=rownames(tab),tab)

pdf(file="output/figures/college_balance_match_rep.pdf")
plot.pval(tab,legend=TRUE,title="After matching with replacement \n \n")
dev.off()











