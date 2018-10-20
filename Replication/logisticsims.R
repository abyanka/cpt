set.seed(1493241)
library(tidyverse)
library(ggplot2)
library(lmtest)

iterN = 10000
sampNs = c(25, 35,50,75,100,150,200)
results = matrix(NA,length(sampNs),4)
colnames(results) = c("sample size", "mean p0", "mean p1", "type1")
results[,1] = sampNs*2
for (j in 1:length(sampNs))
{
  nullp = nullp.cpt = rep(NA,iterN)
  sampN = sampNs[j]
  tr = rep(0:1, each=sampN)
  means = matrix(NA,iterN,2)
  for (i in 1:iterN)
  {
    nulldf = data.frame(matrix(rnorm(2*sampN*20),2*sampN,20))    
    glm.fit <-  glm(tr~(.),family=binomial(link="logit"),data=nulldf)
    glm.fit0 <- glm(tr~1,  family=binomial(link="logit"),data=nulldf)
    nullp[i] = lrtest(glm.fit0, glm.fit)[2,5]
    means[i,1] = mean(glm.fit$fit[1:sampN])
    means[i,2] = mean(glm.fit$fit[(sampN+1):(2*sampN)])
  }
  results[j,2:3]=apply(means,2,mean)
  results[j,4]=mean(nullp<0.05)
}

df = gather(data.frame(results), "variable", "fraction", -1)
varfactor = as.factor(df$variable)
levels(varfactor) = c("Mean P-score, control", "Mean P-score, treatment", "Type 1 error rate")
df$variable = varfactor
p = ggplot(df, aes(x=sample.size, y=fraction, 
                   color=variable, linetype=variable))+
      geom_line(size=1.5)+
      xlim(0,400)+
      ylim(0,1)+
      geom_hline(yintercept = 0.05,lwd=0.5,col="black",linetype="dotted")+
      ylab("") + 
      xlab("Sample Size") +
      labs(color="", linetype="")+
      scale_linetype_manual(values=c("dashed", "dashed", "solid"))+
      theme_bw() + 
      theme(legend.position="top")+
      theme(legend.key.width = unit(1.5,"cm"))+
      theme(axis.text = element_text(size=12))+
      theme(axis.title = element_text(size=16))
pdf("output/figures/logistic.pdf")
print(p)
dev.off()
