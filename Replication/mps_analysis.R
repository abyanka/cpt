##############################################################################################################
# Description: The program replicates all the figures, tables and results
# of the MPs for Sale: Eggers and Hainmueller (2009) data application
#
#
##############################################################################################################

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)


##########################################################
# RDD plot figure: CPT vs energy and crossmatch
##########################################################

load("MPs.rda")

### Function for extracting the results:

getresults = function(analysis)
{
  results = matrix(NA,length(analysis$windowsizes),5)
  colnames(results) = c("window", "observations", "CrossMatch", "Energy", "CPT")
  for (i in 1:length(analysis$windowsizes))
  {
    results[i,1] = analysis$windowsizes[i]
    results[i,2] = analysis$observations[i]
    results[i,3] = analysis$analyses["crossmatch",i]
    results[i,4] = analysis$analyses["energy",i]
    results[i,5] = analysis$analyses["forest",i]
  }
  return(results)
}

### MPs for sale: ###

dp = as.data.frame(getresults(MPs.analysis))
df=gather(dp[,-1], value = pv ,key = test.type, -observations )

p.cpt <- ggplot(df,aes(x=observations,y=pv,col=test.type, shape=test.type))+
  geom_point(size=2.5)+geom_line()+
  labs(
    title = "",
    y = "P-value \n" ,
    x = "Number of Obs. in window \n")+
  theme_bw()+ 
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0.05, lty=2,col="black")+
  scale_colour_grey()+
  guides(col=guide_legend(title="Test type: "),shape=guide_legend(title="Test type: "))+ # adding legend title
  theme(legend.position="bottom") # legend position

pdf("output/figures/MPs_rdd_window.pdf",width=8,height=7)
p.cpt+geom_vline(xintercept = 164, col="grey",lwd=1,lty=2)+
  annotate("text", x = 205, y = 1, 
           label = paste("EH choosen window")
           ,col="black")
dev.off()  


















