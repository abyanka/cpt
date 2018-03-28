library(tidyverse)

n = 100
treat = as.factor(c(rep(0,n/2), rep(1,n/2)))
length.beta.vec = 8
beta.vec = seq(0,0.7,length=length.beta.vec)
length.rho.vec = 9
rho.vec = seq(0,0.8,length=length.rho.vec)

load("sim1.rda")
load("sim1d.rda")
load("sim2.rda")
load("sim2d.rda")

powerplot = function(sim, x, xlab, thetitle="Significance level = 0.05")
{
  results05 = apply(sim <= 0.05, c(1,2), mean)
  dp = as.data.frame(results05)
  dp$x = x
  df=gather(dp, value="power", key="test.type", -x)
  df$test.type = as.factor(df$test.type)
  df$test.type = factor(df$test.type, levels=c("combined", "ensemble", "logistic", "logistic2", "forest", "energy", "crossmatch", "hotelling"))
  levels(df$test.type) <- c("CPT (Combined)","CPT (Ensemble)","CPT (Logistic)","CPT (Logistic2)","CPT (Forest)","Energy","Cross-Match","Hotelling")
  p = ggplot(df,aes(y=power,x=x,col=test.type, linetype=test.type))+
  geom_line(size=1.2)+
  labs(
    title=thetitle,
    y = "Power",
    x = xlab)+
  theme_bw()+ 
  theme(axis.text.x = element_text(colour="grey20",size=16,angle=90,hjust=.5,vjust=.5))+
  theme(axis.text.y = element_text(colour="grey20",size=16))+
  theme(axis.title.y = element_text(colour="grey20",size=17))+
  theme(axis.title.x = element_text(colour="grey20",size=22))+
  theme(panel.border = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.8),
        axis.line.y = element_line(color="black", size = 0.8))+
  theme(title = element_text(size=15))+
  guides(col=guide_legend(title=""),linetype=guide_legend(title=""))+ 
  theme(legend.position=c(.2,.82),
        legend.text = element_text(size=11),
        legend.key.width=unit(5,"line")) 
  return(p)
}

p = powerplot(sim1d, beta.vec, expression(beta), thetitle="")
pdf("output/figures/sim1dpower.pdf")
print(p)
dev.off()

p = powerplot(sim2d, rho.vec, expression(rho), thetitle="")
pdf("output/figures/sim2dpower.pdf")
print(p)
dev.off()

metricplot = function(sima, simd, x, xlab)
{
  thetitle="Significance level = 0.05"
  results = rbind(
      data.frame(metric="rate"    , x=x, apply(sima <= 0.05, c(1,2), mean)[,1:3]),
      data.frame(metric="logscore", x=x, apply(simd <= 0.05, c(1,2), mean)[,1:3])      
  )
  df = gather(results, value=power, key="test.type", -metric, -x)
  p = ggplot(df, aes(x=x, y=power)) + 
    geom_line(size=1.2) + 
    aes(color=metric) +
    aes(shape=metric) +
    labs(
      y = "Power",
      x = xlab)+
    theme_bw()+ 
    facet_wrap(~test.type) +
    theme(axis.text.x = element_text(colour="grey20",size=16,angle=90,hjust=.5,vjust=.5))+
    theme(axis.text.y = element_text(colour="grey20",size=16))+
    theme(axis.title.y = element_text(colour="grey20",size=17))+
    theme(axis.title.x = element_text(colour="grey20",size=22))+
    guides(col=guide_legend(title=""),linetype=guide_legend(title="")) +
    theme(legend.text = element_text(size=15)) +
    theme(strip.text.x = element_text(size = 12))
  return(p)
}

p = metricplot(sim1, sim1d, beta.vec, expression(beta))
pdf("output/figures/sim1metric.pdf", width=12, height=7)
print(p)
dev.off()

p = metricplot(sim2, sim2d, rho.vec, expression(rho))
pdf("output/figures/sim2metric.pdf", width=12, height=7)
print(p)
dev.off()
