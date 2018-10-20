library(tidyverse)
library(scales)

n = 100
treat = as.factor(c(rep(0,n/2), rep(1,n/2)))
length.beta.vec = 8
beta.vec = seq(0,0.7,length=length.beta.vec)
length.rho.vec = 9
rho.vec = seq(0,0.8,length=length.rho.vec)


load("sim4d.rda")

sim = sim4d
x = beta.vec*10/7
xlab = expression(beta)
thetitle = ""
  results05 = apply(sim <= 0.05, c(1,2), mean)
  dp = as.data.frame(results05)
  dp$x = x
  df=gather(dp, value="power", key="test.type", -x)
  df$test.type = as.factor(df$test.type)
  df$test.type = factor(df$test.type, levels=c("combined", "ensemble", "glmnet", "forest", "energy", "crossmatch", "hotelling"))
  levels(df$test.type) <- c("CPT (Combined)","CPT (Ensemble)","CPT (Elastic Net)","CPT (Forest)","Energy","Cross-Match","Hotelling")
  df = df[!(df$test.typ=="Hotelling"),]
  df = df[!(df$test.typ=="Cross-Match"),]  
  p = ggplot(df,aes(y=power,x=x,col=test.type, linetype=test.type))+
  geom_line(size=1.2)+
  labs(
    title=thetitle,
    y = "Power",
    x = xlab)+
  theme_bw()+ 
  ylim(0,0.8)+
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
        legend.key.width=unit(5,"line")) +
  scale_color_manual(
      name=    c("CPT (Combined)","CPT (Ensemble)","CPT (Elastic Net)","CPT (Forest)","Energy"),
      values = hue_pal()(8)[c(1,2,3,5,6)]   ) +
  scale_linetype_manual(
      name=    c("CPT (Combined)","CPT (Ensemble)","CPT (Elastic Net)","CPT (Forest)","Energy"),
      values = c("solid", "22", "42", "44", "13", "1343", "73", "2262", "12223242",
                 "F282", "F4448444", "224282F2", "F1")[c(1,2,3,5,6)] )            
pdf("output/figures/sim4dpower.pdf")
print(p)
dev.off()


load("sim3d.rda")

sim = sim3d
x = beta.vec
xlab = expression(beta)
thetitle = ""
  results05 = apply(sim <= 0.05, c(1,2), mean)
  dp = as.data.frame(results05)
  dp$x = x
  df=gather(dp, value="power", key="test.type", -x)
  df$test.type = as.factor(df$test.type)
  df$test.type = factor(df$test.type, levels=c("combined", "ensemble", "glmnet", "forest", "energy", "crossmatch", "hotelling"))
  levels(df$test.type) <- c("CPT (Combined)","CPT (Ensemble)","CPT (Elastic Net)","CPT (Forest)","Energy","Cross-Match","Hotelling")
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
        legend.key.width=unit(5,"line")) +
  scale_color_manual(
      name=    c("CPT (Combined)","CPT (Ensemble)","CPT (Elastic Net)","CPT (Forest)","Energy","Cross-Match","Hotelling"),
      values = hue_pal()(8)[c(1,2,3,5,6,7,8)]   ) +
  scale_linetype_manual(
      name=    c("CPT (Combined)","CPT (Ensemble)","CPT (Elastic Net)","CPT (Forest)","Energy","Cross-Match","Hotelling"),
      values = c("solid", "22", "42", "44", "13", "1343", "73", "2262", "12223242",
                 "F282", "F4448444", "224282F2", "F1")[c(1,2,3,5,6,7,8)] )    
pdf("output/figures/sim3dpower.pdf")
print(p)
dev.off()




