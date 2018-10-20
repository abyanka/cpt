
rm(list=ls())
set.seed(1493241)

library(cpt)
library(cpt.paper)

library(ggplot2)
library(lmtest)
library(nnet)

library(xtable)
library(stargazer)

library(dplyr)
library(tidyr)

library(randomForest)

library(energy)
library(crossmatch)

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
perm = eval( parse(text=args[1]) )
cat("\n Number of permutations is: ",perm,"\n \n")

data(list="cpt.judges")
labels = cbind(attributes(cpt.judges$rawdata)$var.labels,attributes(cpt.judges$rawdata)$names)

Z = cpt.judges$Z
Z = Z[,!colnames(Z) %in% c("pcp","female")]
tr  = cpt.judges$T
perm.N=perm

##################################################################
### Multinomial logit and liklihood ratio test
##################################################################

multi.fit <- multinom(tr ~ (.), data = data.frame(Z,tr))
multi.fit0 <- multinom(tr ~ 1, data = data.frame(Z,tr))

anova(multi.fit0,multi.fit)

####################################################################################################################################
### The Energy and Cross-Match tests for whether defendants have been randomly allocated to judge calendars
####################################################################################################################################

# Energy
energy.test.results = eqdist.etest(Z, sizes = tapply(tr,tr,length), R=999)
print(energy.test.results)

# The cross-match test is designed for a binary vector - see ?crossmatchtest


##################################################################
### Assignment to all Judge calendars
##################################################################

# CPT ( logistic )
cpt.result0 <- cpt(Z=Z,T= tr, class.methods = "logistic", metric = "logscore", perm.N = perm)
print( cpt.result0$pval )


# CPT (forest )
cpt.result1 <- cpt(Z=Z,T= tr, class.methods = "forest", metric = "logscore",  perm.N = perm)
print( cpt.result1$pval )



##################################################################
### Assignment to a specific judge
##################################################################

pv.lrt <- pv.lrt2 <- coef.number <- coef.number2 <- rep(NA,9) 
for (i in c(1:9)){
  tr.judge <- (tr==paste(i))*1
  
  glm.fit <- glm(tr.judge~(.),family=binomial(link="logit"),data=data.frame(Z))
  coef.number <- length(coef(glm.fit)[!is.na(coef(glm.fit))])
  
  glm.fit2 <- glm(tr.judge~(.)^2,family=binomial(link="logit"),data=data.frame(Z))
  coef.number2 <- length(coef(glm.fit2)[!is.na(coef(glm.fit2))])
  pv.lrt[i] <-  lrtest(glm.fit,glm(tr.judge~1,family=binomial(link="logit"),data=data.frame(Z)))$P[2]
  pv.lrt2[i] <- lrtest(glm.fit2,glm(tr.judge~1,family=binomial(link="logit"),data=data.frame(Z)))$P[2]
}

### Permutation distribution of the LRT:

mat.pv2 <- mat.pv <- matrix(NA,ncol=9,nrow=perm.N)

for (i in c(1:9)){
  cat("Calendar: ",i,"\n")  
  for (s in c(1:perm.N)){
    if(s%%50==0){cat("Iteration: ",s,"\n")}
    tr0 <- sample(tr,size=length(tr),replace=FALSE)
    tr.judge <- (tr0==paste(i))*1
    
    glm.fit.null <- glm(tr.judge~(.),family=binomial(link="logit"),data=data.frame(Z))
    mat.pv[s,i] <- lrtest(glm.fit.null,glm(tr.judge~1,family=binomial(link="logit"),data=data.frame(Z)))$P[2]
    
    glm.fit2.null <- glm(tr.judge~(.)^2,family=binomial(link="logit"),data=data.frame(Z))
    mat.pv2[s,i] <- lrtest(glm.fit2.null,glm(tr.judge~1,family=binomial(link="logit"),data=data.frame(Z)))$P[2]
  }
}

### Summarizing table:
tab <- rbind(pv.lrt,coef.number,matrix(apply(mat.pv,2,function(x){sum(x<0.05)/dim(mat.pv)[1]}),ncol=9),
             pv.lrt2,coef.number2,matrix(apply(mat.pv2,2,function(x){sum(x<0.05)/dim(mat.pv)[1]}),ncol=9))
colnames(tab) <- c(1:9)
rownames(tab) <- c(c("P-value","Number of coefficients","Type-I error rate"),
                   c("P-value","Number of coefficients","Type-I error rate"))
tab <- t(tab)
tab <- round(tab,dig=3)
#rownames(tab) <- c("Main effects only", "All two-way interactions")

stargazer(tab,title="The Likelihood Ratio Test P-values and Type-I error rates for each judge calendar dummy",
          label="tab: type-I error rate judge calendar")

######################################################
# Illustration figure of Type-I error rate
######################################################

dp <- gather(data.frame(mat.pv2),key ="calendar")
dp$calendar = as.factor(dp$calendar)
levels(dp$calendar) <- paste("Judge calendar",c(1:9))
dp <- dplyr::rename(dp,  P.value = value)

fig_type1 <- ggplot(dp,aes(x=P.value))+
  geom_density(fill="grey", bw=.03)+
  facet_wrap(~ calendar)+
  labs(
  title="",
  y = " ",
  x = "\n LRT P-value under the null hypothesis")+
  theme_bw()+ 
  ylim(0,10)+
  theme( panel.border = element_blank(),
         axis.line.x = element_line(color = "black", size = 0.5),
         axis.line.y = element_line(color = "black", size = 0.5 ),
         axis.title.x = element_text( size = 12),
         axis.title.y = element_text(size = 12 ),
         axis.text.x = element_text( size = 14),
         axis.text.y = element_text(size = 12 )
  )+
  scale_fill_grey()+
  theme(strip.background = element_rect(fill="white"))+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  geom_hline(yintercept  = 1,lwd=0.5,col="black")

ggsave(fig_type1, file = "output/figures/type1_judge_lrt.pdf")


##################################################################
### The case of Judge calendar 1 and over-fitting
##################################################################

tr.judge <- (tr==1)*1
glm.fit <- glm(tr.judge~(.),family=binomial(link="logit"),data=data.frame(Z))
glm.fit2 <- glm(tr.judge~(.)^2,family=binomial(link="logit"),data=data.frame(Z))

# CPT
cpt.result <- cpt(Z=Z,T= tr.judge, class.methods = "logistic", metric = "logscore")
print( cpt.result$pval )

cpt.result2 <- cpt(Z=Z,T= tr.judge, class.methods = "logistic2", metric = "logscore")
print( cpt.result2$pval )

dp = data.frame(ps=glm.fit$fit,tr=factor(tr.judge, levels = c(0,1), labels = c("Not Judge 1","Judge 1")))
fig.box1 <- ggplot(dp,aes(x=ps,fill=tr))+
  geom_density(alpha=0.75)+
  labs(
    title="Main effects only",
    y = "",
    x = "P-score")+
  theme_bw()+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+  
  scale_fill_grey()+  
  guides(fill=guide_legend(title="Calendar: "))+  
  theme( panel.border = element_blank(),
         axis.line.x = element_line(color = "black", size = 0.5),
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank())+
  theme(axis.text.x = element_text(colour="grey20",size=14 ))+
  theme(axis.text.y = element_blank())+
  theme(axis.title.y = element_text(colour="grey20",size=16))+
  theme(axis.title.x = element_text(colour="grey20",size=16))+         
  theme(legend.position="bottom") # legend position
  
ggsave(fig.box1, file = "output/figures/fig_pscore_overfit1.pdf")

####

dp = data.frame(ps=glm.fit2$fit,tr=factor(tr.judge, levels = c(0,1), labels = c("Not Judge 1","Judge 1")))
fig.box2 <- ggplot(dp,aes(x=ps,fill=tr))+
  geom_density(alpha=0.75)+
  labs(
    title="All two-way interactions",
    y = "",
    x = "P-score")+
  theme_bw()+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+  
  scale_fill_grey()+  
  guides(fill=guide_legend(title="Calendar: "))+  
  theme( panel.border = element_blank(),
         axis.line.x = element_line(color = "black", size = 0.5),
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank())+
  theme(axis.text.x = element_text(colour="grey20",size=14 ))+
  theme(axis.text.y = element_blank())+
  theme(axis.title.y = element_text(colour="grey20",size=16))+
  theme(axis.title.x = element_text(colour="grey20",size=16))+         
  theme(legend.position="bottom") # legend position

ggsave(fig.box2, file = "output/figures/fig_pscore_overfit2.pdf")

save.image(file = "green_winik_2010_results.rda")

##################################################################
### Continuous treatment
##################################################################
rmse = function(x1,x2){
  return(sqrt(mean((x1-x2)^2)))
}

Z = as.data.frame(Z)
Z$incarcerate = cpt.judges$rawdata$incarcerate
rmse.vec <- matrix(NA,ncol=1,nrow=perm.N)

# Observed values
Z$tr = cpt.judges$T
Z = Z %>% group_by( tr ) %>%
  mutate(
    obs = n(),
    tot = sum(incarcerate),
    loo = (tot - incarcerate)/(obs - 1)
  ) %>%
  ungroup() %>%
  select(-tot,-obs,-incarcerate, -tr)

rf = randomForest(x = select(Z,-loo), y = Z$loo)
rmse.obs = rmse( Z$loo, predict(rf, newdata = select(Z,-loo)) )

for (s in c(1:perm.N)){
  if(s%%50==0){cat("Iteration: ",s,"\n")}
  Z$loo0 =  sample(Z$loo, length(Z$loo), replace = FALSE)
  rf = randomForest(x = select(Z,-loo,-loo0), y = Z$loo0)
  rmse.vec[s] = rmse( Z$loo0, predict(rf, newdata = select(Z,-loo,-loo0)) )
}

xmax = max(c(rmse.obs,rmse.vec)) + min(rmse.vec)/20
xmin = min(c(rmse.obs,rmse.vec)) - max(rmse.vec)/20
p <- ggplot(data.frame(rmse = rmse.vec),aes(x=rmse))+
  xlim(0.081, 0.1)+
  ylim(0,350)+
#  geom_histogram(
#    fill="grey"
#  )+
  geom_density(alpha=0.75, fill="grey")+
  xlim(xmin,xmax)+
  labs(
    title="",
    y = " ",
    x = "\n RMSE")+
  geom_vline(xintercept = rmse.obs, col = "black", lwd=1)+
  theme_bw()+ 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),           
        axis.title = element_text(size = 12))+
  annotate("text", x = 0.0966, y = 320, label = "Observed test statistic",col="black", size=5)        

ggsave(p, file = "output/figures/rf_null.pdf")

save.image(file = "green_winik_2010_results_rf_rmse.rda")



















