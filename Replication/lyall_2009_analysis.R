########################################################
# Layell (2009) data application
# Date: March 19, 2017
#
########################################################

rm(list=ls())
set.seed(1297493)

# Check that relevent packages are installed
lib_list = c("ggplot2","cpt","xtable","BalanceCheck","dplyr","asbio","stargazer","rpart","ggdendro","crossmatch","Hotelling","nnet","randomForest","glmnet")
for (ii in  lib_list){
  cat("Working on package: ",ii, "\n \n")
  tmp = library(ii, logical.return = TRUE, character.only = TRUE)
  if (tmp==FALSE) {
    install.packages(ii)
    tmp = library(ii, logical.return = TRUE, character.only = TRUE)
    stopifnot(tmp == TRUE)
  }
}
tmp = library(cpt.paper, logical.return = TRUE)
if (tmp==FALSE){
  install.packages("/accounts/grad/shemtov/Documents/cpt/packages/cpt.paper_1.0.tar.gz", repos = NULL, type="source" )  
  library(cpt.paper)
}

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
perm = eval( parse(text=args[1]) )
cat("\n Number of permutations is: ",perm,"\n \n")

########################################
# Load data and data descriptives
########################################

data("cpt.violence")

cov = c("lpop2000","poverty","tariq","lelev","iso",
        "lnn","garrison","reb")
# checks:
stopifnot( sum( ! cpt.violence$rawdata[, cov] == cpt.violence$Z) == 0 )

# These covariates have been chosen to illustrate a scenario in which the CPT can detect imbalance in observed characteristics 
# that is not observed by either looking on the marginal distributions or by conducting a standard F-test.

### Descriptives

#### Balance in each covariate

f.stat = function(x1,x0){
  #sd.bias = summary(x1-x0)[4]/sd(c(x1,x0),na.rm=T)
  ks = ks.test(x1,x0)$p.value
  ttest = t.test(x1,x0)$p.value
  wilcox = wilcox.test(x1,x0)$p.value
  return(c(mean(x1,na.rm=T),mean(x0,na.rm=T),ttest,wilcox,ks))
}

tab = t(mapply(f.stat,as.list(data.frame(cpt.violence$Z[cpt.violence$T=="1",])),as.list(data.frame(cpt.violence$Z[cpt.violence$T=="0",]))))
tab = round(tab,dig=3)
colnames(tab) = c("Ave. Treat","Ave. control","T-test","Wilcoxon","KS")
rownames(tab) = c("log-Population","Poverty","Tariqa","log-Elevation","Isolation","log distance to Neighbor","Garrison","Rebel")
xtable(tab,caption="Balance table", dig=3)

# F-test

dp = as.data.frame(cpt.violence$Z)
dp$T = cpt.violence$T=="1"

ols1=lm(T~(.),data=dp)
summary(ols1)

########################################
# Estimation
########################################

violence.analysis = cpt(
  Z = cpt.violence$Z,
  T = cpt.violence$T,
  class.methods = "forest",
  metric = "logscore",
  perm.N = perm
)

violence.analysis2 = cpt(
  Z = cpt.violence$Z,
  T = cpt.violence$T,
  class.methods = "logistic2",
  comb.methods = "logistic2",
  metric = "logscore",
  perm.N = perm
)

########################################
# Figures Lyall (violance example)
########################################


### Random Forest

teststat.obs <- violence.analysis$teststat 
teststat.null <- violence.analysis$nulldist

range = c(
  min(teststat.obs, teststat.null),
  max(teststat.obs, teststat.null)
          )

fig<-ggplot(data.frame(statistic.null=teststat.null),aes(x=statistic.null))+
  xlim(range)+
  ylim(0,200)+
  geom_histogram(
    fill="grey",
    breaks=seq(range[1],range[2],length=40)
  )+
  labs(
    title="Classifier: Random forest \n",
    x="\n Test statistic under the null",
    y="\n Frequency"
  )+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank() )+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
  )+
  geom_vline(xintercept=teststat.obs,lty=1,col="black",size=1)+
  annotate("text", x = -.72, y = 160, label = "Observed test statistic",col="black", size = 5)
ggsave(fig, file = "output/figures/lyall-forest.pdf", width = 5, height=4)


### Logistic with interactions

teststat.obs <- violence.analysis2$teststat 
teststat.null <- violence.analysis2$nulldist

range = c(
  min(teststat.obs, teststat.null),
  max(teststat.obs, teststat.null)
          )

fig<-ggplot(data.frame(statistic.null=teststat.null),aes(x=statistic.null))+
  xlim(range)+
  ylim(0,200)+
  geom_histogram(
    fill="grey",
    breaks=seq(range[1],range[2],length=40)
  )+
  labs(
    title="Classifier: Logistic regression with interactions \n",
    x="\n Test statistic under the null",
    y="\n Frequency"
  )+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank() )+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
  )+
  geom_vline(xintercept=teststat.obs,lty=1,col="black",size=1)+
  annotate("text", x = -.578, y = 160, label = "Observed test statistic",col="black", size=5)

ggsave(fig, file = "output/figures/lyall-log2.pdf", width = 5, height=4)

##############################################
# Tree models: decomposing the difference in 
# joint distribution
##############################################

dp = data.frame(cpt.violence$Z, treat = as.numeric(cpt.violence$T=="1") )

rpart0 = rpart(treat~(.),data=dp)
rpart0 = prune(rpart0,cp=.025) # prune the tree for display purposes.
rpart.data <- dendro_data(rpart0)

dp1 = rpart.data$labels
dp1$label = as.character(dp1$label)
dp1$label = gsub("lnn","Log distance to Neighbor",dp1$label)
dp1$label = gsub("lelev","Log-Elevation",dp1$label)
dp1$label = gsub("lpop2000","Log-Population",dp1$label)
dp1$y = dp1$y + 0.005

dp2 = rpart.data$leaf_labels
dp2$y = dp2$y - 0.005
dp2$label = round(as.numeric(as.character(dp2$label)), dig=3)

fig = ggplot() +
  geom_segment(data=rpart.data$segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  #geom_point(data=rpart.data$segments, aes(x=xend, y=yend), size=4, alpha = 0.6, col = "grey")+
  geom_text(data=dp1, aes(x=x, y=y, label=label)) +
  geom_text(data=dp2, aes(x=x, y=y, label=label),size = 5) +
  theme_dendro()
ggsave(fig, file = "output/figures/fig_DIV_tree1.pdf", width = 5, height=4)


dp$nn.elev = dp$lelev*dp$lnn
# boxplot(dp$nn.elev~dp$treat)
# ks.test(dp$nn.elev[dp$treat==0],dp$nn.elev[dp$treat==1])
# t.test(dp$nn.elev[dp$treat==0],dp$nn.elev[dp$treat==1])

dp$nn.elev.pop = dp$lelev*dp$lnn*dp$lpop2000
# boxplot(dp$nn.elev.pop~dp$treat)
# ks.test(dp$nn.elev.pop[dp$treat==0],dp$nn.elev.pop[dp$treat==1])
# wilcox.test(dp$nn.elev.pop[dp$treat==0],dp$nn.elev.pop[dp$treat==1])
# t.test(dp$nn.elev.pop[dp$treat==0],dp$nn.elev.pop[dp$treat==1])

summary(lm1 <- glm(treat ~ lelev + lpop2000 +  lnn, data = dp, family = binomial(link ="logit")))
summary(lm2 <- glm(treat ~ nn.elev + lelev + lpop2000 +  lnn, data = dp, family = binomial(link ="logit")))
summary(lm3 <- glm(treat ~ nn.elev + nn.elev.pop + lelev + lpop2000 +  lnn, data = dp, family = binomial(link ="logit")))

stargazer(lm1,lm2,lm3,
          covariate.labels = c("Log distance to Neighbor $\\times$ Log-Elevation","Log distance to Neighbor   $\\times$ Log-Elevation  $\\times$ Log-Population",
                               "Log-Elevation","Log-Population","Log distance to Neighbor"))

##############################################
# Random-Forests variable importance figure: 
# decomposing the difference in joint distribution
##############################################

fit = randomForest(x = cpt.violence$Z, y = cpt.violence$T)
vi_rf = importance(fit)

dp = data.frame(
  name = as.factor(rownames(vi_rf)),
  MeanDecreaseGini = vi_rf[,1]
  )
levels(dp$name) = list(
  Garrison = "garrison", 
  Isolation = "iso",
  "log-Elevation" = "lelev",
  "log-Distance to neighbor" = "lnn",
   "log-Population" = "lpop2000",
    "Poverty" = "poverty",
    "Rebel" = "reb",
    "Tariqa" = "tariq")
stopifnot( sum(is.na(dp$name))==0 )

p <-ggplot(dp, aes(name, MeanDecreaseGini))+
  geom_bar(stat = "identity")+
  labs(x = "Covariate name",
         y = "Mean Decrease in Gini")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank() )+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)
  )
p = p+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, file= "output/figures/rf_variable_importance.pdf", width = 4, height = 5)

######################################
# Energy, Hotelling and cross-match tests
######################################

### Checking balance using multivariate tests:
violence.analysis2 = cpt.analysis(
  Z = cpt.violence$Z,
  T = cpt.violence$T,
  perm.N = 10,
  do_5fold = FALSE, 
  do_logistic2 = FALSE, 
  do_forest = FALSE
)

cat("Energy test P-value: ", violence.analysis2$energy$p.value,"\n")
cat("Crossmatch test P-value: ", violence.analysis2$crossmatch$approxpval,"\n")
cat("Hotelling test P-value: ",  
    hotelling.test(cpt.violence$Z[as.numeric(cpt.violence$T)==1,], 
                   cpt.violence$Z[as.numeric(cpt.violence$T)!=1,])$pval,"\n")





