## Estimation of Number of Interventions on Ecosystem Service Multifunctionality 
## Using Mixed Effects Model and Visualization


library(nlme)
library(lme4)
library(lmerTest)
library(stringr)
library(xtable)
library(ggplot2)
library(gridExtra)


## Load the data
test.new <- read.csv("mfdata_combined.csv", header = T, stringsAsFactors = F)
head(test.new)
names(test.new)

test.new$Block <- factor(test.new$Block)
test.new$Tetrad <- factor(test.new$Tetrad)
test.new$Seed <- factor(test.new$Seed, labels=c("absence", "with"))
test.new$Clover <- factor(test.new$Clover, labels=c("absence", "with"))
test.new$NPK <- factor(test.new$NPK, labels=c("absence", "with"))
test.new$FYM <- factor(test.new$FYM, labels=c("absence", "with"))

Nooftreat <- apply(sapply(test.new[3:6], FUN=as.numeric), MARGIN=1, FUN=sum) - 4
test.new$Nooftreatment <- Nooftreat

## List of response and explanatory variables
all.y <- names(test.new)[17:42]           # 26 metrics
all.index <- names(test.new)[10:16]       # 6 combined indices
all.treat <- names(test.new)[3:6]         # 4 treatments



## --------------------------------------------------
## (I) Model the impact of number of treatment
## --------------------------------------------------
## Three types of models

## (1) The baseline model
# The goal is to estimate the effect of number of treatment on the index
# There is supposed to be random variation brought by the block and 
# the different treatment combinations.
# At this stage, we ignore the exact content of the combinations. 
# We only assume that there are variations caused by them and 
# we wish to account for that so that the impact of 
# the number of treatments can be estimated more accurately

lmer0 <- lmer(eq.mf ~ Nooftreatment + (1|Block) + (1|Treatment), 
              data=test.new, REML=TRUE)
summary(lmer0)

sum0 <- summary(lmer0)
var.re <- sum(unlist(sum0$varcor))   # random effect variance 
var.resi <- sum0$sigma^2             # residual variance
var.re / (var.re + var.resi)         # proportion of variance over total variance
sum0$coefficients[2, 1] / sqrt(var.re + var.resi)     # fixed effect over total variance


## (2) The model for investigating the impact from the treatments
# Here we want to investigate the influence from the number of treatments,
# and whether they are a result of different treatments compliment each other
# or a result of a few high impact treatments dominating the effect.
# Therefore, the random effect of the treatment combinations is replaced by
# the random effect of the treatments themselves. 

lmer1 <- lmer(eq.mf ~ Nooftreatment + (1|Block) + 
                (1|NPK) + (1|FYM) + (1|Seed) + (1|Clover), 
              data=test.new, REML=TRUE)
# control = lmerControl(optimizer ="Nelder_Mead")
summary(lmer1)

sum1 <- summary(lmer1)
unlist(sum1$varcor) > 1e-4    # check which one has a larger impact
var.re <- sum(unlist(sum1$varcor))
var.resi <- sum1$sigma^2
var.re / (var.re + var.resi)
sum1$coefficients[2, 1] / sqrt(var.re + var.resi)


## (3) The model for investigating the impact from the interaction of treatments
# Here we further investigate the impact of the interactions of treatments

lmer2 <- lmer(eq.mf ~ Nooftreatment + (1|Block) + 
                (1|NPK) + (1|FYM) + (1|Seed) + (1|Clover) + 
                (1|NPK:FYM) + (1|NPK:Seed) + (1|NPK:Clover) + (1|FYM:Seed) + 
                (1|FYM:Clover) + (1|Seed:Clover), 
              data=test.new, REML=TRUE)
summary(lmer2)

sum2 <- summary(lmer2)
unlist(sum2$varcor) > 1e-4
var.re <- sum(unlist(sum2$varcor))
var.resi <- sum2$sigma^2
var.re / (var.re + var.resi)
sum2$coefficients[2, 1] / sqrt(var.re + var.resi)



## ----------------------------------------
## (II) Fit mixed model for all indices
## ----------------------------------------
index.List <- vector(mode="list", length=length(all.index))    # a list for storing the models
random.List <- vector(mode="list", length=length(all.index))   # contribution of the random effect
effect.List <- vector(mode="list", length=length(all.index))   # effect of fixed and random effect, AIC
fstats.mat <- matrix(NA, nrow=2, ncol=length(all.index))       # The F-statistics

for (i in 1:length(all.index)) {
  yname <- all.index[i] 
  
  ## the baseline model
  formula0 <- formula(paste(yname, "~ Nooftreatment + (1|Block) + (1|Treatment)"))
  lmer0 <- lmer(formula0, data=test.new, REML=TRUE)
  # influence from the fixed and random effect
  sum0 <- summary(lmer0)
  lmer0.rand <- unlist(sum0$varcor)
  var.re <- sum(unlist(sum0$varcor))    # variance of each random component
  var.resi <- sum0$sigma^2
  lmer0.var <- var.re + var.resi        # variance not explained by the fixed effect
  lmer0.re <- var.re / (var.re + var.resi)        # contribution from the random effect
  lmer0.coef <- sum0$coefficients[2, 1]
  lmer0.fe <- sum0$coefficients[2, 1] / sqrt(var.re + var.resi)       # effect size of fixed effect
  lmer0.aic <- sum0$AICtab
  lmer0.p <- sum0$coefficients[2, 5]
  
  fstats.mat[1, i] <- anova(lmer0)$`F value`
  
  ## the model with treatment identity
  formula1 <- formula(paste(yname, "~ Nooftreatment + (1|Block) +  
                            (1|Seed) + (1|Clover) + (1|NPK) + (1|FYM)"))
  lmer1 <- lmer(formula1, data=test.new, REML=TRUE)
  
  sum1 <- summary(lmer1)
  lmer1.rand <- unlist(sum1$varcor)   
  var.re <- sum(unlist(sum1$varcor))
  var.resi <- sum1$sigma^2
  lmer1.var <- var.re + var.resi
  lmer1.re <- var.re / (var.re + var.resi)
  lmer1.coef <- sum1$coefficients[2, 1]
  lmer1.fe <- sum1$coefficients[2, 1] / sqrt(var.re + var.resi)
  lmer1.aic <- sum1$AICtab
  lmer1.p <- sum1$coefficients[2, 5]
  
  fstats.mat[2, i] <- anova(lmer1)$`F value`
  
  ## the model with treatments interaction
  formula2 <- formula(paste(yname, "~ Nooftreatment + (1|Block) + 
                            (1|Seed) + (1|Clover) + (1|NPK) + (1|FYM) + 
                            (1|NPK:FYM) + (1|NPK:Seed) + (1|NPK:Clover) + 
                            (1|FYM:Seed) + (1|FYM:Clover) + (1|Seed:Clover)"))
  lmer2 <- lmer(formula2, data=test.new, REML=TRUE)
  
  sum2 <- summary(lmer2)
  lmer2.rand <- unlist(sum2$varcor)
  var.re <- sum(unlist(sum2$varcor))
  var.resi <- sum2$sigma^2
  lmer2.var <- var.re + var.resi
  lmer2.re <- var.re / (var.re + var.resi)
  lmer2.coef <- sum2$coefficients[2, 1]
  lmer2.fe <- sum2$coefficients[2, 1] / sqrt(var.re + var.resi)
  lmer2.aic <- sum2$AICtab
  lmer2.p <- sum2$coefficients[2, 5]
  
  ## save the result
  index.List[[i]] <- list(lmer0=lmer0, lmer1=lmer1, lmer2=lmer2)
  random.List[[i]] <- list(lmer0=lmer0.rand, lmer1=lmer1.rand, lmer2=lmer2.rand)
  effect.term <- cbind(c(lmer0.var, lmer1.var, lmer2.var),
                       c(lmer0.re, lmer1.re, lmer2.re),
                       c(lmer0.coef, lmer1.coef, lmer2.coef),
                       c(lmer0.fe, lmer1.fe, lmer2.fe),
                       c(lmer0.p, lmer1.p, lmer2.p),
                       c(lmer0.aic, lmer1.aic, lmer2.aic))
  colnames(effect.term) <- c("variation", "re.size", "nooftreatment", "fe.size", "pvalue",  "aic")
  rownames(effect.term) <- c("lmer0", "lmer1", "lmer2")
  effect.List[[i]] <- effect.term
}


## Extract information on the contribution of different model components
# summat.List <- vector(mode="list", length=length(all.index))
# for (i in 1:length(all.index)) {
#   basic <- effect.List[[i]][, c(3, 5)]
#   aic <- effect.List[[i]][, 6]
#   re0 <- c(NA, NA)
#   identity1.var <- sum(random.List[[i]]$lmer1[-1])
#   re1 <- c(identity1.var / effect.List[[i]][2, 1], NA)
#   identity2.var <- sum(random.List[[i]]$lmer2[8:11])
#   interact2.var <- sum(random.List[[i]]$lmer2[1:6])
#   re2 <- c(identity2.var / effect.List[[i]][3, 1], 
#            interact2.var / effect.List[[i]][3, 1])
#   re <- rbind(re0, re1, re2)
#   colnames(re) <- c("identity", "interact")
#   
#   summat.List[[i]] <- cbind(basic, re, aic)
# }
# names(summat.List) <- all.index
# The results above is for comparing the three types of models


## Extract the coefficient of nooftreatment (for reporting)
coef0 <- sapply(1:length(all.index), FUN=function(x) {effect.List[[x]][1, 3]})
coef1 <- sapply(1:length(all.index), FUN=function(x) {effect.List[[x]][2, 3]})
pvalue0 <- sapply(1:length(all.index), FUN=function(x) {effect.List[[x]][1, 5]})
pvalue1 <- sapply(1:length(all.index), FUN=function(x) {effect.List[[x]][2, 5]})

coef.notreat.index <- rbind(coef0, coef1)
pvalue.notreat.index <- rbind(pvalue0, pvalue1)
colnames(coef.notreat.index) <- colnames(pvalue.notreat.index) <- all.index
rownames(coef.notreat.index) <- rownames(pvalue.notreat.index) <- c("Nooftreatment_baseline", "Nooftreatment_mixed")


## F-statistic of Nooftreatment
fstats.index <- matrix(NA, nrow=2, ncol=length(all.index))
for (i in 1:length(all.index)) {
  lmer0 <- index.List[[i]]$lmer0
  lmer1 <- index.List[[i]]$lmer1
  aov0 <- anova(lmer0)
  aov1 <- anova(lmer1)
  fstats.index[1, i] <- aov0$`F value`
  fstats.index[2, i] <- aov1$`F value`
}
rownames(fstats.index) <- c("Nooftreatment_baseline", "Nooftreatment_mixed")
colnames(fstats.index) <- all.index



## ------------------------------------------
## (III) Fit mixed model for all 26 metrics
## ------------------------------------------
ally.List <- vector(mode="list", length=length(all.y))    # a list for storing the models
ally.random.List <- vector(mode="list", length=length(all.y))   # contribution of the random effect
ally.effect.List <- vector(mode="list", length=length(all.y))   # effect of fixed and random effect, AIC

for (i in 1:length(all.y)) {
  yname <- all.y[i] 
  
  ## the baseline model
  formula0 <- formula(paste(yname, "~ Nooftreatment + (1|Block) + (1|Treatment)"))
  lmer0 <- lmer(formula0, data=test.new, REML=TRUE)
  # influence from the fixed and random effect
  sum0 <- summary(lmer0)
  var.re <- sum(unlist(sum0$varcor))
  var.resi <- sum0$sigma^2
  lmer0.re <- var.re / (var.re + var.resi)
  lmer0.notreat <- sum0$coefficients[2, 1] / sqrt(var.re + var.resi)
  lmer0.aic <- sum0$AICtab
  lmer0.p <- sum0$coefficients[2, 5]
  
  ## the model with treatment identity
  formula1 <- formula(paste(yname, "~ Nooftreatment + (1|Block) +  
                            (1|Seed) + (1|Clover) + (1|NPK) + (1|FYM)"))
  lmer1 <- lmer(formula1, data=test.new, REML=TRUE)
  
  sum1 <- summary(lmer1)
  lmer1.large <- unlist(sum1$varcor) > 1e-4   
  # need to pick a threshold, maybe 0.1 * sum1$sigma^2?
  var.re <- sum(unlist(sum1$varcor))
  var.resi <- sum1$sigma^2
  lmer1.re <- var.re / (var.re + var.resi)
  lmer1.notreat <- sum1$coefficients[2, 1] / sqrt(var.re + var.resi)
  lmer1.aic <- sum1$AICtab
  lmer1.p <- sum1$coefficients[2, 5]
  
  ## save the result
  ally.List[[i]] <- list(lmer0=lmer0, lmer1=lmer1)
  ally.random.List[[i]] <- list(lmer1=lmer1.large)
  ally.effect.term <- cbind(c(lmer0.re, lmer1.re),
                            c(lmer0.notreat, lmer1.notreat),
                            c(lmer0.p, lmer1.p),
                            c(lmer0.aic, lmer1.aic))
  colnames(ally.effect.term) <- c("random", "notreatment", "pvalue",  "aic")
  rownames(ally.effect.term) <- c("lmer0", "lmer1")
  ally.effect.List[[i]] <- ally.effect.term
  
}

# save(index.List, random.List, effect.List, ally.List, ally.random.List, 
#      ally.effect.List, file="Colt_park_model2.RData")


## Extract the coefficient of number of treatment (for reporting)
coef.notreat.y <- sapply(1:length(all.y), FUN=function(x) { ally.effect.List[[x]][, 2] })
pvalue.notreat.y <- sapply(1:length(all.y), FUN=function(x) { ally.effect.List[[x]][, 3] })
colnames(coef.notreat.y) <- colnames(pvalue.notreat.y) <- all.y
rownames(coef.notreat.y) <- rownames(pvalue.notreat.y) <- c("Nooftreatment_baseline", "Nooftreatment_mixed")

coef.notreat.extend <- cbind(coef.notreat.y, coef.notreat.index)
pvalue.notreat.extend <- cbind(pvalue.notreat.y, pvalue.notreat.index)


## F-statistic of Nooftreatment
fstats.ally <- matrix(NA, nrow=2, ncol=length(all.y))
for (i in 1:length(all.y)) {
  lmer0 <- ally.List[[i]]$lmer0
  lmer1 <- ally.List[[i]]$lmer1
  aov0 <- anova(lmer0)
  aov1 <- anova(lmer1)
  fstats.ally[1, i] <- aov0$`F value`
  fstats.ally[2, i] <- aov1$`F value`
}
rownames(fstats.ally) <- c("Nooftreatment_baseline", "Nooftreatment_mixed")
colnames(fstats.ally) <- all.y

fstats.notreat.extend <- cbind(fstats.ally, fstats.index)

# save(coef.notreat.extend, pvalue.notreat.extend, fstats.notreat.extend,
#      file="Modelling_result2.RData")


## Finally, create the following matrices for visualization
# First run the code on R_script_3 to get the coefficient and P-value matrices
# If coef.mat and pvalue.mat have been saved, then load the result
load("Modelling_result.RData")
coef.mat.extend <- rbind(coef.mat, coef.notreat.y)
pvalue.mat.extend <- rbind(pvalue.mat, pvalue.notreat.y)
rownames(coef.mat.extend) <- c(all.terms, "Nooftreatment_baseline", "Nooftreatment_mixed")
rownames(pvalue.mat.extend) <- c(all.terms, "Nooftreatment_baseline", "Nooftreatment_mixed")



## -----------------------------------------------------
## (IV) Effect of increasing the number of treatments
## -----------------------------------------------------
## Three types of effect may exist as the number of treatments increases.
## They are: sampling effect, synergistic effect and antagonistic effect
## Maybe able to distinguish them by comparing the observed effect 
## to the predicted dominant effect and additive effect
## The code below are motivated by the method in Rillig et al (2019)

## create a new variable called "Treatment ID"
NoT <- test.new[1:16, c(2, 7)]
orderNoT <- NoT$Treatment[order(NoT$Nooftreatment)]
temp <- sapply(orderNoT, FUN=function(x) {which(test.new$Treatment == x)} )
newID <- as.vector(temp)
NoTID <- rep(1:16, each=3)
test.new <- test.new[newID, ]
test.new$newTreatment <- NoTID
test.new$TreatmentID <- factor(test.new$newTreatment, labels=paste0("treat", 1:16))


yname <- all.index[1]
# yname <- all.y[1] 
B <- 200   # bootstrap number


## (1) The control effect and the single treatment effect
## control data
y.ct <- test.new[test.new$Nooftreatment==0, yname]
n.ct <- length(y.ct)   # sample size

## single treatment data
y.tr <- vector(mode="list", length=length(all.treat))
for (j in 1:length(all.treat)) {
  xname <- all.treat[j]
  y.tr[[j]] <- test.new[(test.new$Nooftreatment==1) & (test.new[xname] == "with"), yname]
}
names(y.tr) <- all.treat
n.tr <- sapply(y.tr, FUN=length)

# Bootstrap sampling for effect size of identity treatment (used in prediction)
ES.ct.all <- ES.tr.all <- as.numeric(0)
for (b in 1:B) {
  ES.ct <- mean(sample(y.ct, n.ct, replace=TRUE))
  ES.tr <- sapply(1:length(all.treat), FUN=function(x) { 
    mean(sample(y.tr[[x]], size=length(y.tr[[x]]), replace=TRUE)) })
  
  ES.ct.all <- c(ES.ct.all, ES.ct)                # control
  ES.tr.all <- rbind(ES.tr.all, ES.tr)            # single treatment
}


## (2) The combination effect (observed and predicted)
ES.joint.pred.all <- vector(mode="list", length=4)
ES.joint.obs.all <- vector(mode="list", length=4)

## Single treatment
# the observed effect
spID <- sample(1:B, size=B, replace=FALSE)
ES.add <- ES.tr.all - ES.ct.all
ES.identity.obs1 <- ES.add[-1, 1][spID[1:(B/4)]]
ES.identity.obs2 <- ES.add[-1, 2][spID[(1:(B/4)) + (B/4)]]
ES.identity.obs3 <- ES.add[-1, 3][spID[(1:(B/4)) + 2*(B/4)]]
ES.identity.obs4 <- ES.add[-1, 4][spID[(1:(B/4)) + 3*(B/4)]]
ES.joint.obs.all[[1]] <- c(ES.identity.obs1, ES.identity.obs2,
                           ES.identity.obs3, ES.identity.obs4)

# the predicted effect (use different samples)
spID <- sample(1:B, size=B, replace=FALSE)
ES.identity.pred1 <- ES.add[-1, 1][spID[1:(B/4)]]
ES.identity.pred2 <- ES.add[-1, 2][spID[(1:(B/4)) + (B/4)]]
ES.identity.pred3 <- ES.add[-1, 3][spID[(1:(B/4)) + 2*(B/4)]]
ES.identity.pred4 <- ES.add[-1, 4][spID[(1:(B/4)) + 3*(B/4)]]
pred.temp <- c(ES.identity.pred1, ES.identity.pred2,
               ES.identity.pred3, ES.identity.pred4)
ES.joint.pred.all[[1]] <- list(Additive=pred.temp, 
                               # Multiplicative=pred.temp,
                               Dominative=pred.temp)

## two to four treatments combined
for (i in 2:4) {
  ## Get the information about treatment ID
  temp <- test.new[test.new$Nooftreatment == i, ]
  combID <- apply(temp[temp$Block == 1, 3:6] == "with", 
                  MARGIN=1, FUN=which)
  treatID <- temp$Treatment[temp$Block == 1]
  ncb <- length(unique(temp$TreatmentID))       # number of combinations
  
  ## Effect size of the treatment combination
  B0 <- floor(B / ncb)
  ES.joint.obs.temp <- numeric(0)
  ES.joint.add.temp <- ES.joint.mult.temp <- ES.joint.dom.temp <- numeric(0)
  
  # sample from each treatment combination
  for (k in 1:ncb) {
    trID <- treatID[k]
    y.cb <- test.new[test.new$Treatment == trID, yname]
    n.cb <- length(y.cb)
    
    # the observed effect of combination k
    ES.comb <- numeric(0)
    for (b in 1:B0) {                    
      ES.ct <- ES.ct.all[-1][b + (k-1)*B0]                     # not from re-sampling
      # ES.ct <- mean(sample(y.ct, size=n.ct, replace=TRUE))   # from re-sampling
      ES.cb <- mean(sample(y.cb, size=n.cb, replace=TRUE))
      
      # compare to the control 
      ES.comb <- c(ES.comb, (ES.cb - ES.ct))
    }
    ES.joint.obs.temp <- c(ES.joint.obs.temp, ES.comb)
    
    # the predicted effect of combination k
    cbID <- combID[, k]
    spID <- sample(1:B, size=B0, replace=FALSE)
    ES.ct.temp <- ES.ct.all[-1][spID]
    ES.tr.temp <- ES.tr.all[-1, cbID][spID, ]
    
    # Additive
    ES.add.temp <- ES.tr.temp - ES.ct.temp
    ES.joint.temp <- apply(ES.add.temp, MARGIN=1, FUN=sum)
    ES.joint.add.temp <- c(ES.joint.add.temp, ES.joint.temp)
    
    # Multiplicative
    # ES.mult.temp <- ES.tr.temp / ES.ct.temp
    # ES.joint.temp <- sapply(1:length(spID), FUN=function(x) {
    #   (prod(ES.mult.temp[x, ]) - 1) * ES.ct.temp[x]
    # }) 
    # ES.joint.mult.temp <- c(ES.joint.mult.temp, ES.joint.temp)
    
    # Dominative  
    ES.dom.temp <- ES.tr.temp - ES.ct.temp
    ES.joint.temp <- apply(ES.dom.temp, MARGIN=1, FUN=max) 
    ES.joint.dom.temp <- c(ES.joint.dom.temp, ES.joint.temp)
    
  }   # end of the iteration of treatment combinations
  
  # save the result
  ES.joint.obs.all[[i]] <- ES.joint.obs.temp
  ES.joint.pred.all[[i]] <- list(Additive=ES.joint.add.temp, 
                                 # Multiplicative=ES.joint.mult.temp,
                                 Dominative=ES.joint.dom.temp)
}

# In the above analysis, the multiplicative effect was not computed. 
# This is because that the control effect is zero for some metrics,
# resulting in a zero multiplicative effect, which does not make sense.


## (3) Visualise the result
mean.obs <- sapply(ES.joint.obs.all, FUN=mean, na.rm=TRUE)
se.obs <- sqrt(sapply(ES.joint.obs.all, FUN=var, na.rm=TRUE))
upper.obs <- mean.obs + 1.96*se.obs
lower.obs <- mean.obs - 1.96*se.obs

ES.joint.additive <- sapply(1:4, FUN=function(x) { ES.joint.pred.all[[x]]$Additive })
ES.joint.dominative <- sapply(1:4, FUN=function(x) { ES.joint.pred.all[[x]]$Dominative })

ymin <- min(c(range(ES.joint.pred.all), range(ES.joint.obs.all)))
ymax <- max(c(range(ES.joint.pred.all), range(ES.joint.obs.all)))


df <- data.frame(NoT=rep(1:4, times=sapply(ES.joint.obs.all, FUN=length)),
                  Observed=unlist(ES.joint.obs.all), 
                  Additive=unlist(ES.joint.additive),
                  Dominative=unlist(ES.joint.dominative))
df$NoT <- as.factor(df$NoT)

plot1 <- ggplot(data=df) +
  geom_boxplot(mapping=aes(x=NoT, y=Observed), colour="steelblue4", fill="steelblue",
               size=0.3, outlier.shape=19, outlier.size=0.2) +
  theme(panel.background = element_rect(
    fill = "grey95",
    colour = "grey75",
    linewidth = 0.2
  ), 
  plot.margin = margin(4.5, 4.5, 2.5, 2.5, unit="pt"),
  text = element_text(size=4.8),
  axis.text = element_text(size=4.8),
  axis.title = element_text(size=5),
  plot.title = element_text(hjust = 1)) + 
  labs(x = 'Number of treatments',
       y = 'Observed effect',
       title = "(a)")


plot2 <- ggplot(data=df) +
  geom_boxplot(mapping=aes(x=NoT, y=Additive), colour="steelblue4", fill="steelblue",
               size=0.3, outlier.shape=19, outlier.size=0.2) +
  theme(panel.background = element_rect(
    fill = "grey95",
    colour = "grey75",
    linewidth = 0.2
  ), 
  plot.margin = margin(4.5, 4.5, 2.5, 2.5, unit="pt"),
  text = element_text(size=4.8),
  axis.text = element_text(size=4.8),
  axis.title = element_text(size=5),
  plot.title = element_text(hjust = 1)) + 
  labs(x = 'Number of treatments',
       y = 'Bootstrap additive effect',
       title = "(b)")


plot3 <- ggplot(data=df) +
  geom_boxplot(mapping=aes(x=NoT, y=Dominative), colour="steelblue4", fill="steelblue",
               size=0.3, outlier.shape=19, outlier.size=0.2) +
  theme(panel.background = element_rect(
    fill = "grey95",
    colour = "grey75",
    linewidth = 0.2
  ), 
  plot.margin = margin(4.5, 4.5, 2.5, 2.5, unit="pt"),
  text = element_text(size=4.8),
  axis.text = element_text(size=4.8),
  axis.title = element_text(size=5),
  plot.title = element_text(hjust = 1)) + 
  labs(x = 'Number of treatments',
       y = 'Bootstrap dominative effect',
       title = "(c)")


plot.all <- grid.arrange(plot1, plot2, plot3, nrow = 1)

ggsave(filename="Effect_NoT.png",  plot.all, device="png", 
       width=1500, height=500, units="px")


