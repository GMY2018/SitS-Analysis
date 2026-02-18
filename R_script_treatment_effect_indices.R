## -------------------------------
## Script 5 for data analysis
## -------------------------------

# Multiple Comparison of the Effect of Different Treatment Combinations 
# on Ecosystem Service Multifunctionality


library(nlme)
library(lme4)
library(lmerTest)
library(stringr)
library(xtable)
library(multcomp)
library(multcompView)


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



## -----------------------------------------------
## Compare the effect of treatment combinations
## -----------------------------------------------
## Create a new variable called "Treatment ID"
NoT <- test.new[1:16, c(2, 7)]
orderNoT <- NoT$Treatment[order(NoT$Nooftreatment)]
temp <- sapply(orderNoT, FUN=function(x) {which(test.new$Treatment == x)} )
newID <- as.vector(temp)
NoTID <- rep(1:16, each=3)
test.new <- test.new[newID, ]
test.new$newTreatment <- NoTID
test.new$TreatmentID <- factor(test.new$newTreatment, labels=paste0("treat", 1:16))

# The treatment names (labels)
temp <- sapply(1:48, FUN=function(x) { names(test.new)[3:6][test.new[x, 3:6] == "with"] } )
temp[sapply(temp, FUN=length) == 0] <- "Control"
TreatName <- sapply(1:48, FUN=function(x) {str_flatten(temp[[x]], collapse="+")} )
test.new$TreatName <- factor(TreatName, levels=unique(TreatName))
test.new[, c(2, 7, 43, 44, 45)]


## (1) Fit models using treatment combinations as the main effect
## Iterate over all variables
all.y.List <- vector(mode="list", length=length(all.y))        # a list for storing the models
all.y.coef.List <- vector(mode="list", length=length(all.y))   # coefficient of treatments

for (i in 1:length(all.y)) {
  yname <- all.y[i] 
  
  # the model
  formula0 <- formula(paste(yname, "~ -1 + TreatmentID + (1|Block)"))
  lmer0 <- lmer(formula0, data=test.new, REML=TRUE)
  
  # influence from the fixed and random effect
  sum0 <- summary(lmer0)
  lmer0.coef <- sum0$coefficients[, c(1:2, 5)]
  
  # save the result
  all.y.List[[i]] <- lmer0
  colnames(lmer0.coef) <- c("coef", "std", "pvalue")
  all.y.coef.List[[i]] <- lmer0.coef
}


## Iterate over the multifunctional indices
all.index.List <- vector(mode="list", length=length(all.index))        # a list for storing the models
all.index.coef.List <- vector(mode="list", length=length(all.index))   # coefficient of treatments

for (i in 1:length(all.index)) {
  yname <- all.index[i] 
  
  # the baseline model
  formula0 <- formula(paste(yname, "~ -1 + TreatmentID + (1|Block)"))
  lmer0 <- lmer(formula0, data=test.new, REML=TRUE)
  
  # influence from the fixed and random effect
  sum0 <- summary(lmer0)
  lmer0.coef <- sum0$coefficients[, c(1:2, 5)]
  
  # save the result
  all.index.List[[i]] <- lmer0
  colnames(lmer0.coef) <- c("coef", "std", "pvalue")
  all.index.coef.List[[i]] <- lmer0.coef
}


## (2) Pairwise comparison plot (letter based display)
## plot for all 26 metrics
cld.list.ally <- list()
for (i in 1:26) {
  mod <- all.y.List[[i]]
  coef.mat <- all.y.coef.List[[i]]
  
  glht.obj <- glht(model=mod, linfct=mcp("TreatmentID" = "Tukey"))
  cld.obj <- cld(glht.obj)
  cld.list.ally[[i]] <- cld.obj
  
  # visualization
  df <- data.frame(x=cld.obj$x, y=cld.obj$y, treat=test.new$TreatName,
                   letter=rep(cld.obj$mcletters$Letters, each=3))
  
  plot1 <- ggplot(data=df) + 
    geom_boxplot(mapping=aes(x=treat, y=y),  colour="grey30", fill="grey60",
                 size=0.3, outlier.shape=19, outlier.size=0.2) +
    annotate("text", x=1:16, y=rep(1.1*max(cld.obj$y), 16), 
             label=as.vector(cld.obj$mcletters$Letters), 
             size=2, colour="grey30") +
    ylim(0.9*min(cld.obj$y), 1.1*max(cld.obj$y)) + 
    theme(panel.background = element_rect(
      fill = "grey95",
      colour = "grey75",
      linewidth = 0.2
    ), 
    plot.margin = margin(4.5, 4.5, 2.5, 2.5, unit="pt"),
    text = element_text(size=4.5),
    axis.text = element_text(size=4.5),
    axis.title = element_text(size=5),
    plot.title = element_text(hjust = 1),
    axis.text.x = element_text(angle = 50, hjust=1)) + 
    labs(x = 'Treatment combination',
         y = 'Effect from treatment')
  
  ggsave(filename=paste0("cld_", all.y[i], ".png"),  plot1, device="png", 
         width=1500, height=750, units="px")
  
}


## Plot for all indices
cld.list.index <- list()
for (i in 1:7) {
  mod <- all.index.List[[i]]
  coef.mat <- all.index.coef.List[[i]]

  glht.obj <- glht(model=mod, linfct=mcp("TreatmentID" = "Tukey"))
  cld.obj <- cld(glht.obj)
  cld.list.index[[i]] <- cld.obj
  
  df <- data.frame(x=cld.obj$x, y=cld.obj$y, treat=test.new$TreatName,
                   letter=rep(cld.obj$mcletters$Letters, each=3))
  
  plot2 <- ggplot(data=df) + 
    geom_boxplot(mapping=aes(x=treat, y=y),  colour="grey30", fill="grey60",
                 size=0.3, outlier.shape=19, outlier.size=0.2) +
    annotate("text", x=1:16, y=rep(1.1*max(cld.obj$y), 16), 
             label=as.vector(cld.obj$mcletters$Letters), 
             size=2, colour="grey30") +
    ylim(0.9*min(cld.obj$y), 1.1*max(cld.obj$y)) + 
    theme(panel.background = element_rect(
      fill = "grey95",
      colour = "grey75",
      linewidth = 0.2
    ), 
    plot.margin = margin(4.5, 4.5, 2.5, 2.5, unit="pt"),
    text = element_text(size=4.5),
    axis.text = element_text(size=4.5),
    axis.title = element_text(size=5),
    plot.title = element_text(hjust = 1),
    axis.text.x = element_text(angle = 50, hjust=1)) + 
    labs(x = 'Treatment combination',
         y = 'Effect from treatment')
  
  ggsave(filename=paste0("cld_", all.index[i], ".png"),  plot2, device="png", 
         width=1500, height=750, units="px")
}

# save(all.y.List, all.y.coef.List, all.index.List, all.index.coef.List,
#      cld.list.ally, cld.list.index, file="Colt_park_model3.RData")


