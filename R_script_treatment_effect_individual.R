## Estimation of Treatment Effects on Individual Ecosystem Service Indicators 
## Using Mixed Effects Model and Visualization

library(nlme)
library(lme4)
library(lmerTest)
library(stringr)
library(xtable)


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
## (I) Fit mixed effect model for 26 metrics
## -----------------------------------------------
## Define the scope of the linear effect and interactions
main.var <- c("Seed", "Clover", "NPK", "FYM")
twoway.var <- c("Seed*Clover", "Seed*NPK", "Seed*FYM", "Clover*NPK", 
                "Clover*FYM", "NPK*FYM")

## Determine the significance level used in variable selection
sig.thresh <- 0.1 

## List to save the result
lmer.List <- vector(mode="list", length=length(all.y))   # list of models
beta.List <- vector(mode="list", length=length(all.y))   # list of coefficients, etc.
mixedTF <- rep(1, length(all.y))    # whether the model has random effect or not

for (i in 1:length(all.y)) {
  ## The baseline model
  yname <- all.y[i] 
  formula0 <- formula(paste(yname, "~", 
                            str_flatten(main.var, collapse=" + "), 
                            "+ (1|Block)"))
  
  lmer0 <- lmer(formula0, data=test.new, REML=TRUE)
  # backward selection (dropping variables one by one)
  drop.lmer0 <- drop1(lmer0, scope=main.var, data=test.new)
  
  # keep the main effect that has a significant level of "sig.thresh"
  keep.main <- main.var[-which(drop.lmer0$`Pr(>F)` > sig.thresh)]
  if (length(keep.main) == 0) {
    formula.update <- formula(paste(yname, "~", "(1|Block)"))
  } else {
    formula.update <- formula(paste(yname, "~", 
                                    str_flatten(keep.main, collapse=" + "), 
                                    "+ (1|Block)"))
  }
  
  # the selected baseline model is 
  lmer00 <- lmer(formula.update, data=test.new, REML=TRUE)
  if (isSingular(lmer00)) {
    # drop the random effect if it has no effect on the model
    if (length(keep.main) == 0) {
      formula.update <- formula(paste(yname, "~ 1"))
    } else {
      formula.update <- formula(paste(yname, "~", 
                                      str_flatten(keep.main, collapse=" + ")))
    }
    lmer00 <- lm(formula.update, data=test.new)
    mixedTF[i] <- 0 
    
    # adding interaction to the simple linear model
    lmer1 <- lmer00 
    form.null <- formula.update
    form.full <- formula(paste(str_flatten(as.character(form.null)[c(2, 1, 3)], collapse = " "),
                               "+", str_flatten(twoway.var, collapse = " + ")))
    # use step wise selection to determine which interactions to include
    lmer1 <- step(lmer1, direction="both", 
                  scope=list(lower=form.null, upper=form.full), 
                  data=test)
    lmer.List[[i]] <- lmer1
    
    # extract the coefficients, standard deviation and P-values
    lmer1.summary <- summary(lmer1)
    lmer1.coef <- lmer1.summary$coefficients[, 1]
    lmer1.coef.sd <- lmer1.summary$coefficients[, 2]
    ci.upper <- lmer1.coef + 1.96*lmer1.coef.sd
    ci.lower <- lmer1.coef - 1.96*lmer1.coef.sd
    pvalue <- lmer1.summary$coefficients[, 4]
    
    beta.df <- data.frame(beta_est = lmer1.coef,  CI_lower = ci.lower, 
                          CI_upper = ci.upper, P_value=pvalue)
    
    beta.List[[i]] <- beta.df
    
  } else {
    # adding the interactions to the mixed effect model
    # there is no equivalence to step() for lmer(), so we do it manually
    lmer1 <- lmer00 
    keep.twoway <- twoway.var
    while (!is.null(keep.twoway)) {
      add.twoway <- add1(lmer1, data=test.new, scope=keep.twoway) 
      addID <- which.min(add.twoway$AIC) - 1   # don't count the intercept, so -1
      if (addID == 0) {
        keep.twoway <- NULL
        # make the scope empty and exit the loop
      } else {
        # update the model
        formula.update <- formula(paste(str_flatten(as.character(formula.update)[c(2, 1, 3)], collapse = " "), 
                                        "+", keep.twoway[addID]))
        lmer1 <- lmer(formula.update, data=test.new, REML=TRUE)
        
        # update the scope
        keep.twoway <- keep.twoway[-addID]
      }
    }
    lmer.List[[i]] <- lmer1
    
    # Effect size calculation 
    lmer1.summary <- summary(lmer1)
    lmer1.coef <- lmer1.summary$coefficients[, 1]
    
    lmer1.randvar <- as.numeric(VarCorr(lmer1))
    lmer1.resdvar <- lmer1.summary$sigma^2
    all.var <- sum(lmer1.randvar, lmer1.resdvar)
    effect.size <- lmer1.coef / sqrt(all.var)
    
    # Confidence interval and P-values of the fixed effect coefficients
    lmer1.coef.sd <- lmer1.summary$coefficients[, 2]
    ci.upper <- lmer1.coef + 1.96*lmer1.coef.sd
    ci.lower <- lmer1.coef - 1.96*lmer1.coef.sd
    pvalue <- lmer1.summary$coefficients[, 5]
    
    beta.df <- data.frame(beta_est = lmer1.coef,  CI_lower = ci.lower, CI_upper = ci.upper, 
                          P_value=pvalue, effect_size=effect.size)
    
    beta.List[[i]] <- beta.df
  }
  
}

names(lmer.List) <- all.y
names(beta.List) <- all.y

# save(lmer.List, beta.List, mixedTF, file="Colt_Park_model1.RData")



## ----------------------------------------------------------
## (II) Coefficient and P-value matrices for visualisation
## ----------------------------------------------------------
all.terms <- c("Seedwith", "Cloverwith", "NPKwith", "FYMwith", 
               "Seedwith:Cloverwith", "Seedwith:NPKwith", "Seedwith:FYMwith", 
               "Cloverwith:NPKwith", "Cloverwith:FYMwith", "NPKwith:FYMwith")
terms.name <- str_split(all.terms, pattern=":")

coef.mat <- pvalue.mat <- matrix(0, nrow=length(all.terms), ncol=length(all.y))
for (i in 1:length(all.y)) {
  beta.mat <- beta.List[[i]]
  beta.name <- rownames(beta.mat)[-1]
  beta.names <- str_split(beta.name, pattern=":")
  nt <- sapply(beta.names, FUN=length)
  mID <- rep(0, length(beta.name))
  for (j in which(nt==1)) {
    mID[j] <- match(beta.name[j], all.terms)
  }
  for (j in which(nt==2)) {
    mID[j] <- which(sapply(1:length(terms.name), FUN=function(k) {
      all(beta.names[[j]] %in% terms.name[[k]]) } ))
  }
  
  coef.mat[mID, i] <- beta.mat[-1, 1]
  pvalue.mat[mID, i] <- beta.mat[-1, 4]
}

rownames(coef.mat) <- all.terms
colnames(coef.mat) <- all.y
rownames(pvalue.mat) <- all.terms
colnames(pvalue.mat) <- all.y

# save(coef.mat, pvalue.mat, file="Modelling_result.RData")
