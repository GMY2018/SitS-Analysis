## ------------------------------------------
## Application: soil moisture signature
## ------------------------------------------

library(nlmrt)
library(matrixcalc)
library(numDeriv)
library(zoo)
library(lubridate)

source("PELT_soil_function_nlm.R")


## Output from the changepoint method
PELT.output = readRDS("Application_result/PELT_9sites_output_model.rds")
PELT.data = readRDS("Application_result/PELT_9sites_data.rds")


## ----------------------------------
## Get the coefficient matrix
## ----------------------------------
## Using SRER as an example

PELT = PELT.output$PELT.srer
sm = PELT.data$srer$sm
st = PELT.data$srer$time

lastpar <- PELT$lastchangecoef[PELT$cpt+1]
lastalpha0 <- sapply(1:length(lastpar), FUN=function(i) lastpar[[i]][1])
lastalpha1 <- sapply(1:length(lastpar), FUN=function(i) lastpar[[i]][2])
lastgamma <- sapply(1:length(lastpar), FUN=function(i) lastpar[[i]][3])
lastdrate <- sapply(1:length(lastpar), FUN=function(i) exp(-exp(lastpar[[i]][3])))
cmat <- data.frame(lastalpha0, lastalpha1, lastgamma)

# cmat.srer <- cmat


## remove the winter period (if necessary)
# cut.date <- which(st >= as.POSIXct("2022-01-01"))[1]
# cut.tp <- rev(PELT$cpt[PELT$cpt < cut.date])[1]
# cut.cpt <- which(PELT$cpt >= cut.date)[1] - 1
cut.cpt = length(lastalpha0)  # no cut off

last.alpha0 <- lastalpha0[1:cut.cpt]
last.alpha1 <- lastalpha1[1:cut.cpt]
last.gamma <- lastgamma[1:cut.cpt]
last.drate <- lastdrate[1:cut.cpt]

bad1 <- which(last.drate < 0.5)
bad2 <- which(last.alpha1 == 0.001)
bad3 <- which((last.alpha0 == 0) | (last.alpha0 > quantile(sm, probs=0.95)))
bad <- unique(c(bad1, bad2, bad3, cut.cpt:length(PELT$cpt)))

par(mfrow=c(1, 2))
hist(lastdrate, xlab="decay rate", ylab="", main="")
hist(lastdrate[-bad], xlab="decay rate", ylab="", main="")

hist(lastalpha0, xlab="alpha_0", ylab="", main="")
hist(lastalpha0[-bad], xlab="alpha_0", ylab="", main="")

hist(lastalpha1, xlab="alpha_1", ylab="", main="")
hist(lastalpha1[-bad], xlab="alpha_1", ylab="", main="")



## ---------------------------------
## Extract the standard errors
## ---------------------------------
## Using SRER as an example

data = PELT.data$srer$sm

sm.pelt = PELT.output$PELT.srer
cpt = c(0, sm.pelt$cpt)
lastchangecoef = sm.pelt$lastchangecoef
lastchangelike = sm.pelt$lastchangelike
thresh = 0.001
ini.asym = 0.01
ini.par = NULL


## Extracting the uncertainty
## Here we extracting the uncertainty by re-fitting the models and 
## record the standard deviations of the estimates

fmodel <- sm ~ asym + alpha0*exp(-exp(lgamma)*t)
fpar <- c(0.01, 0.1, -5)
names(fpar) <- c("asym", "alpha0", "lgamma")
exp.ss2 <- model2ssfun(fmodel, fpar, funname="exp.ss2")
exp.gr2 <- model2grfun(fmodel, fpar, funname="exp.gr2")
exp.hess2 <- function(par, data){
  Hessian <- array(0, dim=c(length(par), length(par), nrow(data)))
  eexp.gamma <- exp(-exp(par[3])*data$t)
  exp.gamma <- -exp(par[3])*data$t
  Hessian[2, 3, ] <- Hessian[3, 2, ] <- eexp.gamma * exp.gamma
  Hessian[3, 3, ] <- par[2] * eexp.gamma * exp.gamma * (exp.gamma + 1)
  return(-Hessian)
}


Jacobian.list <- vector(mode="list", length=length(cpt)-1)  # for the cross product
Deriv.list <- vector(mode="list", length=length(cpt)-1)
Resid.list <- Resid2.list <- vector(mode="list", length=length(cpt)-1)
sighat <- sighat2 <- rep(NA, length(cpt)-1)
SE.mat <- SE2.mat <- matrix(NA, nrow=length(cpt)-1, ncol=3)
coef.mat <- coef2.mat <- matrix(NA, nrow=length(cpt)-1, ncol=3)

for (i in 1:(length(cpt)-1)) {
  sID <- cpt[i] + 1
  eID <- cpt[i+1]
  y <- data[sID:eID]
  nt <- length(y)
  segsm <- data.frame(t=1:nt, sm=y)
  
  # the parameters
  coef.vec <- lastchangecoef[[eID+1]]
  
  # the asymptotic parameter
  if (is.null(ini.asym)) ini.asym <- round(0.8 * min(segsm$sm, na.rm=TRUE), 3)
  # the jump and drying parameters
  if (is.null(ini.par)) {
    ini.alpha0 <- round(max(segsm$sm, na.rm=TRUE), 3) - ini.asym
    ini.R <- try(ar(segsm$sm, order.max=1, na.action=na.pass)$ar, silent=TRUE)
    if (length(ini.R) == 1) {
      if ((ini.R < 1) & (ini.R > 0)) {
        ini.lgamma <- round(log(-log(ini.R)), 3)
      } else {
        ini.lgamma <- -5.3  # this is roughly log(-log(0.995))
      }
    } else if (length(ini.R) == 0) {
      ini.lgamma <- -5.3
    }
  } else {
    ini.alpha0 <- ini.par[1]
    ini.lgamma <- ini.par[2]
  }
  
  # re-fit the time series and get the Jacobian
  fit.nls <- try(
    nlfb(start=list(asym=ini.asym, alpha0=ini.alpha0, lgamma=ini.lgamma), 
         resfn=exp.res2, jacfn=exp.jac2, trace=FALSE, 
         lower=c(0, thresh, -20), upper=c(0.4, 0.5, 3),
         data=segsm),
    silent = TRUE
  )
  if (is.list(fit.nls)) {
    Resid.list[[i]] <- fit.nls$resid
    coef.mat[i, ] <- fit.nls$coefficients
    Jacobian.list[[i]] <- crossprod(fit.nls$jacobian)
    nls.sum <- summary(fit.nls)
    SE.mat[i, ] <- nls.sum$SEs
    ssnew <- exp.ss2(fit.nls$coefficients, sm=segsm$sm, t=segsm$t)
    sighat[i] <- ssnew / (nrow(segsm) - 3)
  }
  
  # or using nls() function
  fit.nls2 <- try(
    nls(sm ~ asym + alpha0*exp(-exp(lgamma)*t), data=segsm,
        start=list(asym=ini.asym, alpha0=ini.alpha0, lgamma=ini.lgamma),
        trace=FALSE, algorithm="port",
        lower=c(0, thresh, -20), upper=c(0.4, 0.5, 3)),
    silent = TRUE
  )
  if (is.list(fit.nls2)) {
    Resid2.list[[i]] <- residuals(fit.nls2)
    nls.sum2 <- summary(fit.nls2)
    coef2.mat[i, ] <- nls.sum2$coefficients[, 1]
    SE2.mat[i, ] <- nls.sum2$coefficients[, 2]
    ssnew <- exp.ss2(nls.sum2$coefficients[, 1], sm=segsm$sm, t=segsm$t)
    sighat2[i] <- ssnew / (nrow(segsm) - 3)
  }
}


## Sometimes the fit does not have a converged standard deviation
## Therefore, we combine the result from two different fitting functions to 
## maximise the chance that we have a converged result

diff.mat <- coef.mat - srer.coef[, 1:3]
dID1 <- which(apply(abs(diff.mat) > 0.05, MARGIN=1, FUN=any))
diff.mat <- coef2.mat - srer.coef[, 1:3]
dID2 <- which(apply(abs(diff.mat) > 0.05, MARGIN=1, FUN=any))

stdev.mat <- SE.mat
stdev.mat[dID1, ] <- NA

mID <- which(apply(is.na(stdev.mat), MARGIN=1, FUN=any))
coef.mat[mID, ]    # boundary solutions
cbind(coef.mat[mID, 3], SE2.mat[mID, 3])
rID <- mID[!(mID %in% dID2)]
stdev.mat[rID, ] <- SE2.mat[rID, ]


std.alpha0 <- stdev.mat[,1]
std.alpha1 <- stdev.mat[,2]
std.gamma <- stdev.mat[,3]
smat <- data.frame(std.alpha0, std.alpha1, std.gamma)

# smat.srer <- smat



## ----------------------------------
## comparing results from 9 sites
## ----------------------------------
cmat.list = readRDS("Application_result/PELT_9sites_output_par.rds")
smat.list = readRDS("Application_result/PELT_9sites_output_se.rds")

drate.list <- alpha0.list <- drate.list2 <- vector(mode="list", length=length(cmat.list))
for (i in 1:length(cmat.list)) {
  temp1 <- cmat.list[[i]]
  temp2 <- smat.list[[i]]
  alpha0.list[[i]] <- temp1[(!is.na(temp2[, 1]) & (temp2[, 1] <= 0.005) & (temp1[, 1] > 0)), 1]
  drate.list[[i]] <- exp(-exp(temp1[!is.na(temp2[, 1]) & (temp2[, 3] <= 0.5), 3]))
  drate.list2[[i]] <- temp1[!is.na(temp2[, 1]) & (temp2[, 3] <= 0.5), 3]
}


## Some plotting and summary tables
sapply(drate.list, FUN=length)
sapply(alpha0.list, FUN=length)

site.names <- c("SRER", "TALL", "OSBS", "UNDE", "CPER", "SCBI", 
                "ONAQ", "GUAN", "ORNL")
soil.types <- c("Entisol", "Ultisol", "Entisol", "Spodisol", "Mollisol", "Alfisol",
                "Aridisol", "Aridisol", "Ultisol")
soil.types <- as.factor(soil.types)
as.numeric(soil.types)

site.names[order(soil.types)]
newID <- order(soil.types)


png(file="9sites_good_alpha0.png", width=1000, height=550, pointsize=16)
par(mfrow=c(1, 1), mar=c(7, 4.5, 2.5, 2.5), cex=1)
boxplot(alpha0.list[newID], pch=20, names=paste(soil.types[newID], site.names[newID]),
        ylab="alpha_0", las=2)
dev.off()


png(file="9sites_good_decayrate.png", width=1000, height=550, pointsize=16)
par(mfrow=c(1, 1), mar=c(7, 4.5, 2.5, 2.5), cex=1)
boxplot(drate.list[newID], pch=20, names=paste(soil.types[newID], site.names[newID]),
        ylab="decay rate", las=2, ylim=c(0.8, 1))
dev.off()


png(file="9sites_good_params_soiltype.png", width=1000, height=700, pointsize=16)
par(mfrow=c(2, 1))
boxplot(alpha0.list[newID], pch=20, names=soil.types[newID], ylab="alpha_0")
boxplot(drate.list[newID], pch=20, names=soil.types[newID], ylab="decay rate", ylim=c(0.8, 1))
dev.off()


library(xtable)
atab <- t(sapply(alpha0.list, FUN=summary))
atab <- t(sapply(drate.list, FUN=summary))
row.names(atab) <- site.names
xtable(atab, caption="Summary statistics", digits=4)





