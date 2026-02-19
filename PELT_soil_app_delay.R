## ------------------------------
## Soil moisture model results
## ------------------------------

## Soil water in different depths and travel time

library(nlmrt)
library(matrixcalc)
library(numDeriv)
library(zoo)



## soil moisture from 3 depths
# data.srer = readRDS("Application_data/SRER_Dep123.rds")
PELT.srer.3dep = readRDS("Application_result/PELT_SRER_3dep_output.rds")

# data.osbs = readRDS("Application_data/OSBS_Dep123.rds")
PELT.osbs.3dep = readRDS("Application_result/PELT_OSBS_3dep_output.rds")

# data.tall = readRDS("Application_data/TALL_Dep123.rds")
PELT.tall.3dep = readRDS("Application_result/PELT_TALL_3dep_output.rds")


## Matching the changepoints (e.g., TALL, location 4)
cpt.d1 <- PELT.tall.3dep$PELT.h4d1$cpt
cpt.d2 <- PELT.tall.3dep$PELT.h4d2$cpt
cpt.d3 <- PELT.tall.3dep$PELT.h4d3$cpt

# match the changepoints
cpt1 <- cpt.d1   # shallow
cpt2 <- cpt.d2   # deep
# first the exact match
exact.cpt <- cpt2[cpt2 %in% cpt1]
# then the relaxed match
relax.temp <- cpt2[!(cpt2 %in% cpt1)]
cpt1.temp <- cpt1[!(cpt1 %in% exact.cpt)]

relax.cpt <- NULL
remove.cpt <- NULL
for (j in 1:length(relax.temp)) {
  cpt.diff <- relax.temp[j] - cpt1.temp
  bin.select <- (cpt.diff > 0) & (cpt.diff < 24)   # condition
  cpt.match <- any(bin.select)
  if (cpt.match == TRUE) {
    relax.cpt <- c(relax.cpt, relax.temp[j])
    remove.cpt <- c(remove.cpt,  cpt1.temp[(which(bin.select)[1])])
    cpt1.temp <- cpt1.temp[-(which(bin.select)[1])]
  }
}

cpt1.match <- sort(c(exact.cpt, remove.cpt))
cpt2.match <- sort(c(exact.cpt, relax.cpt))

cpt1m <- cpt1.match
cpt2m <- cpt2.match

# add another depths
exact.cpt <- cpt.d3[cpt.d3 %in% cpt2m]
relax.temp <- cpt.d3[!(cpt.d3 %in% cpt2m)]
cpt2.temp <- cpt2m[!(cpt2m %in% exact.cpt)]

relax.cpt <- NULL
remove.cpt <- NULL
for (j in 1:length(relax.temp)) {
  cpt.diff <- relax.temp[j] - cpt2.temp
  bin.select <- (cpt.diff > 0) & (cpt.diff < 72)   # condition
  cpt.match <- any(bin.select)
  if (cpt.match == TRUE) {
    relax.cpt <- c(relax.cpt, relax.temp[j])
    remove.cpt <- c(remove.cpt,  cpt2.temp[(which(bin.select)[1])])
    cpt2.temp <- cpt2.temp[-(which(bin.select)[1])]
  } 
}

cpt3.mtemp <- sort(c(exact.cpt, relax.cpt))
cpt2.mtemp <- sort(c(exact.cpt, remove.cpt))
cpt3.match <- rep(NA, times=length(cpt2.match))
cpt3.match[cpt2.match %in% cpt2.mtemp] <- cpt3.mtemp

cpt.match <- cbind(cpt1.match, cpt2.match, cpt3.match)


## divide by the depths to get the travel time per cm
d12 <- (cpt.match[, 2] - cpt.match[, 1]) / 10
d23 <- (cpt.match[, 3] - cpt.match[, 2]) / 10
d13 <- (cpt.match[, 3] - cpt.match[, 1]) / 20

# divide by 2 if it is 30-min interval data
d12 = d12/2
d23 = d23/2
d13 = d13/2


## save the result
tall.delay.h4 <- cpt.match
tall.delay.h4 <- cbind(tall.delay.h4, d12=d12, d23=d23, d13=d13)


## Analyse the result
delay.list = readRDS("Application_result/PELT_3sites_delayed_peaks.rds")

osbs.delay = delay.list$osbs
srer.delay = delay.list$srer
tall.delay = delay.list$tall

# combine information by site
sitenames <- c("SRER", "OSBS", "TALL")

i <- 4
site.d12 <- list(srer=c(srer.delay$hloc2[, i], srer.delay$hloc4[, i]),
                 osbs=c(osbs.delay$hloc1[, i], osbs.delay$hloc2[, i], osbs.delay$hloc3[, i]),
                 tall=c(tall.delay$hloc2[, i], tall.delay$hloc3[, i], tall.delay$hloc4[, i]))

i <- 5
site.d23 <- list(srer=c(srer.delay$hloc2[, i], srer.delay$hloc4[, i]),
                 osbs=c(osbs.delay$hloc1[, i], osbs.delay$hloc2[, i], osbs.delay$hloc3[, i]),
                 tall=c(tall.delay$hloc2[, i], tall.delay$hloc3[, i], tall.delay$hloc4[, i]))

i <- 6
site.d13 <- list(srer=c(srer.delay$hloc2[, i], srer.delay$hloc4[, i]),
                 osbs=c(osbs.delay$hloc1[, i], osbs.delay$hloc2[, i], osbs.delay$hloc3[, i]),
                 tall=c(tall.delay$hloc2[, i], tall.delay$hloc3[, i], tall.delay$hloc4[, i]))


png(filename="Travel_time_hist_sites.png", width=1200, height=700, pointsize=16)
par(mfcol=c(2, 3), mar=c(4.5, 3.5, 3.5, 2.5), cex=1)
brk.list <- c(5, 5, 3)
for (i in 1:3) {
  hist(site.d12[[i]], breaks=brk.list[i], xlab="travel time d1-d2 (hours / cm)", 
       ylab="", main=sitenames[i], cex.main=1, xlim=c(0, 7))
  hist(site.d23[[i]], breaks=10, xlab="travel time d2-d3 (hours / cm)", 
       ylab="", main=sitenames[i], cex.main=1, xlim=c(0, 7))
}
dev.off()


