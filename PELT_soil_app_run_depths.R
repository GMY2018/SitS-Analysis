## -----------------------------------------
## Soil moisture, depth 1, 2, and 3
## -----------------------------------------

source("PELT_soil_function_nlm.R")


## Site SRER, 3 depths in locations 2, 3 and 4
## The hourly data
data.list = readRDS("Application_data/SRER_Dep123.rds")

data.l2.1h = data.list$hloc2
data.l3.1h = data.list$hloc3
data.l4.1h = data.list$hloc4

## Run the model using different penalty parameters
pen.list <- c(50, 100, 150, 200, 250)

## Fitting using the new function
PELT.result.h2d1 <- PELT.result.h2d2 <- PELT.result.h2d3 <- vector(mode="list", length=length(pen.list))
timing.h2d1 <- timing.h2d2 <- timing.h2d3 <- rep(0, length(pen.list))

PELT.result.h3d1 <- PELT.result.h3d2 <- PELT.result.h3d3 <- vector(mode="list", length=length(pen.list))
timing.h3d1 <- timing.h3d2 <- timing.h3d3 <- rep(0, length(pen.list))

PELT.result.h4d1 <- PELT.result.h4d2 <- PELT.result.h4d3 <- vector(mode="list", length=length(pen.list))
timing.h4d1 <- timing.h4d2 <- timing.h4d3 <- rep(0, length(pen.list))

# Model: y_t = asym + alpha0*exp(-exp(lgamma)*t) + e_t
for (i in 1:length(pen.list)) {
  # choose a penalty
  pen.i <- pen.list[i]

  # depth 1
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l2.1h$hloc_2_dep1, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h2d1[[i]] <- sm.pelt
  timing.h2d1[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l3.1h$hloc_3_dep1, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h3d1[[i]] <- sm.pelt
  timing.h3d1[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l4.1h$hloc_4_dep1, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h4d1[[i]] <- sm.pelt
  timing.h4d1[i] <- time2[1]
  

  # depth 2
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l2.1h$hloc_2_dep2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h2d2[[i]] <- sm.pelt
  timing.h2d2[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l3.1h$hloc_3_dep2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h3d2[[i]] <- sm.pelt
  timing.h3d2[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l4.1h$hloc_4_dep2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h4d2[[i]] <- sm.pelt
  timing.h4d2[i] <- time2[1]

  
  # depth 3
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l2.1h$hloc_2_dep3, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1

  PELT.result.h2d3[[i]] <- sm.pelt
  timing.h2d3[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l3.1h$hloc_3_dep3, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h3d3[[i]] <- sm.pelt
  timing.h3d3[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l4.1h$hloc_4_dep3, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h4d3[[i]] <- sm.pelt
  timing.h4d3[i] <- time2[1]

  print(i)

  save(PELT.result.h2d1, PELT.result.h2d2, PELT.result.h2d3,
       timing.h2d1, timing.h2d2, timing.h2d3, 
       PELT.result.h3d1, PELT.result.h3d2, PELT.result.h3d3,
       timing.h3d1, timing.h3d2, timing.h3d3,
       PELT.result.h4d1, PELT.result.h4d2, PELT.result.h4d3,
       timing.h4d1, timing.h4d2, timing.h4d3, 
       file="PELT_SRER_msl_3dep.RData")
}

# location 2, pen = 150, 150, 150
# location 4, pen = 200, 200, 200
## -----------------------------------------------------------------------------



## Site TALL, 3 depths in locations2, 3, 4
data.list = readRDS("Application_data/TALL_Dep123.rds")

data.l2.30m = data.list$hloc2
data.l3.30m = data.list$hloc3
data.l4.30m = data.list$hloc4

## Run the model using different penalty parameters
pen.list <- c(50, 100, 150, 200, 250)

## Fitting using the new function
PELT.result.h2d1 <- PELT.result.h2d2 <- PELT.result.h2d3 <- vector(mode="list", length=length(pen.list))
timing.h2d1 <- timing.h2d2 <- timing.h2d3 <- rep(0, length(pen.list))

PELT.result.h3d1 <- PELT.result.h3d2 <- PELT.result.h3d3 <- vector(mode="list", length=length(pen.list))
timing.h3d1 <- timing.h3d2 <- timing.h3d3 <- rep(0, length(pen.list))

PELT.result.h4d1 <- PELT.result.h4d2 <- PELT.result.h4d3 <- vector(mode="list", length=length(pen.list))
timing.h4d1 <- timing.h4d2 <- timing.h4d3 <- rep(0, length(pen.list))


# Model: y_t = asym + alpha0*exp(-exp(lgamma)*t) + e_t
for (i in 1:length(pen.list)) {
  # choose a penalty
  pen.i <- pen.list[i]

  # depth 1
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l2.30m$hloc_2_dep1, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h2d1[[i]] <- sm.pelt
  timing.h2d1[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l3.30m$hloc_3_dep1, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1

  PELT.result.h3d1[[i]] <- sm.pelt
  timing.h3d1[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l4.30m$hloc_4_dep1, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h4d1[[i]] <- sm.pelt
  timing.h4d1[i] <- time2[1]

  
  # depth 2
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l2.30m$hloc_2_dep2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h2d2[[i]] <- sm.pelt
  timing.h2d2[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l3.30m$hloc_3_dep2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h3d2[[i]] <- sm.pelt
  timing.h3d2[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l4.30m$hloc_4_dep2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h4d2[[i]] <- sm.pelt
  timing.h4d2[i] <- time2[1]
  
  
  # depth 3
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l2.30m$hloc_2_dep3, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1

  PELT.result.h2d3[[i]] <- sm.pelt
  timing.h2d3[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l3.30m$hloc_3_dep3, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h3d3[[i]] <- sm.pelt
  timing.h3d3[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l4.30m$hloc_4_dep3, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h4d3[[i]] <- sm.pelt
  timing.h4d3[i] <- time2[1]
  

  print(i)

  save(PELT.result.h2d1, PELT.result.h2d2, PELT.result.h2d3,
       timing.h2d1, timing.h2d2, timing.h2d3, 
       PELT.result.h3d1, PELT.result.h3d2, PELT.result.h3d3,
       timing.h3d1, timing.h3d2, timing.h3d3,
       PELT.result.h4d1, PELT.result.h4d2, PELT.result.h4d3,
       timing.h4d1, timing.h4d2, timing.h4d3,
       file="PELT_TALL_msl_3dep.RData")
}

# location 2, pen = 250, 150, 150
# location 3, pen = 200, 150, 250
# location 4, pen = 200, 150, 200
## ------------------------------------------------------------------



## Site OSBS, 3 depths in locations 1, 2, 3
data.list = readRDS("Application_data/OSBS_Dep123.rds")

data.l1.1h = data.list$hloc1
data.l2.1h = data.list$hloc2
data.l3.1h = data.list$hloc3

## Run the model using different penalty parameters
pen.list <- c(50, 100, 150, 200, 250)

## Fitting using the new function
PELT.result.h1d1 <- PELT.result.h1d2 <- PELT.result.h1d3 <- vector(mode="list", length=length(pen.list))
timing.h1d1 <- timing.h1d2 <- timing.h1d3 <- rep(0, length(pen.list))

PELT.result.h2d1 <- PELT.result.h2d2 <- PELT.result.h2d3 <- vector(mode="list", length=length(pen.list))
timing.h2d1 <- timing.h2d2 <- timing.h2d3 <- rep(0, length(pen.list))

PELT.result.h3d1 <- PELT.result.h3d2 <- PELT.result.h3d3 <- vector(mode="list", length=length(pen.list))
timing.h3d1 <- timing.h3d2 <- timing.h3d3 <- rep(0, length(pen.list))


# Model: y_t = asym + alpha0*exp(-exp(lgamma)*t) + e_t
for (i in 1:length(pen.list)) {
  # choose a penalty
  pen.i <- pen.list[i]
  
  # depth 1
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l1.1h$hloc_1_dep1, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h1d1[[i]] <- sm.pelt
  timing.h1d1[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l2.1h$hloc_2_dep1, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h2d1[[i]] <- sm.pelt
  timing.h2d1[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l3.1h$hloc_3_dep1, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h3d1[[i]] <- sm.pelt
  timing.h3d1[i] <- time2[1]
  
  
  # depth 2
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l1.1h$hloc_1_dep2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h1d2[[i]] <- sm.pelt
  timing.h1d2[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l2.1h$hloc_2_dep2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h2d2[[i]] <- sm.pelt
  timing.h2d2[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l3.1h$hloc_3_dep2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h3d2[[i]] <- sm.pelt
  timing.h3d2[i] <- time2[1]
  
  
  # depth 3
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l1.1h$hloc_1_dep3, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h1d3[[i]] <- sm.pelt
  timing.h1d3[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l2.1h$hloc_2_dep3, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h2d3[[i]] <- sm.pelt
  timing.h2d3[i] <- time2[1]
  
  
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=data.l3.1h$hloc_3_dep3, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  
  PELT.result.h3d3[[i]] <- sm.pelt
  timing.h3d3[i] <- time2[1]
  
  
  print(i)
  
  save(PELT.result.h1d1, PELT.result.h1d2, PELT.result.h1d3,
       timing.h1d1, timing.h1d2, timing.h1d3, 
       PELT.result.h2d1, PELT.result.h2d2, PELT.result.h2d3,
       timing.h2d1, timing.h2d2, timing.h2d3,
       PELT.result.h3d1, PELT.result.h3d2, PELT.result.h3d3,
       timing.h3d1, timing.h3d2, timing.h3d3,
       file="PELT_OSBS_msl_3dep.RData")
}

# location 1, pen = 150, 150, 150
# location 2, pen = 150, 150, 250
# location 3, pen = 200, 150, 150
## ------------------------------------------------------------------

