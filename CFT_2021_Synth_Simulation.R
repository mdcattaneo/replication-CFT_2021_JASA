# Replication: simulation
# Cattaneo, Feng and Rocio (2020)
# Date: Feb 6, 2021

rm(list = ls())
#setwd("~/Dropbox/SC/simulation/")
library(tsDyn)
library(limSolve)
library(optiSolve)
library(Qtools)
library(foreach)
library(doParallel)
library(doRNG)
library(matrixStats)
source("SuppFuns_simul.R")

##########################################3
# pars
rep   <- 2000
T0    <- 100
J     <- 10
M     <- 400 
eq    <- 1
lb    <- 0
alpha <- 0.10
w0.u  <- c(0.3, 0.4, 0.3, rep(0, J-3))
rho.max <- 1

par <- read.csv("SimulModel.csv", header = T, colClasses=c(rep("numeric", 2), "logical", "numeric"))
num <- 1  # no. of dgp
model <- par$model[num]
err   <- par$err[num]
u.order <- par$order[num]

if (model==1) {
  # iid
  ar  <- matrix(0, J, J)
} else if (model==2) {
  # AR
  ar  <- diag(J)*0.5
} else {
  # cointegration
  ar  <- diag(J)
}

if (num %in% c(6, 9, 12)) rho.max <- 0.3 # set a maximum truncation only for model=3 with missp. error to avoid failures

#######################################
### Prepare design ####################
#######################################
set.seed(12345)
X.ls <- list()
for (i in 1:rep) {
   X.ls[[i]] <- dgp.x(T0, ar)
}

#######################################
# get true w0.c, will be used only when err=TRUE, otherwise =w0.u
if (err) {
   w0.c.ls <- sapply(X.ls, function(x) getw0(x[1:T0,], x[(T0+1),,drop=F], T0, J, model, w0.u, eq, lb))
}
#######################################

#######################################
### Conditional Prediction#############
#######################################
# pick a realization in the medium range
range.l <- rowQuantiles(sapply(X.ls, function(x) colQuantiles(x, probs = 0.1)), probs=.2)
range.r <- rowQuantiles(sapply(X.ls, function(x) colQuantiles(x, probs = 0.9)), probs=.8)
within <- sapply(X.ls, function(x) check(x, range.l, range.r))
select <- which(within==T)[1]
data.x <- X.ls[[select]]
y.co.0 <- data.x[1:T0,]
y.co.1 <- data.x[(T0+1),,drop=F]

# consider some evals
sd <- sd(y.co.0[,1])
y.co.eval <- t(sapply(seq(-1,1,0.5), function(j) y.co.1+c(j*sd, rep(0, J-1))))

if (err) {
  w0.c <- w0.c.ls[,select]
} else {
  w0.c <- w0.u
}

# simulate
cl <- makeCluster(32)
registerDoParallel(cl)
#writeLines(c(""), "log.txt")

output <- foreach (i = 1:rep, .options.RNG=12345, .packages=c('tsDyn','optiSolve','Qtools','matrixStats', 'limSolve'),
                   .combine=rbind) %dorng% {
                     #sink("log.txt", append=TRUE)
                     #cat(paste("Starting iteration", i, "\n"))
                     #sink()
                     #data   <- dgp.cond(y.co.0, y.co.eval, T0, J, model, err, w0.u, w0.c)
                     #output <- rbind(data$y.tr.0, data$y.tr.1)
                     output <- sc.sim.cond(i, y.co.0=y.co.0, y.co.1=y.co.eval, 
                                           T0, J, model, err, eq, lb, method.u="all", alpha, M,
                                           w0.u, w0.c, u.order, rho.max=rho.max) 
                     output   # (5*rep) by length(Kseq) matrix
                   }

stopCluster(cl)

###################
write.table(output, paste("rawoutput_cond_dgp", num, "txt", sep = "."), sep = ",", 
            row.names = F, col.names = F)

#########################################
### Unconditional coverage ##############
#########################################
if (err) {
  w0.c <- w0.c.ls
} else {
  w0.c <- matrix(w0.u, J, rep)
}

# simulate
cl <- makeCluster(32)
registerDoParallel(cl)
#writeLines(c(""), "log.txt")

output <- foreach (i = 1:rep, .options.RNG=12345, .packages=c('tsDyn','optiSolve','Qtools','matrixStats', 'limSolve'),
                   .combine=rbind) %dorng% {
                     #sink("log.txt", append=TRUE)
                     #cat(paste("Starting iteration", i, "\n"))
                     #sink()
                     #data.x <- X.ls[[i]]
                     #ty.co.0 <- data.x[1:T0,,drop=F]
                     #ty.co.1 <- data.x[(T0+1),,drop=F]
                     #tw0.c   <- w0.c[,i]
                     #data   <- dgp.cond(ty.co.0, ty.co.1, T0, J, model, err, w0.u, tw0.c)
                     #output <- rbind(data$y.tr.0, data$y.tr.1)
                     
                     output <- sc.sim.unc(i, x.list=X.ls, 
                                          T0, J, model, err, eq, lb, method.u="all", alpha, M,
                                          w0.u, w0.c.ls=w0.c, u.order=u.order, rho.max=rho.max) 
                     output   # (5*rep) by length(Kseq) matrix
                   }

stopCluster(cl)

###################
write.table(output, paste("rawoutput_unc_dgp", num, "txt", sep = "."), sep = ",", 
            row.names = F, col.names = F)


