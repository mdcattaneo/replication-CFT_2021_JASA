# Replication: simulation
# Cattaneo, Feng and Rocio (2020)
# Notes: this file is only used to reproduce tables in the paper
# Date: Feb 6, 2021

library(Hmisc)
setwd("~/Dropbox/SC/simulation/")
source("SuppFuns.R")


####################################
table <- NULL
rej.ind <- 1:7

# Four tables: num=1:3, 4:6, 7:9, 10:12
for (num in 1:3) {

### conditional ####################
output <- as.matrix(read.table(paste("rawoutput/rawoutput_cond_dgp", num, "txt", sep="."), sep = ","))

summary.c <- matrix(colMeans(output), 5)
summary.c[,rej.ind] <- sapply(rej.ind, function(i) 1-summary.c[,i])
colnames(summary.c) <- c("rej.w", "rej.1", "rej.2", "rej.3", "rej.4", "rej.per", "rej.conf", 
                       "len.w", "len.1", "len.2", "len.3", "len.4", "len.per", "len.conf")

summary.c <- format(round(summary.c[,c("rej.w", "len.w", "rej.1", "len.1", "rej.4", "len.4",
                      "rej.2", "len.2", "rej.3", "len.3", "rej.per", "len.per", 
                      "rej.conf", "len.conf")], 3), 3)
#summary.c <- cbind(1:5, summary.c)



###########################################
### unconditional #########################
output <- as.matrix(read.table(paste("rawoutput/rawoutput_unc_dgp", num, "txt", sep="."), sep = ","))

summary.u <- t(colMeans(output))
summary.u[,rej.ind] <- sapply(rej.ind, function(i) 1-summary.u[,i])

colnames(summary.u) <- c("rej.w", "rej.1", "rej.2", "rej.3", "rej.4", "rej.per", "rej.conf", 
                       "len.w", "len.1", "len.2", "len.3", "len.4", "len.per", "len.conf")
summary.u <- format(round(summary.u[,c("rej.w", "len.w", "rej.1", "len.1", "rej.4", "len.4",
                                "rej.2", "len.2", "rej.3", "len.3", "rej.per", "len.per",
                                "rej.conf", "len.conf")], 3), 3)
#summary.u <- cbind("", summary.u)
table <- rbind(table, summary.c, summary.u)    
}

n.cgroup <- rep(2, 6)
colheads <- rep(c("CP", "AL"), 6)
cgroup   <- c("M1", "M1-S", "M2", "M3", "PERM", "CONF")

n.rgroup <- rep(6, 3)
rgroup   <- c("$\\rho=0$", "$\\rho=0.5$", "$\\rho=1$")
rowname  <- rep(c("Cond. 1", "2", "3", "4", "5", "Uncond."), 3)
latex(table[,-c(1,2)], file=paste("Table_Main", num/3, ".txt", sep = ""), rowlabel.just="r",  
      append=FALSE, table.env=FALSE, center="none", title="",
      n.cgroup=n.cgroup, cgroup=cgroup, colheads=colheads,
      n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname
)
