# Replication file: empirical illustration
# Cattaneo, Feng and Titiunik (2021)
# Notes: replication needs to load supporting functions first
# Date: Aug 7, 2021

rm(list=ls())
library(ggplot2)
library(optiSolve)
library(Qtools)
library(latex2exp)
library(Hmisc)
#setwd("~/Dropbox/SC/Empirical_Illustration/")
source("CFT_2021_JASA_empapp-fun.R")
set.seed(1234)

##############################
###### Germany ###############
data <- read.csv("Germany.csv")

### 1st differenced data ####
t <- data$year[data$index==1][-1]
y <- data$gdp[data$index==7]
x <- matrix(data$gdp[data$index!=7], 44)
y <- diff(log(data$gdp[data$index==7]))*100
x <- apply(log(x), 2, diff)*100

A <- y[t<=1990]
B <- x[t<=1990,]
C <- NULL
L <- 7
x.T <- x[(t>=1990+1&t<=1990+L),]
eq <- 1
lb <- 0
method.u <- "all"
alpha <- 0.1
M <- 100
result <- sc.pi(A=A, B=B, C=C, x.T=x.T, eq=eq, lb=lb, method.u=method.u, alpha=alpha, M=M, model=1, u.order=1, constant=F)

sc <- c(result$fit.pr, result$fit.po)
na <- rep(NA, nrow(B))
sc.l.0  <- c(na, result$sc.l.0); sc.r.0  <- c(na, result$sc.r.0)
sc.l.1  <- c(na, result$sc.l.1); sc.r.1  <- c(na, result$sc.r.1)
sc.l.2  <- c(na, result$sc.l.2); sc.r.2  <- c(na, result$sc.r.2)
sc.l.3  <- c(na, result$sc.l.3); sc.r.3  <- c(na, result$sc.r.3)

# Plot
dat    <- data.frame(t=t[t<=1990+L], y=y[t<=1990+L], sname="Germany")
dat.sc <- data.frame(t=t[t<=1990+L], sc=sc, 
                     lb0=sc.l.0, ub0=sc.r.0, 
                     lb1=sc.l.1, ub1=sc.r.1,
                     lb2=sc.l.2, ub2=sc.r.2,
                     lb3=sc.l.3, ub3=sc.r.3,
                     sname="SC unit")
pdf.width <- 6
pdf.height <- 4.5

plot <- ggplot() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="year", y="growth in GDP per capita (%)") +
  theme(legend.position = c(0.15,0.15), 
        legend.background = element_rect(fill = "transparent"))
plot <- plot + geom_line(data=dat, aes(x=t, y=y, colour=sname), linetype=1) +
  geom_point(data=dat, aes(x=t, y=y, colour=sname), shape=1) +
  geom_line(data=dat.sc, aes(x=t, y=sc, colour=sname), linetype=5) +
  geom_point(data=dat.sc, aes(x=t, y=sc, colour=sname), col="blue", shape=19) +
  geom_vline(xintercept = 1990, linetype="dotted") +
  geom_text(aes(x=1990, label="\nreunification", y=10), angle=90, size=4) +
  scale_x_continuous(breaks=c(seq(1960, 1990, 10), 1997))

plot0 <- plot + scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                     linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/Germany_SC.pdf", width=pdf.width, height=pdf.height)

plot.w <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                               col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/Germany_SCwithPI_w.pdf", width=pdf.width, height=pdf.height)

plot1 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb1, ymax=ub1), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/Germany_SCwithPI_1.pdf", width=pdf.width, height=pdf.height)

plot2 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb2, ymax=ub2), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/Germany_SCwithPI_2.pdf", width=pdf.width, height=pdf.height)

plot3 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb3, ymax=ub3), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/Germany_SCwithPI_3.pdf", width=pdf.width, height=pdf.height)


#### Sensitivity analysis #######
# Gaussian bound:
multi <- c(0.25, 0.5, 1, 1.5, 2)
time <- c(1993)

for (l in 1:length(time)) {
  ssc.l.1 <- ssc.r.1 <- c()
  u.mean <- result$u.1[time[l]-1990]
  sig <- sqrt(result$u.2[time[l]-1990])
  sig.seq <- multi*sig
  
  for (s in 1:length(sig.seq)) {
    eps  <- sqrt(-log(alpha/4)*2*(sig.seq[s]^2))
    ssc.l.1[s] <- sc.l.0[t==time[l]] + u.mean - eps
    ssc.r.1[s] <- sc.r.0[t==time[l]] + u.mean + eps
  }
  
  sen.dat <- data.frame(t=c(1:5), lb1=ssc.l.1, ub1=ssc.r.1, 
                        lb=rep(sc.l.0[t==time[l]], 5),
                        ub=rep(sc.r.0[t==time[l]], 5), 
                        lab=as.factor(multi))
  plot <- ggplot() + theme_bw() + 
                     theme(panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank()) +
                           labs(x="sd. of e", y="growth in GDP per capita (%)")
  plot <- plot + geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb1, ymax=ub1), 
                               col="maroon", width=0.2, linetype=5) +
                 geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb, ymax=ub), 
                               col="blue", width=0.2, linetype=1) +
                 geom_hline(yintercept = y[t==time[l]], linetype=1, size=0.3, alpha=0.8) +
                 annotate("text", x=5.4, y=y[t==time[l]]-.2,label="Y(1)", size=3.5) +
                 scale_x_discrete(labels=c(parse(text=TeX("$0.25\\hat{\\sigma}$")), 
                                           parse(text=TeX("$0.5\\hat{\\sigma}$")),
                                           parse(text=TeX("$\\hat{\\sigma}$")), 
                                           parse(text=TeX("$1.5\\hat{\\sigma}$")),                                
                                           parse(text=TeX("$2\\hat{\\sigma}$"))))
  ggsave(paste("figures/Germany_sens_", l, ".pdf", sep=""), width=pdf.width, height=pdf.height)
}


#########################
### Level ###############
data <- read.csv("Germany.csv")
t <- data$year[data$index==1]
y <- data$gdp[data$index==7]/1000
x <- matrix(data$gdp[data$index!=7], 44)/1000

A <- y[t<=1990]
B <- x[t<=1990,]
C <- NULL
L <- 7
x.T <- x[(t>=1990+1&t<=1990+L),]
eq <- 1
lb <- 0
method.u <- "all"
alpha <- 0.1
M <- 100
result <- sc.pi(A=A, B=B, C=C, x.T=x.T, eq=eq, lb=lb, method.u=method.u, alpha=alpha, M=M, model=3, u.order=1, constant=F)

sc <- c(result$fit.pr, result$fit.po)
na <- rep(NA, nrow(B))
sc.l.0  <- c(na, result$sc.l.0); sc.r.0  <- c(na, result$sc.r.0)
sc.l.1  <- c(na, result$sc.l.1); sc.r.1  <- c(na, result$sc.r.1)
sc.l.2  <- c(na, result$sc.l.2); sc.r.2  <- c(na, result$sc.r.2)
sc.l.3  <- c(na, result$sc.l.3); sc.r.3  <- c(na, result$sc.r.3)

# Plot
dat    <- data.frame(t=t[t<=1990+L], y=y[t<=1990+L], sname="Germany")
dat.sc <- data.frame(t=t[t<=1990+L], sc=sc, 
                     lb0=sc.l.0, ub0=sc.r.0, 
                     lb1=sc.l.1, ub1=sc.r.1,
                     lb2=sc.l.2, ub2=sc.r.2,
                     lb3=sc.l.3, ub3=sc.r.3,
                     sname="SC unit")
pdf.width <- 6
pdf.height <- 4.5

plot <- ggplot() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="year", y="GDP per capita (thousand US dollars)") +
  theme(legend.position = c(0.1, 0.95), 
        legend.background = element_rect(fill = "transparent"))
plot <- plot + geom_line(data=dat, aes(x=t, y=y, colour=sname), linetype=1) +
  geom_point(data=dat, aes(x=t, y=y, colour=sname), shape=1) +
  geom_line(data=dat.sc, aes(x=t, y=sc, colour=sname), linetype=5) +
  geom_point(data=dat.sc, aes(x=t, y=sc, colour=sname), col="blue", shape=19) +
  geom_vline(xintercept = 1990, linetype="dotted") +
  geom_text(aes(x=1990, label="\nreunification", y=10), angle=90, size=4) +
  scale_x_continuous(breaks=c(seq(1960, 1990, 10), 1997))

plot0 <- plot + scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                     linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/LevelGermany_SC.pdf", width=pdf.width, height=pdf.height)

plot.w <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                               col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/LevelGermany_SCwithPI_w.pdf", width=pdf.width, height=pdf.height)

plot1 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb1, ymax=ub1), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/LevelGermany_SCwithPI_1.pdf", width=pdf.width, height=pdf.height)

plot2 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb2, ymax=ub2), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/LevelGermany_SCwithPI_2.pdf", width=pdf.width, height=pdf.height)

plot3 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb3, ymax=ub3), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/LevelGermany_SCwithPI_3.pdf", width=pdf.width, height=pdf.height)


#### Sensitivity analysis #######
# Gaussian bound:
multi <- c(0.25, 0.5, 1, 1.5, 2)
time <- c(1993)

for (l in 1:length(time)) {
  ssc.l.1 <- ssc.r.1 <- c()
  u.mean <- result$u.1[time[l]-1990]
  sig <- sqrt(result$u.2[time[l]-1990])
  sig.seq <- multi*sig
  
  for (s in 1:length(sig.seq)) {
    eps  <- sqrt(-log(alpha/4)*2*(sig.seq[s]^2))
    ssc.l.1[s] <- sc.l.0[t==time[l]] + u.mean - eps
    ssc.r.1[s] <- sc.r.0[t==time[l]] + u.mean + eps
  }
  
  sen.dat <- data.frame(t=c(1:5), lb1=ssc.l.1, ub1=ssc.r.1, 
                        lb=rep(sc.l.0[t==time[l]], 5),
                        ub=rep(sc.r.0[t==time[l]], 5), 
                        lab=as.factor(multi))
  plot <- ggplot() + theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x="sd. of e", y="GDP per capita (thousand US dollars)")
  plot <- plot + geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb1, ymax=ub1), 
                               col="maroon", width=0.2, linetype=5) +
    geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb, ymax=ub), 
                  col="blue", width=0.2, linetype=1) +
    geom_hline(yintercept = y[t==time[l]], linetype=1, size=0.3, alpha=0.8) +
    annotate("text", x=5.4, y=y[t==time[l]]-.05,label="Y(1)", size=3.5) +
    scale_x_discrete(labels=c(parse(text=TeX("$0.25\\hat{\\sigma}$")), 
                              parse(text=TeX("$0.5\\hat{\\sigma}$")),
                              parse(text=TeX("$\\hat{\\sigma}$")), 
                              parse(text=TeX("$1.5\\hat{\\sigma}$")),                                
                              parse(text=TeX("$2\\hat{\\sigma}$"))))
  ggsave(paste("figures/LevelGermany_sens_", l, ".pdf", sep=""), width=pdf.width, height=pdf.height)
}



############################################################################
######### California Tobacco Control #######################################
############################################################################
data <- read.csv("California.csv")

t <- data$year[data$state=="California"][-1]
y <- data$cigsale[data$state=="California"]
x <- data[data$state!="California",1:3]
x <- as.matrix(reshape(x, timevar = "state", idvar="year", direction = "wide"))[,-1]
y <- diff(log(y))*100
x <- apply(log(x), 2, diff)*100

# cutoff: 1988
A <- y[t<=1988]
B <- x[t<=1988,]
C <- NULL
L <- 7
x.T <- x[(t>=1988+1&t<=1988+L),]
eq <- 1
lb <- 0
method.u <- "all"
alpha <- 0.1
M <- 100
result <- sc.pi(A=A, B=B, C=C, x.T=x.T, eq=eq, lb=lb, method.u=method.u, alpha=alpha, M=M, model=1, u.order=1, constant=F)

sc <- c(result$fit.pr, result$fit.po)
na <- rep(NA, nrow(B))
sc.l.0  <- c(na, result$sc.l.0); sc.r.0  <- c(na, result$sc.r.0)
sc.l.1  <- c(na, result$sc.l.1); sc.r.1  <- c(na, result$sc.r.1)
sc.l.2  <- c(na, result$sc.l.2); sc.r.2  <- c(na, result$sc.r.2)
sc.l.3  <- c(na, result$sc.l.3); sc.r.3  <- c(na, result$sc.r.3)


# Plot
dat    <- data.frame(t=t[t<=1988+L], y=y[t<=1988+L], sname="California")
dat.sc <- data.frame(t=t[t<=1988+L], sc=sc, 
                     lb0=sc.l.0, ub0=sc.r.0, 
                     lb1=sc.l.1, ub1=sc.r.1,
                     lb2=sc.l.2, ub2=sc.r.2,
                     lb3=sc.l.3, ub3=sc.r.3,
                     sname="SC unit")
pdf.width <- 6
pdf.height <- 4.5

plot <- ggplot() + theme_bw() + 
                   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                   labs(x="year", y="growth in cigarette sales per capita (%)") +
                   theme(legend.position = c(0.15,0.15), 
                   legend.background = element_rect(fill = "transparent"))
plot <- plot + geom_line(data=dat, aes(x=t, y=y, colour=sname), linetype=1) +
               geom_point(data=dat, aes(x=t, y=y, colour=sname), shape=1) +
               geom_line(data=dat.sc, aes(x=t, y=sc, colour=sname), linetype=5) +
               geom_point(data=dat.sc, aes(x=t, y=sc, colour=sname), col="blue", shape=19) +
               geom_vline(xintercept=1988, linetype="dotted") +
               geom_text(aes(x=1988, label="\nProposition 99", y=1.5), angle=90, size=3.5)
plot0 <- plot + scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SC.pdf", width=pdf.width, height=pdf.height)

plot.w <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                               col="blue", width=0.5, linetype=1) +
                 scale_color_manual(name="", values = c("black", "blue"),
                                    guide=guide_legend(override.aes = list(
                                    linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SCwithPI_w.pdf", width=pdf.width, height=pdf.height)

plot1 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb1, ymax=ub1), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SCwithPI_1.pdf", width=pdf.width, height=pdf.height)

plot2 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb2, ymax=ub2), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SCwithPI_2.pdf", width=pdf.width, height=pdf.height)

plot3 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb3, ymax=ub3), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SCwithPI_3.pdf", width=pdf.width, height=pdf.height)


#### Sensitivity analysis #######
# Gaussian bound:
multi <- c(0.25, 0.5, 1, 1.5, 2)
time <- c(1989)

for (l in 1:length(time)) {
  ssc.l.1 <- ssc.r.1 <- c()
  u.mean <- result$u.1[time[l]-1988]
  sig <- sqrt(result$u.2[time[l]-1988])
  sig.seq <- multi*sig
  
  for (s in 1:length(sig.seq)) {
    eps  <- sqrt(-log(alpha/4)*2*(sig.seq[s]^2))
    ssc.l.1[s] <- sc.l.0[t==time[l]] + u.mean - eps
    ssc.r.1[s] <- sc.r.0[t==time[l]] + u.mean + eps
  }
  
  sen.dat <- data.frame(t=c(1:5), lb1=ssc.l.1, ub1=ssc.r.1, 
                        lb=rep(sc.l.0[t==time[l]], 5),
                        ub=rep(sc.r.0[t==time[l]], 5), 
                        lab=as.factor(multi))
  plot <- ggplot() + theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x="sd. of e", y="growth in cigarette sales per capita (%)")
  plot <- plot + geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb1, ymax=ub1), 
                               col="maroon", width=0.2, linetype=5) +
    geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb, ymax=ub), 
                  col="blue", width=0.2, linetype=1) +
    geom_hline(yintercept = y[t==time[l]], linetype=1, size=0.3, alpha=0.8) +
    annotate("text", x=5.4, y=y[t==time[l]]-.4,label="Y(1)", size=3.5) +
    scale_x_discrete(labels=c(parse(text=TeX("$0.25\\hat{\\sigma}$")), 
                              parse(text=TeX("$0.5\\hat{\\sigma}$")),
                              parse(text=TeX("$\\hat{\\sigma}$")), 
                              parse(text=TeX("$1.5\\hat{\\sigma}$")),                                
                              parse(text=TeX("$2\\hat{\\sigma}$"))))
  ggsave(paste("figures/CA_sens_", l, ".pdf", sep=""), width=pdf.width, height=pdf.height)
}


######## Level ############
data <- read.csv("California.csv")

t <- data$year[data$state=="California"]
y <- data$cigsale[data$state=="California"]
x <- data[data$state!="California",1:3]
x <- as.matrix(reshape(x, timevar = "state", idvar="year", direction = "wide"))[,-1]
y <- y/100
x <- x/100

# cutoff: 1988
A <- y[t<=1988]
B <- x[t<=1988,]
C <- NULL
L <- 7
x.T <- x[(t>=1988+1&t<=1988+L),]
eq <- 1
lb <- 0
method.u <- "all"
alpha <- 0.1
M <- 100
result <- sc.pi(A=A, B=B, C=C, x.T=x.T, eq=eq, lb=lb, method.u=method.u, alpha=alpha, M=M, model=3, u.order=1, constant=F)

sc <- c(result$fit.pr, result$fit.po)*100
na <- rep(NA, nrow(B))
sc.l.0  <- c(na, result$sc.l.0)*100; sc.r.0  <- c(na, result$sc.r.0)*100
sc.l.1  <- c(na, result$sc.l.1)*100; sc.r.1  <- c(na, result$sc.r.1)*100
sc.l.2  <- c(na, result$sc.l.2)*100; sc.r.2  <- c(na, result$sc.r.2)*100
sc.l.3  <- c(na, result$sc.l.3)*100; sc.r.3  <- c(na, result$sc.r.3)*100


# Plot
dat    <- data.frame(t=t[t<=1988+L], y=y[t<=1988+L]*100, sname="California")
dat.sc <- data.frame(t=t[t<=1988+L], sc=sc, 
                     lb0=sc.l.0, ub0=sc.r.0, 
                     lb1=sc.l.1, ub1=sc.r.1,
                     lb2=sc.l.2, ub2=sc.r.2,
                     lb3=sc.l.3, ub3=sc.r.3,
                     sname="SC unit")
pdf.width <- 6
pdf.height <- 4.5

plot <- ggplot() + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="year", y="cigarette sales per capita (in packs)") +
  theme(legend.position = c(0.15,0.15), 
        legend.background = element_rect(fill = "transparent"))
plot <- plot + geom_line(data=dat, aes(x=t, y=y, colour=sname), linetype=1) +
  geom_point(data=dat, aes(x=t, y=y, colour=sname), shape=1) +
  geom_line(data=dat.sc, aes(x=t, y=sc, colour=sname), linetype=5) +
  geom_point(data=dat.sc, aes(x=t, y=sc, colour=sname), col="blue", shape=19) +
  geom_vline(xintercept=1988, linetype="dotted") +
  geom_text(aes(x=1988, label="\nProposition 99", y=120), angle=90, size=3.5)
plot0 <- plot + scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                     linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/LevelCA-Tobacco_SC.pdf", width=pdf.width, height=pdf.height)

plot.w <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                               col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/LevelCA-Tobacco_SCwithPI_w.pdf", width=pdf.width, height=pdf.height)

plot1 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb1, ymax=ub1), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/LevelCA-Tobacco_SCwithPI_1.pdf", width=pdf.width, height=pdf.height)

plot2 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb2, ymax=ub2), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/LevelCA-Tobacco_SCwithPI_2.pdf", width=pdf.width, height=pdf.height)

plot3 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb3, ymax=ub3), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb0, ymax=ub0), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/LevelCA-Tobacco_SCwithPI_3.pdf", width=pdf.width, height=pdf.height)


#### Sensitivity analysis #######
# Gaussian bound:
multi <- c(0.25, 0.5, 1, 1.5, 2)
time <- c(1989)

for (l in 1:length(time)) {
  ssc.l.1 <- ssc.r.1 <- c()
  u.mean <- result$u.1[time[l]-1988]
  sig <- sqrt(result$u.2[time[l]-1988])
  sig.seq <- multi*sig
  
  for (s in 1:length(sig.seq)) {
    eps  <- sqrt(-log(alpha/4)*2*(sig.seq[s]^2))
    ssc.l.1[s] <- sc.l.0[t==time[l]] + (u.mean - eps)*100
    ssc.r.1[s] <- sc.r.0[t==time[l]] + (u.mean + eps)*100
  }
  
  sen.dat <- data.frame(t=c(1:5), lb1=ssc.l.1, ub1=ssc.r.1, 
                        lb=rep(sc.l.0[t==time[l]], 5),
                        ub=rep(sc.r.0[t==time[l]], 5), 
                        lab=as.factor(multi))
  plot <- ggplot() + theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x="sd. of e", y="cigarette sales per capita (in packs)")
  plot <- plot + geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb1, ymax=ub1), 
                               col="maroon", width=0.2, linetype=5) +
    geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb, ymax=ub), 
                  col="blue", width=0.2, linetype=1) +
    geom_hline(yintercept = y[t==time[l]]*100, linetype=1, size=0.3, alpha=0.8) +
    annotate("text", x=5.4, y=y[t==time[l]]*100-.4, label="Y(1)", size=3.5) +
    scale_x_discrete(labels=c(parse(text=TeX("$0.25\\hat{\\sigma}$")), 
                              parse(text=TeX("$0.5\\hat{\\sigma}$")),
                              parse(text=TeX("$\\hat{\\sigma}$")), 
                              parse(text=TeX("$1.5\\hat{\\sigma}$")),                                
                              parse(text=TeX("$2\\hat{\\sigma}$"))))
  ggsave(paste("figures/LevelCA_sens_", l, ".pdf", sep=""), width=pdf.width, height=pdf.height)
}