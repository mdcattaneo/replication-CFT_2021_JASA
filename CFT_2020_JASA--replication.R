# Replication file
# Cattaneo, Feng and Rocio (2019)
library(kernlab)
library(pracma)
library(ggplot2)
library(optiSolve)
library(Qtools)
library(latex2exp)
library(Hmisc)
setwd("~/Dropbox/SC/Empirical_Illustration/")


##############################################
############# Summary ########################
data <- read.csv("Germany.csv")
t <- data$year[data$index==1][-1]
y <- data$gdp[data$index==7]
x <- matrix(data$gdp[data$index!=7], 44)
y <- diff(log(data$gdp[data$index==7]))*100
x <- apply(log(x), 2, diff)*100
pre <- as.vector(x[t<=1990,]); pos <- as.vector(x[t>1990&t<=1997,])
pre <- c(summary(pre), sd(pre))
pos <- c(summary(pos), sd(pos))
Germany <- rbind(pre, pos)

data <- read.csv("Israel.csv")
t <- 1:40
y <- data$Israel
x <- as.matrix(data[,-c(1, 13)])
pre <- as.vector(x[t<=19,]); pos <- as.vector(x[t>19&t<=40,])
pre <- c(summary(pre), sd(pre))
pos <- c(summary(pos), sd(pos))
ISR <- rbind(pre, pos)

data <- read.csv("California.csv")
t <- data$year[data$state=="California"][-1]
y <- data$cigsale[data$state=="California"]
x <- data[data$state!="California",1:3]
x <- as.matrix(reshape(x, timevar = "state", idvar="year", direction = "wide"))[,-1]
y <- diff(log(y))*100
x <- apply(log(x), 2, diff)*100
pre <- as.vector(x[t<=1988,]); pos <- as.vector(x[t>1988&t<=1995,])
pre <- c(summary(pre), sd(pre))
pos <- c(summary(pos), sd(pos))
CA  <- rbind(pre, pos)

table <- round(rbind(Germany, ISR, CA), 2)
colheads  <- c("Min.", "1st qu.", "Median", "Mean", "3rd qu.", "Max.", "Std. dev.")
rowname   <- rep(c("pre-treament", "post-treatment"), 3)
n.rgroup  <- rep(2, 3)
rgroup    <- c("Germany", "Israel", "California")
latex(table, file=paste("figures/summary", ".txt", sep = ""), 
      append=FALSE, table.env=FALSE, center="none", title="",
      n.cgroup=NULL, cgroup=NULL, colheads=colheads,
      n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname)
############################################################



###################################################
########## Germany ################################
rm(list=ls())

data <- read.csv("Germany.csv")
t <- data$year[data$index==1][-1]
y <- data$gdp[data$index==7]
x <- matrix(data$gdp[data$index!=7], 44)
y <- diff(log(data$gdp[data$index==7]))*100
x <- apply(log(x), 2, diff)*100

X <- C <- x[t<=1990,]; Y <- d <- y[t<=1990]
n <- nrow(C); p <- ncol(C)
xbar <- colMeans(C); ybar <- mean(d)
C <- scale(C, scale = F); d <- d - ybar
Aeq <- rep(1, p); beq <- 1

H <- crossprod(C); c <- -t(d) %*% C
opt <- ipop(c=c, H=H, A=Aeq, b=1, r=0, l=rep(0, p), u=rep(1, p), sigf=5)
w.hat <- w <- primal(opt)
res <- d - C %*% w.hat
fit.pre <- C %*% w + ybar

# Shrink w on purpose
# rule of thumb threshold:
eta <- sqrt(mean(res^2)*log(p)/n)/min(apply(C, 2, sd))
w[w < eta] <- 0; index <- which(w!=0); s <- sum(w!=0)
Gram  <- crossprod(C)/n

xu <- C * c(res); gamma <- scale(xu, scale = F)
#eta.l <- (log(p)/n)^(3/2)
eta.l <- 1/n
sc0 <- sc0.l <- sc0.u <- sc0.ll.3 <- sc0.uu.3 <- sc0.ll.4 <- sc0.uu.4 <- c()

# set level
level <- 0.68; alpha <- (1-level)/4
L <- 7

# conditional variance  
ufit   <- lm(res~X[,index]); u2fit <- lm(log((res-ufit$fit)^2)~X[,index])
res.st <- (res-ufit$fit)/sqrt(exp(u2fit$fit))
# quantile reg
qfit <- rrq(res~X[,index], tau=c(alpha, 1-alpha))

for (l in 1:L) {
  x.new <- x[t==1990+l,]
  x.T <- (x.new-xbar)
  lb <- ub <- c()
  
  # define a linear objective fn
  obj <- linfun(a=x.T, d=-sum(x.T * w))
  # boot NEGATIVE support
  for (i in 1:100) {
    zeta    <- rnorm(n)
    gamma.b <- colMeans(gamma * zeta)
    
    # define quadratic constraint
    qcon <- quadcon(Q=Gram, a=-2*gamma.b-2*c(t(w)%*%Gram), 
                    d=2*sum(gamma.b*w)+sum(w*(Gram %*% w)), val=eta.l)
    # define optimization; min
    co <- cop(f=obj, max=F, lb=lbcon(rep(0,p)), 
              lc=lincon(A=t(rep(1, p)), val=sum(w), name="eq"), qc=qcon)
    result <- solvecop(co, solver="cccp", quiet=T)
    ub[i] <- -validate(co, result, quiet=T)$obj.fun
    
    # define optimization; max
    co <- cop(f=obj, max=T, lb=lbcon(rep(0,p)), 
              lc=lincon(A=t(rep(1, p)), val=sum(w), name="eq"), qc=qcon)
    result <- solvecop(co, solver="cccp", quiet=T)
    lb[i] <- -validate(co, result, quiet=T)$obj.fun
  }
  sc0[l]   <- sum(x.T * w.hat) + ybar
  sc0.l[l] <- sc0[l] + quantile(lb, alpha)
  sc0.u[l] <- sc0[l] + quantile(ub, 1-alpha)
  
  # Adjust error u.T
  # 3rd approach: conditional mean and variance
  u.T.mean <- sum(c(1, x.new[index]) * ufit$coeff)
  u.T.sig  <- sqrt(exp(sum(c(1,x.new[index])*u2fit$coeff)))
  sc0.ll.3[l] <- sc0.l[l] + u.T.mean + u.T.sig*quantile(res.st, alpha)
  sc0.uu.3[l] <- sc0.u[l] + u.T.mean + u.T.sig*quantile(res.st, 1-alpha)
  
  # 4th approach: quantile reg
  sc0.ll.4[l] <- sc0.l[l] + sum(qfit$coefficients[,1]*c(1, x.new[index]))
  sc0.uu.4[l] <- sc0.u[l] + sum(qfit$coefficients[,2]*c(1, x.new[index]))
}

# Adjust error terms
# When there is an intercept, we should consider u_T-mean(u)
# Gaussian bound:
eps <- sqrt(-log(2*alpha/2)*2*mean(res^2)*(1+1/n))
sc0.ll.1 <- sc0.l - eps
sc0.uu.1 <- sc0.u + eps
# Alternative: use residual quantile
sc0.ll.2 <- sc0.l + quantile(res, alpha)
sc0.uu.2 <- sc0.u + quantile(res, 1-alpha)

sc <- c(fit.pre, sc0)
sc.l <- c(rep(NA, n), sc0.l);       sc.u  <- c(rep(NA, n),sc0.u)
sc.ll.1 <- c(rep(NA, n), sc0.ll.1); sc.uu.1 <- c(rep(NA, n),sc0.uu.1)
sc.ll.2 <- c(rep(NA, n), sc0.ll.2); sc.uu.2 <- c(rep(NA, n),sc0.uu.2)
sc.ll.3 <- c(rep(NA, n), sc0.ll.3); sc.uu.3 <- c(rep(NA, n),sc0.uu.3)
sc.ll.4 <- c(rep(NA, n), sc0.ll.4); sc.uu.4 <- c(rep(NA, n),sc0.uu.4)

# Plot
dat    <- data.frame(t=t[t<=1990+L], y=y[t<=1990+L], sname="Germany")
dat.sc <- data.frame(t=t[t<=1990+L], sc=sc, 
                     lb=sc.l,    ub=sc.u, 
                     lb1=sc.ll.1, ub1=sc.uu.1,
                     lb2=sc.ll.2, ub2=sc.uu.2,
                     lb3=sc.ll.3, ub3=sc.uu.3,
                     lb4=sc.ll.4, ub4=sc.uu.4,
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

plot.w <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                               col="blue", width=0.5, linetype=1) +
                 scale_color_manual(name="", values = c("black", "blue"),
                                    guide=guide_legend(override.aes = list(
                                    linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/Germany_SCwithPI_w.pdf", width=pdf.width, height=pdf.height)

plot1 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb1, ymax=ub1), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/Germany_SCwithPI_1.pdf", width=pdf.width, height=pdf.height)

plot2 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb2, ymax=ub2), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/Germany_SCwithPI_2.pdf", width=pdf.width, height=pdf.height)

plot3 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb3, ymax=ub3), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/Germany_SCwithPI_3.pdf", width=pdf.width, height=pdf.height)

plot4 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb4, ymax=ub4), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/Germany_SCwithPI_4.pdf", width=pdf.width, height=pdf.height)


#### Sensitivity analysis #######
# Gaussian bound:
sig2    <- sqrt(mean(res^2)); multi <- c(0.25, 0.5, 1, 1.5, 2)
sig.seq <- multi*sig2

time <- c(1, 3)

for (l in 1:length(time)) {
  ssc.ll.1 <- ssc.uu.1 <- c()
  for (s in 1:length(sig.seq)) {
    eps  <- sqrt(-log(2*alpha/2)*2*(sig.seq[s]^2)*(1+1/n))
    ssc.ll.1[s] <- sc0.l[time[l]] - eps
    ssc.uu.1[s] <- sc0.u[time[l]] + eps
  }
  
  sen.dat <- data.frame(t=c(1:5), lb1=ssc.ll.1, ub1=ssc.uu.1, 
                        lb=rep(sc0.l[time[l]], 5),
                        ub=rep(sc0.u[time[l]], 5), 
                        lab=as.factor(multi))
  plot <- ggplot() + theme_bw() + 
                     theme(panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank()) +
                           labs(x="sd. of u", y="growth in GDP per capita (%)")
  plot <- plot + geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb1, ymax=ub1), 
                               col="maroon", width=0.2, linetype=5) +
                 geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb, ymax=ub), 
                               col="blue", width=0.2, linetype=1) +
                 geom_hline(yintercept = y[t==1990+time[l]], linetype=1, size=0.3, alpha=0.8) +
                 annotate("text", x=5.4, y=y[t==1990+time[l]]-.4,label="Y(1)", size=3.5) +
                 scale_x_discrete(labels=c(parse(text=TeX("$0.25\\hat{\\sigma}$")), 
                                           parse(text=TeX("$0.5\\hat{\\sigma}$")),
                                           parse(text=TeX("$\\hat{\\sigma}$")), 
                                           parse(text=TeX("$1.5\\hat{\\sigma}$")),                                
                                           parse(text=TeX("$2\\hat{\\sigma}$"))))
  ggsave(paste("figures/Germany_sens_", l, ".pdf", sep=""), width=pdf.width, height=pdf.height)
}


##########################################################
##################### Quarterly Israel ######################
##########################################################
rm(list=ls())
data <- read.csv("Israel.csv")
t <- 1:40
y <- data$Israel
x <- as.matrix(data[,-c(1, 13)])
# cutoff: 19
X <- C <- x[t<=19,]; Y <- d <- y[t<=19]
n <- nrow(C); p <- ncol(C)
xbar <- colMeans(C); ybar <- mean(d)
C <- scale(C, scale = F); d <- d - ybar
Aeq <- rep(1, p); beq <- 1

H <- crossprod(C); c <- -t(d) %*% C
opt <- ipop(c=c, H=H, A=Aeq, b=1, r=0, l=rep(0, p), u=rep(1, p), sigf=5)
w.hat <- w <- primal(opt)
res <- d - C %*% w.hat
fit.pre <- scale(x[t<=19,], scale = F) %*% w + ybar

# Shrink w on purpose
# rule of thumb threshold:
eta <- 1/sqrt(n)
w[w < eta] <- 0; index <- which(w!=0); s <- sum(w!=0)
Gram  <- crossprod(C)/n

xu <- C * c(res); gamma <- scale(xu, scale = F)
eta.l <- 1/n
sc0 <- sc0.l <- sc0.u <- sc0.ll.3 <- sc0.uu.3 <- sc0.ll.4 <- sc0.uu.4 <- c()

# set level
level <- 0.68; alpha <- (1-level)/4
L <- 21

# conditional mean and variance
ufit <- lm(res~X[,index]);   u2fit <- lm(log((res-ufit$fit)^2)~X[,index])
res.st  <- (res-ufit$fit)/sqrt(exp(u2fit$fit))
# quantile reg
qfit <- rrq(res~X[,index], tau=c(alpha, 1-alpha))

for (l in 1:L) {
  x.new <- x[t==19+l] 
  x.T   <- x.new-xbar
  lb <- ub <- c()
  
  # define a linear objective fn
  obj <- linfun(a=x.T, d=-sum(x.T * w))
  
  # boot NEGATIVE support
  for (i in 1:100) {
    zeta    <- rnorm(n)
    gamma.b <- colMeans(gamma * zeta)
    
    # define quadratic constraint
    qcon <- quadcon(Q=Gram, a=-2*gamma.b-2*c(t(w)%*%Gram), 
                    d=2*sum(gamma.b*w)+sum(w*(Gram %*% w)), val=eta.l)
    # define optimization; min
    co <- cop(f=obj, max=F, lb=lbcon(rep(0,p)), 
              lc=lincon(A=t(rep(1, p)), val=sum(w), name="eq"), qc=qcon)
    result <- solvecop(co, solver="cccp", quiet=T)
    ub[i] <- -validate(co, result, quiet=T)$obj.fun
    
    # define optimization; max
    co <- cop(f=obj, max=T, lb=lbcon(rep(0,p)), 
              lc=lincon(A=t(rep(1, p)), val=sum(w), name="eq"), qc=qcon)
    result <- solvecop(co, solver="cccp", quiet=T)
    lb[i] <- -validate(co, result, quiet=T)$obj.fun
  }
  sc0[l]   <- sum(x.T * w.hat) + ybar
  sc0.l[l] <- sc0[l] + quantile(lb, alpha)
  sc0.u[l] <- sc0[l] + quantile(ub, 1-alpha)
  
  # Adjust error u.T
  # 3rd approach: conditional
  # conditional mean of u
  u.T.mean <- sum(c(1, x.new[index]) * ufit$coeff)
  u.T.sig <- sqrt(exp(sum(c(1,x.new[index])*u2fit$coeff)))
  sc0.ll.3[l] <- sc0.l[l] + u.T.mean + u.T.sig*quantile(res.st, alpha)
  sc0.uu.3[l] <- sc0.u[l] + u.T.mean + u.T.sig*quantile(res.st, 1-alpha)
  
  # 4th approach: quantile reg
  sc0.ll.4[l] <- sc0.l[l] + sum(qfit$coefficients[,1]*c(1, x.new[index]))
  sc0.uu.4[l] <- sc0.u[l] + sum(qfit$coefficients[,2]*c(1, x.new[index]))
}

# Adjust error terms
# When there is an intercept, we should consider u_T-mean(u)
# Gaussian bound:
eps <- sqrt(-log(2*alpha/2)*2*mean(res^2)*(1+1/n))
sc0.ll.1 <- sc0.l - eps
sc0.uu.1 <- sc0.u + eps
# Alternative: use residual quantile
sc0.ll.2 <- sc0.l + quantile(res, alpha)
sc0.uu.2 <- sc0.u + quantile(res, 1-alpha)

sc <- c(fit.pre, sc0)
sc.l <- c(rep(NA, n), sc0.l);       sc.u  <- c(rep(NA, n),sc0.u)
sc.ll.1 <- c(rep(NA, n), sc0.ll.1); sc.uu.1 <- c(rep(NA, n),sc0.uu.1)
sc.ll.2 <- c(rep(NA, n), sc0.ll.2); sc.uu.2 <- c(rep(NA, n),sc0.uu.2)
sc.ll.3 <- c(rep(NA, n), sc0.ll.3); sc.uu.3 <- c(rep(NA, n),sc0.uu.3)
sc.ll.4 <- c(rep(NA, n), sc0.ll.4); sc.uu.4 <- c(rep(NA, n),sc0.uu.4)

# Plot
dat    <- data.frame(t=t[t<=19+L], y=y[t<=19+L], sname="Israel")
dat.sc <- data.frame(t=t[t<=19+L], sc=sc, 
                     lb=sc.l,      ub=sc.u, 
                     lb1=sc.ll.1,  ub1=sc.uu.1,
                     lb2=sc.ll.2,  ub2=sc.uu.2,
                     lb3=sc.ll.3,  ub3=sc.uu.3,
                     lb4=sc.ll.4,  ub4=sc.uu.4,
                     sname="SC unit")
pdf.width <- 6
pdf.height <- 4.5

plot <- ggplot() + theme_bw() + 
  labs(x="quarter", y="growth in quarterly GDP (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.15,0.15), 
        legend.background = element_rect(fill = "transparent"))
plot <- plot + geom_line(data=dat, aes(x=t, y=y, colour=sname), linetype=1) +
  geom_point(data=dat, aes(x=t, y=y, colour=sname), shape=1) +
  geom_line(data=dat.sc, aes(x=t, y=sc, colour=sname), linetype=5) +
  geom_point(data=dat.sc, aes(x=t, y=sc, colour=sname), col="blue", shape=19)+
  geom_vline(xintercept=19, linetype="dotted") +
  geom_text(aes(x=19, label="\nPalestinian Intifada", y=-1), angle=90, size=3)
plot0 <- plot + scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                     linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/QGDP_Israel_SC.pdf", width=pdf.width, height=pdf.height)

plot.w <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                               col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/QGDP_Israel_SCwithPI_w.pdf", width=pdf.width, height=pdf.height)

plot1 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb1, ymax=ub1), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/QGDP_Israel_SCwithPI_1.pdf", width=pdf.width, height=pdf.height)

plot2 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb2, ymax=ub2), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/QGDP_Israel_SCwithPI_2.pdf", width=pdf.width, height=pdf.height)

plot3 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb3, ymax=ub3), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/QGDP_Israel_SCwithPI_3.pdf", width=pdf.width, height=pdf.height)

plot4 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb4, ymax=ub4), 
                              col="maroon", width=0.5, linetype=5) +
  geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                col="blue", width=0.5, linetype=1) +
  scale_color_manual(name="", values = c("black", "blue"),
                     guide=guide_legend(override.aes = list(
                       linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/QGDP_Israel_SCwithPI_4.pdf", width=pdf.width, height=pdf.height)

#### Sensitivity analysis #######
# Gaussian bound:
sig2    <- sqrt(mean(res^2)); multi <- c(0.25, 0.5, 1, 1.5, 2)
sig.seq <- multi*sig2

time <- c(2, 5)

for (l in 1:length(time)) {
  ssc.ll.1 <- ssc.uu.1 <- c()
  for (s in 1:length(sig.seq)) {
    eps  <- sqrt(-log(2*alpha/2)*2*(sig.seq[s]^2)*(1+1/n))
    ssc.ll.1[s] <- sc0.l[time[l]] - eps
    ssc.uu.1[s] <- sc0.u[time[l]] + eps
  }
  
  sen.dat <- data.frame(t=c(1:5), lb1=ssc.ll.1, ub1=ssc.uu.1,
                        lb=rep(sc0.l[time[l]], 5),
                        ub=rep(sc0.u[time[l]], 5), 
                        lab=as.factor(multi))
  plot <- ggplot() + theme_bw() + 
                     theme(panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank()) +
                           labs(x="sd. of u", y="growth in quarterly GDP (%)")
  plot <- plot + geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb1, ymax=ub1), 
                               col="maroon", width=0.2, linetype=5) +
                 geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb, ymax=ub), 
                               col="blue", width=0.2, linetype=1) +
                 geom_hline(yintercept = y[t==19+time[l]], linetype=1, size=.3, alpha=.8) +
                 annotate("text", x=5.45, y=y[t==19+time[l]]-.4,label="Y(1)", size=3.5) +
                 scale_x_discrete(labels=c(parse(text=TeX("$0.25\\hat{\\sigma}$")), 
                                           parse(text=TeX("$0.5\\hat{\\sigma}$")),
                                           parse(text=TeX("$\\hat{\\sigma}$")), 
                                           parse(text=TeX("$1.5\\hat{\\sigma}$")),                                
                                           parse(text=TeX("$2\\hat{\\sigma}$"))))
  ggsave(paste("figures/Israel_sens_", l, ".pdf", sep=""), width=pdf.width, height=pdf.height)
}


############################################################################
######### California Tobacco Control #######################################
############################################################################
rm(list=ls())
data <- read.csv("California.csv")

t <- data$year[data$state=="California"][-1]
y <- data$cigsale[data$state=="California"]
x <- data[data$state!="California",1:3]
x <- as.matrix(reshape(x, timevar = "state", idvar="year", direction = "wide"))[,-1]
y <- diff(log(y))*100
x <- apply(log(x), 2, diff)*100

X <- C <- x[t<=1988,]; Y <- d <- y[t<=1988]; n <- nrow(C); p <- ncol(C)
xbar <- colMeans(C); ybar <- mean(d)
C <- scale(C, scale = F); d <- d - ybar
Aeq <- rep(1, p); beq <- 1

H <- crossprod(C); c <- -t(d) %*% C
opt <- ipop(c=c, H=H, A=Aeq, b=1, r=0, l=rep(0, p), u=rep(1, p), sigf=5)
w.hat <- w <- primal(opt)
res <- d - C %*% w.hat
fit.pre <- C %*% w + ybar

# Shrink w on purpose
# rule of thumb threshold:
eta <- sqrt(mean(res^2)*log(p)/n)/min(apply(C, 2, sd))
w[w < eta] <- 0; index <- which(w!=0); s <- sum(w!=0)
Gram  <- crossprod(C)/n

xu <- C * c(res); gamma <- scale(xu, scale = F)
eta.l <- 1/n
sc0 <- sc0.l <- sc0.u <- sc0.ll.3 <- sc0.uu.3 <- sc0.ll.4 <- sc0.uu.4 <- c()

# set level
level <- 0.68; alpha <- (1-level)/4
L <- 7

# conditional variance
ufit <- lm(res~X[,index]); u2fit <- lm(log((res-ufit$fit)^2)~X[,index])
res.st <- (res-ufit$fit)/sqrt(exp(u2fit$fit))
# quantile reg
qfit <- rrq(res~X[,index], tau=c(alpha, 1-alpha))


for (l in 1:L) {
  x.new <- x[t==1988+l]
  x.T <- (x.new-xbar)
  lb <- ub <- c()
  
  # define a linear objective fn
  obj <- linfun(a=x.T, d=-sum(x.T * w))
  
  # boot NEGATIVE support
  for (i in 1:100) {
    zeta    <- rnorm(n)
    gamma.b <- colMeans(gamma * zeta)
    
    # define quadratic constraint
    qcon <- quadcon(Q=Gram, a=-2*gamma.b-2*c(t(w)%*%Gram), 
                    d=2*sum(gamma.b*w)+sum(w*(Gram %*% w)), val=eta.l)
    # define optimization; min
    co <- cop(f=obj, max=F, lb=lbcon(rep(0,p)), 
              lc=lincon(A=t(rep(1, p)), val=sum(w), name="eq"), qc=qcon)
    result <- solvecop(co, solver="cccp", quiet=T)
    ub[i] <- -validate(co, result, quiet=T)$obj.fun
    
    # define optimization; max
    co <- cop(f=obj, max=T, lb=lbcon(rep(0,p)), 
              lc=lincon(A=t(rep(1, p)), val=sum(w), name="eq"), qc=qcon)
    result <- solvecop(co, solver="cccp", quiet=T)
    lb[i] <- -validate(co, result, quiet=T)$obj.fun
  }
  sc0[l]   <- sum(x.T * w.hat) + ybar
  sc0.l[l] <- sc0[l] + quantile(lb, alpha)
  sc0.u[l] <- sc0[l] + quantile(ub, 1-alpha)
  
  # Adjust error u.T
  # 3rd approach: conditional mean and variance of u
  u.T.mean <- sum(c(1, x.new[index]) * ufit$coeff)
  u.T.sig <- sqrt(exp(sum(c(1,x.new[index])*u2fit$coeff)))
    
  sc0.ll.3[l] <- sc0.l[l] + u.T.mean + u.T.sig*quantile(res.st, alpha)
  sc0.uu.3[l] <- sc0.u[l] + u.T.mean + u.T.sig*quantile(res.st, 1-alpha)
  
  # 4th approach: quantile reg
  sc0.ll.4[l] <- sc0.l[l] + sum(qfit$coefficients[,1]*c(1, x.new[index]))
  sc0.uu.4[l] <- sc0.u[l] + sum(qfit$coefficients[,2]*c(1, x.new[index]))
}

# Adjust error terms
# When there is an intercept, we should consider u_T-mean(u)
# Gaussian bound:
eps <- sqrt(-log(2*alpha/2)*2*mean(res^2)*(1+1/n))
sc0.ll.1 <- sc0.l - eps
sc0.uu.1 <- sc0.u + eps
# Alternative: use residual quantile
sc0.ll.2 <- sc0.l + quantile(res, alpha)
sc0.uu.2 <- sc0.u + quantile(res, 1-alpha)

sc <- c(fit.pre, sc0)
sc.l <- c(rep(NA, n), sc0.l);       sc.u  <- c(rep(NA, n),sc0.u)
sc.ll.1 <- c(rep(NA, n), sc0.ll.1); sc.uu.1 <- c(rep(NA, n),sc0.uu.1)
sc.ll.2 <- c(rep(NA, n), sc0.ll.2); sc.uu.2 <- c(rep(NA, n),sc0.uu.2)
sc.ll.3 <- c(rep(NA, n), sc0.ll.3); sc.uu.3 <- c(rep(NA, n),sc0.uu.3)
sc.ll.4 <- c(rep(NA, n), sc0.ll.4); sc.uu.4 <- c(rep(NA, n),sc0.uu.4)

# Plot
dat    <- data.frame(t=t[t<=1988+L], y=y[t<=1988+L], sname="California")
dat.sc <- data.frame(t=t[t<=1988+L], sc=sc, 
                     lb=sc.l,     ub=sc.u, 
                     lb1=sc.ll.1, ub1=sc.uu.1,
                     lb2=sc.ll.2, ub2=sc.uu.2,
                     lb3=sc.ll.3, ub3=sc.uu.3,
                     lb4=sc.ll.4, ub4=sc.uu.4,
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
               geom_text(aes(x=1988, label="\nProposition 99", y=.5), angle=90, size=3.5)
plot0 <- plot + scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SC.pdf", width=pdf.width, height=pdf.height)

plot.w <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                               col="blue", width=0.5, linetype=1) +
                 scale_color_manual(name="", values = c("black", "blue"),
                                    guide=guide_legend(override.aes = list(
                                    linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SCwithPI_w.pdf", width=pdf.width, height=pdf.height)

plot1 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb1, ymax=ub1), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SCwithPI_1.pdf", width=pdf.width, height=pdf.height)

plot2 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb2, ymax=ub2), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SCwithPI_2.pdf", width=pdf.width, height=pdf.height)

plot3 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb3, ymax=ub3), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SCwithPI_3.pdf", width=pdf.width, height=pdf.height)

plot4 <- plot + geom_errorbar(data=dat.sc, aes(x=t, ymin=lb4, ymax=ub4), 
                              col="maroon", width=0.5, linetype=5) +
                geom_errorbar(data=dat.sc, aes(x=t, ymin=lb, ymax=ub), 
                              col="blue", width=0.5, linetype=1) +
                scale_color_manual(name="", values = c("black", "blue"),
                                   guide=guide_legend(override.aes = list(
                                   linetype=c(1,5), shape=c(1, 19))))
ggsave("figures/CA-Tobacco_SCwithPI_4.pdf", width=pdf.width, height=pdf.height)

#### Sensitivity analysis #######
# Gaussian bound:
sig2    <- sqrt(mean(res^2)); multi <- c(0.25, 0.5, 1, 1.5, 2)
sig.seq <- multi*sig2

time <- c(2, 3)

for (l in 1:length(time)) {
  ssc.ll.1 <- ssc.uu.1 <- c()
  for (s in 1:length(sig.seq)) {
    eps  <- sqrt(-log(2*alpha/2)*2*(sig.seq[s]^2)*(1+1/n))
    ssc.ll.1[s] <- sc0.l[time[l]] - eps
    ssc.uu.1[s] <- sc0.u[time[l]] + eps
  }
  
  sen.dat <- data.frame(t=c(1:5), lb1=ssc.ll.1, ub1=ssc.uu.1, 
                        lb=rep(sc0.l[time[l]], 5),
                        ub=rep(sc0.u[time[l]], 5), 
                        lab=as.factor(multi))
  plot <- ggplot() + theme_bw() + 
                     theme(panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank()) +
                     labs(x="sd. of u", y="growth in cigarette sales per capita (%)")
  plot <- plot + geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb1, ymax=ub1), 
                               col="maroon", width=0.2, linetype=5) +
                 geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb, ymax=ub), 
                               col="blue", width=0.2, linetype=1) +
                 geom_hline(yintercept = y[t==1988+time[l]], linetype=1, size=.3, alpha=.8) +
                 annotate("text", x=5.4, y=y[t==1988+time[l]]-.4,label="Y(1)", size=3.5) +
                 scale_x_discrete(labels=c(parse(text=TeX("$0.25\\hat{\\sigma}$")), 
                                           parse(text=TeX("$0.5\\hat{\\sigma}$")),
                                           parse(text=TeX("$\\hat{\\sigma}$")), 
                                           parse(text=TeX("$1.5\\hat{\\sigma}$")),                                
                                           parse(text=TeX("$2\\hat{\\sigma}$"))))
  ggsave(paste("figures/CA_sens_", l, ".pdf", sep=""), width=pdf.width, height=pdf.height)
}

