# Replication: simulation
# Cattaneo, Feng and Titiunik (2021)
# Supporting functions
# Date: Aug 7, 2021

# DGP, generate design
# T0: sample size; ar: AR(1) coeff mat;
# return pre-treatment control units
dgp.x <- function(T0, ar) {
  y.co <- VAR.sim(B=ar, n=T0+1, include="none")
  return(y.co)
}

# DGP, conditional sampling process
# y.co.0: T0 by J, pretreatment covariates; y.co.1: post-treatment covariates, neval by J
# y.tr.1: neval by 1; w0.u, w0.c should be vectors
dgp.cond <- function(y.co.0, y.co.1, T0, J, model=1, err=F, w0.u, w0.c=NULL) {
  zeta <- rnorm(T0+1, sd=0.5)
  if (is.vector(y.co.1)) {
    y.co.1 <- t(y.co.1)
  }
  
  # gen misspecification error + idn. error
  if (err) {
    if (model==3) {
      u.0  <- 0.9*c(y.co.0[1,1], diff(y.co.0[,1])) + zeta[1:T0]
      u.1  <- 0.9*(y.co.1[,1]-y.co.0[T0,1])        + zeta[T0+1]
    } else {
      u.0  <- 0.2*y.co.0[,1] + zeta[1:T0]
      u.1  <- 0.2*y.co.1[,1] + zeta[T0+1]
    }  
  } else {
    u.0  <- zeta[1:T0]
    u.1  <- zeta[T0+1]
  }
  
  xw.0   <- y.co.0 %*% w0.u
  y.tr.0 <- xw.0 + u.0
  xw.1   <- y.co.1 %*% w0.u
  y.tr.1 <- xw.1 + u.1
  
  if (err) {
    xw0.1 <- y.co.1 %*% w0.c
  } else {
    xw0.1 <- xw.1
  }
  
  return(list(y.co.0=y.co.0, y.co.1=y.co.1, 
              y.tr.0=y.tr.0, y.tr.1=y.tr.1,
              xw0.1=xw0.1))         # xw0.1 true conditional x*w0
}

# SC estimation (always J units + covariates); x.T can be a matrix, neval * (d.B+d.C)
sc.est <- function(A, B, C=NULL, x.T, eq=1, lb=0, ini=0) {
  Z <- cbind(B, C)
  d.B <- ncol(B)
  d.C <- 0
  if (!is.null(C)) d.C <- ncol(C)
  
  Gram <- crossprod(Z); a <- -2 * t(A) %*% Z
  
  # define a quad objective fn
  obj <- quadfun(Q=Gram, a=a)
  
  # Lower bound on w
  lb.opt  <- lbcon(val = c(rep(lb, d.B), rep(-Inf, d.C)),
                   id = c(1:(d.B+d.C)))
  
  # Constraint on the norm of w
  roweq       <- t(c(rep(1, d.B), rep(0, d.C)))
  norm.constr <- lincon(A = roweq, dir = c("=="), val = c(eq), name = "Norm-1")
  
  # roweq   <- c(rep(1, d.B), rep(0, d.C))
  # rowineq <- cbind(diag(d.B), matrix(0, d.B, d.C))
  # Mat.cst <- rbind(roweq, rowineq)
  # lc      <- lincon(A=Mat.cst, dir=c("==", rep(">=", d.B)), 
  #                   val=c(eq, rep(lb, d.B)), name=1:nrow(Mat.cst))
  
  # define optimization; min
  co <- cop(f=obj, max=F, lc=norm.constr, lb=lb.opt)
  
  # optimization
  x0        <- rep(ini, (d.B+d.C))
  names(x0) <- co$id
  result <- suppressMessages(solvecop(co, solver="slsqp", quiet=T, X = x0))
  
  w      <- result$x
  fit.pr <- Z %*% w  
  res    <- A - fit.pr
  fit.po <- x.T %*% w

  return(list(w=w, fit.pr=fit.pr, fit.po=fit.po, res=res))
}


# SC prediction intervals
sc.pi <- function(A, B, C=NULL, x.T, eq=1, lb=0, method.u, alpha, M=200, model, 
                  u.order=1, constant=F, rho.max=1, ini=0, vce="HC0") {
  if (is.vector(x.T)) {
    x.T <- t(x.T)
  }
  T0 <- nrow(B); d.B <- ncol(B); neval <- nrow(x.T)
  
  C.full <- C; x.T.full <- x.T
  
  if (constant) {
    C.full <- cbind(C.full, rep(1,T0))
    x.T.full <- cbind(x.T.full, rep(1, neval))
  }
  
  d.C <- 0
  if (!is.null(C.full)) d.C <- ncol(C.full)
  
  # estimation
  sc.pred <- sc.est(A=A, B=B, C=C.full, x.T=x.T.full, eq=eq, lb=lb, ini=ini)
  w.star  <- w <- sc.pred$w
  fit.po  <- sc.pred$fit.po 
  res     <- sc.pred$res
  #print(w)
  # inference
  # truncate
  Z <- cbind(B, C.full)
  if (model==3) {
    # cointegration
    rho <- sqrt(mean((res-mean(res))^2))*log(T0)/T0/(sqrt(min(colSums(B*B)))/T0)
    #rho <- sqrt(mean((res-mean(res))^2))*log(T0)/T0/(min(apply(Z[,1:d.B], 2, sd))/sqrt(T0))
  } else if (model==1) {
    # iid
    rho <- sqrt(mean((res-mean(res))^2)*log(T0)/T0)/min(apply(Z[,1:d.B], 2, sd))
  } else if (model==2) {
    # AR
    rho <- sqrt(mean((res-mean(res))^2)*log(T0)/T0)/min(apply(Z[,1:d.B], 2, sd))
  }
  # upper bound
  #print(rho)
  rho <- min(rho, rho.max)
  
  w.star[1:d.B][w.star[1:d.B] < rho] <- 0 
  index <- w.star!=0
  #s <- sum(w.star!=0)
  
  if (model == 3) {
    Q <- crossprod(Z[-1,,drop=F])/(T0-1)
  } else {
    Q <- crossprod(Z)/T0
  }
  
  # linear constraint
  # Constraint on the norm of w
  roweq       <- t(c(rep(1, d.B), rep(0, d.C)))
  norm.constr <- lincon(A = roweq, dir = c("=="), val = c(sum(w.star)), name = "Norm-1")
  
  # roweq   <- c(rep(1, d.B), rep(0, d.C))
  # rowineq <- cbind(diag(d.B), matrix(0, d.B, d.C))
  # Mat.cst <- rbind(roweq, rowineq)
  # lc  <- lincon(A=Mat.cst, dir=c("==", rep(">=", d.B)), 
  #               val=c(sum(w.star), rep(lb, d.B)), name=1:nrow(Mat.cst))
  
  # variance estimate
  if (u.order==0) {
    if (model==3) {
      des.0 <- as.matrix(rep(1, T0-1))
    } else {
      des.0 <- as.matrix(rep(1, T0))
    }
    des.1 <- as.matrix(rep(1, neval))
  } else {
    # design without constant, selected
    if (constant) {
      index.sub <- index[-length(index)]
    } else {
      index.sub <- index
    }
    if (model==3) {
      des.0 <- cbind(apply(B, 2, diff), C)[,index.sub,drop=F] 
      des.1 <- cbind(sweep(x.T[,1:d.B,drop=F], 2, B[T0,]), 
                     x.T[,-(1:d.B),drop=F])[,index.sub,drop=F]
    } else {
      des.0 <- cbind(B, C)[,index.sub, drop=F]
      des.1 <- x.T[,index.sub,drop=F]
    }
    
    if (u.order > 1) {
      des.0 <- poly(des.0, degree=u.order, raw=T, simple=T)
      des.1 <- poly(des.1, degree=u.order, raw=T, simple=T)
    }
    
    if (model == 3) {
      des.0 <- cbind(des.0, rep(1, T0-1))
    } else {
      des.0 <- cbind(des.0, rep(1, T0))
    }
    des.1 <- cbind(des.1, rep(1, neval))
  }
  
  if (model==3) {
    res <- res[-1]
  }
  
  u.mean <- lm.fit(des.0, res)$fitted.values
  
  if (vce=="HC0") {
    vce.weights <- 1
  } else if (vce=="HC1") {
    if (model==3) vce.weights <- (T0-1)/(T0-2-sum(w!=0))
    else          vce.weights <- T0/(T0-1-sum(w!=0))
  }
  
  Omega  <- diag(c((res-u.mean)^2)) * vce.weights
  
  if (model==3) { 
    Sigma  <- t(Z[-1,,drop=F]) %*% Omega %*% Z[-1,,drop=F] / ((T0-1)^2)
  } else {
    Sigma  <- t(Z) %*% Omega %*% Z / (T0^2)
  }
  Sigma.root <- sqrtm(Sigma)
  
  # obj fun
  obj.ls <- apply(x.T.full, 1, function(xt) linfun(a=xt, d=-sum(xt * w.star)))     # x *(w-w.star)
  
  # Lower bound on w
  lb.opt  <- lbcon(val = c(rep(lb, d.B), rep(-Inf, d.C)),
                   id = c(1:(d.B+d.C)))
  
  # simulate
  vsig <- matrix(unlist(sapply(obj.ls, function(obj) simulate.w(w.star=w.star, Q=Q, Sigma.root=Sigma.root,
                                                                obj=obj, lc=norm.constr, lb.opt=lb.opt, M=M, alpha=alpha/4, ini=ini))), 2)
  lb.w <- vsig[1,]
  ub.w <- vsig[2,]
  
  # PIs for w
  sc.l.0 <- fit.po + lb.w
  sc.r.0 <- fit.po + ub.w
  len.0  <- sc.r.0 - sc.l.0
  
  # PIs for u
  lb.u.1 <- ub.u.1 <- lb.u.2 <- ub.u.2 <- lb.u.3 <- ub.u.3 <- lb.u.4 <- ub.u.4 <- rep(NA, neval)
  sc.l.1 <- sc.r.1 <- sc.l.2 <- sc.r.2 <- sc.l.3 <- sc.r.3 <- sc.l.4 <- sc.r.4 <- rep(NA, neval) 
  len.1 <- len.2 <- len.3 <- len.4 <- rep(NA, neval)
  
  if (method.u=="gaussian" | method.u=="all") {
    pi.u   <- sc.pi.u(res=res, x=des.0, eval=des.1, method="gaussian", alpha=alpha/4)
    lb.u.1 <- pi.u$lb
    ub.u.1 <- pi.u$ub
    sc.l.1 <- sc.l.0 + pi.u$lb
    sc.r.1 <- sc.r.0 + pi.u$ub
    len.1  <- sc.r.1 - sc.l.1
  }
  
  if (method.u=="ls" | method.u=="all") {
    pi.u   <- sc.pi.u(res=res, x=des.0, eval=des.1, method="ls", alpha=alpha/4)
    lb.u.2 <- pi.u$lb
    ub.u.2 <- pi.u$ub
    sc.l.2 <- sc.l.0 + pi.u$lb
    sc.r.2 <- sc.r.0 + pi.u$ub
    len.2  <- sc.r.2 - sc.l.2
  }
  
  if (method.u=="qreg" | method.u=="all") {
    if (u.order == 0) {
      sc.l.3 <- sc.l.0 + quantile(res, alpha/4)
      sc.r.3 <- sc.r.0 + quantile(res, 1-alpha/4)
      len.3  <- sc.r.3 - sc.l.3
    } else {
      pi.u <- sc.pi.u(res=res, x=des.0, eval=des.1, method="qreg", alpha=alpha/4)
      lb.u.3 <- pi.u$lb
      ub.u.3 <- pi.u$ub
      sc.l.3 <- sc.l.0 + pi.u$lb
      sc.r.3 <- sc.r.0 + pi.u$ub
      len.3  <- sc.r.3 - sc.l.3
    }
  }
  
  if (method.u=="double" | method.u=="all") {
    pi.u   <- sc.pi.u(res=res, x=des.0, eval=des.1, method="double", alpha=alpha/4)
    lb.u.4 <- pi.u$lb
    ub.u.4 <- pi.u$ub
    sc.l.4 <- sc.l.0 + pi.u$lb
    sc.r.4 <- sc.r.0 + pi.u$ub
    len.4  <- sc.r.4 - sc.l.4
  }
  
  
  return(list(fit.pr=sc.pred$fit.pr, fit.po=fit.po,
              sc.l.0=sc.l.0,         sc.r.0=sc.r.0,
              sc.l.1=sc.l.1,         sc.r.1=sc.r.1,
              sc.l.2=sc.l.2,         sc.r.2=sc.r.2,
              sc.l.3=sc.l.3,         sc.r.3=sc.r.3,
              sc.l.4=sc.l.4,         sc.r.4=sc.r.4,
              len.0=len.0,           len.1=len.1, 
              len.2=len.2,           len.3=len.3, 
              len.4=len.4))
}


# permutation-based intervals
# data: (T0+1) by (J+1)
sc.permute <- function(data.ab, data.c=NULL, eq=1, lb=0, alpha, ini=0) {
  C <- NULL
  x.T.c <- NULL
  if (!is.null(data.c)) {
    C <- data.c[1:T0,,drop=F]
    x.T.c <- data.c[-(1:T0),,drop=F]
  }
  post <- sapply(1:J, function (i) sc.est(A=data.ab[1:T0,i,drop=F], B=data.ab[1:T0, -i,drop=F],
                                          C=C, x.T=cbind(data.ab[-(1:T0),-i,drop=F], x.T.c),
                                          eq=eq, lb=lb, ini=ini)$fit.po)
  if (is.vector(post)) post <- t(post)
  range <- rowQuantiles(post, probs=c(alpha/2, 1-alpha/2))
  if (is.vector(range)) range <- t(range)
  
  return(list(per.l=range[,1], per.r=range[,2]))
}



#######################################
### Simulation funs ###################
# sc simulate, conditional
sc.sim.cond <- function(i, y.co.0, y.co.1, T0, J, model, err, eq=1, lb=0, method.u, alpha, 
                        M=200, w0.u, w0.c, u.order, rho.max=1, vce="HC0") {
  data   <- dgp.cond(y.co.0, y.co.1, T0, J, model, err, w0.u, w0.c)
  result <- sc.pi(A=data$y.tr.0, B=data$y.co.0, C=NULL, x.T=data$y.co.1,
                  eq=eq, lb=lb, method.u=method.u, alpha=alpha, M=M, model=model, 
                  u.order=u.order, rho.max=rho.max, vce=vce)
  
  # this paper
  rej.w <- rej.1 <- rej.2 <- rej.3 <- rej.4 <- NA
  rej.w <- ((data$xw0.1  < result$sc.l.0) | (data$xw0.1  > result$sc.r.0)) * 1 
  #rej.0 <- ((data$y.tr.1 < result$sc.l.0) | (data$y.tr.1 > result$sc.r.0)) * 1 
  rej.1 <- ((data$y.tr.1 < result$sc.l.1) | (data$y.tr.1 > result$sc.r.1)) * 1 
  rej.2 <- ((data$y.tr.1 < result$sc.l.2) | (data$y.tr.1 > result$sc.r.2)) * 1 
  rej.3 <- ((data$y.tr.1 < result$sc.l.3) | (data$y.tr.1 > result$sc.r.3)) * 1
  rej.4 <- ((data$y.tr.1 < result$sc.l.4) | (data$y.tr.1 > result$sc.r.4)) * 1
  
  # Abadie et al. (2010)
  result.per <- sc.permute(data.ab=rbind(cbind(data$y.tr.0, y.co.0), 
                                         cbind(data$y.tr.1, y.co.1)), 
                           data.c=NULL, eq=eq, lb=lb, alpha=alpha)
  rej.per <- ((data$y.tr.1 < result.per$per.l) | (data$y.tr.1 > result.per$per.r)) * 1
  len.per <- result.per$per.r-result.per$per.l
  
  # Chernozhukov et al. (2020)
  # define a grid
  grid <- seq(min(data$y.tr.1)-3*sd(data$y.tr.0), max(data$y.tr.1)+3*sd(data$y.tr.0), 
              length.out = 100)
  result.conf <- moving.block.ci(y.tr.0=data$y.tr.0, data.x=rbind(y.co.0, y.co.1),
                                 T0=T0, alpha=alpha, grid=grid)
  rej.conf <- ((data$y.tr.1 < result.conf$lb) | (data$y.tr.1 > result.conf$ub)) * 1
  len.conf <- result.conf$ub - result.conf$lb
  
  return(c(rej.w, rej.1, rej.2, rej.3, rej.4, rej.per, rej.conf,
           result$len.0, result$len.1, result$len.2, 
           result$len.3, result$len.4, len.per, len.conf))
}

# sc simulate, unconditional
sc.sim.unc <- function(i, x.list, T0, J, model, err, eq=1, lb=0, method.u, alpha, 
                       M=200, w0.u, w0.c.ls, u.order, rho.max=1, vce="HC0") {
  data.x <- x.list[[i]]
  y.co.0 <- data.x[1:T0,,drop=F]
  y.co.1 <- data.x[(T0+1),,drop=F]
  w0.c   <- w0.c.ls[,i]
  data   <- dgp.cond(y.co.0, y.co.1, T0, J, model, err, w0.u, w0.c)
  
  # This paper
  result <- sc.pi(A=data$y.tr.0, B=data$y.co.0, C=NULL, x.T=data$y.co.1,
                  eq=eq, lb=lb, method.u=method.u, alpha=alpha, M=M, model=model, 
                  u.order=u.order, rho.max=rho.max, vce=vce)
  
  rej.w <- rej.1 <- rej.2 <- rej.3 <- rej.4 <- NA
  rej.w <- ((data$xw0.1  < result$sc.l.0) | (data$xw0.1  > result$sc.r.0)) * 1 
  #rej.0 <- ((data$y.tr.1 < result$sc.l.0) | (data$y.tr.1 > result$sc.r.0)) * 1 
  rej.1 <- ((data$y.tr.1 < result$sc.l.1) | (data$y.tr.1 > result$sc.r.1)) * 1 
  rej.2 <- ((data$y.tr.1 < result$sc.l.2) | (data$y.tr.1 > result$sc.r.2)) * 1 
  rej.3 <- ((data$y.tr.1 < result$sc.l.3) | (data$y.tr.1 > result$sc.r.3)) * 1
  rej.4 <- ((data$y.tr.1 < result$sc.l.4) | (data$y.tr.1 > result$sc.r.4)) * 1
  
  # Abadie et al. (2010)
  result.per <- sc.permute(data.ab=rbind(cbind(data$y.tr.0, y.co.0), 
                                         cbind(data$y.tr.1, y.co.1)), 
                           data.c=NULL, eq=eq, lb=lb, alpha=alpha)
  rej.per <- ((data$y.tr.1 < result.per$per.l) | (data$y.tr.1 > result.per$per.r)) * 1
  len.per <- result.per$per.r-result.per$per.l
  
  # Chernozhukov et al. (2020)
  grid <- seq(min(data$y.tr.1)-3*sd(data$y.tr.0), max(data$y.tr.1)+3*sd(data$y.tr.0), 
              length.out = 100)
  result.conf <- moving.block.ci(y.tr.0=data$y.tr.0, data.x=rbind(y.co.0, y.co.1),
                                 T0=T0, alpha=alpha, grid=grid)
  rej.conf <- ((data$y.tr.1 < result.conf$lb) | (data$y.tr.1 > result.conf$ub)) * 1
  len.conf <- result.conf$ub - result.conf$lb
  
  return(c(rej.w, rej.1, rej.2, rej.3, rej.4, rej.per, rej.conf,  
           result$len.0, result$len.1, result$len.2, 
           result$len.3, result$len.4, len.per, len.conf))
}


##########################
## Ancillary functions ###
##########################
# vectorized simulate fun for w
simulate.w <- function(w.star, Q, Sigma.root, obj, lc, lb.opt, M, alpha, ini) {
  vsig <- sapply(1:M, function(i) sc.pi.w(i, w.star=w.star, Q=Q, Sigma.root=Sigma.root, obj=obj, lc=lc, lb.opt=lb.opt, ini=ini))
  lb.w <- quantile(vsig[1,], alpha, na.rm=T)   # account for potential NaN's
  ub.w <- quantile(vsig[2,], 1-alpha, na.rm=T)
  
  return(list(lb.w=lb.w, ub.w=ub.w))
}


# PI for (w_hat-w)
sc.pi.w <- function(i, w.star, Q, Sigma.root, obj, lc, lb.opt, ini) {
  zeta <- rnorm(length(w.star))
  G    <- Sigma.root %*% zeta
  
  # quadratic constraint
  qcon <- quadcon(Q=Q, a=-2*G - 2*c(t(w.star) %*% Q), 
                  d=2*sum(G*w.star) + sum(w.star*(Q %*% w.star)), val=0)
  
  # define optimization; min
  co     <- cop(f=obj, max=F, lc=lc, qc=qcon, lb=lb.opt)
  x0        <- rep(ini, length(w.star))
  names(x0) <- co$id
  result <- suppressMessages(solvecop(co, solver="slsqp", quiet=T, X=x0))
  ub     <- -validate(co, result, quiet=T)$obj.fun
  
  # define optimization; max
  co     <- cop(f=obj, max=T, lc=lc, qc=qcon, lb=lb.opt)
  names(x0) <- co$id
  result <- suppressMessages(solvecop(co, solver="slsqp", quiet=T, X=x0))
  lb     <- -validate(co, result, quiet=T)$obj.fun
  
  return(c(lb, ub))
}


# Prediction interval, for u
sc.pi.u <- function(res, x, eval, method="gaussian", alpha) {
  #res.fit <- NA
  if (is.vector(eval)) {
    eval <- matrix(eval, 1)
  }
  neval <- nrow(eval)
  if (method == "gaussian") {
    x.more   <- rbind(eval, x)
    fit      <- predict(y=res, x=x, eval=x.more, type="lm")
    u.mean   <- fit[1:neval]
    res.fit  <- fit[-(1:neval)]
    var.pred <- predict(y=log((res-res.fit)^2), x=x, eval=x.more, type="lm")
    u.sig2   <- exp(var.pred[1:neval])
    
    eps <- sqrt(-log(alpha)*2*u.sig2)
    lb <- u.mean - eps
    ub <- u.mean + eps
    
  } else if (method == "ls") {
    x.more  <- rbind(eval, x)
    fit     <- predict(y=res, x=x, eval=x.more, type="lm")
    u.mean  <- fit[1:neval]
    res.fit <- fit[-(1:neval)]
    
    var.pred <- predict(y=log((res-res.fit)^2), x=x, eval=x.more, type="lm")
    u.sig    <- sqrt(exp(var.pred[1:neval]))
    res.st   <- (res-res.fit)/sqrt(exp(var.pred[-(1:neval)]))
    
    lb <- u.mean + u.sig * quantile(res.st, alpha)
    ub <- u.mean + u.sig * quantile(res.st, 1-alpha)
    
  } else if (method == "qreg") {
    #res.fit <- lm.fit(x, res)$fitted.values
    u.pred  <- predict(y=res, x=x, eval=eval, type="qreg", tau=c(alpha, 1-alpha))
    lb <- u.pred[,1]
    ub <- u.pred[,2]
  } else if (method == "double") {
    x.more   <- rbind(eval, x)
    fit      <- predict(y=res, x=x, eval=x.more, type="lm")
    u.mean   <- fit[1:neval]
    res.fit  <- fit[-(1:neval)]
    var.pred <- predict(y=log((res-res.fit)^2), x=x, eval=x.more, type="lm")
    u.sig2   <- exp(var.pred[1:neval])
    
    eps <- 1.5*sqrt(-log(alpha)*2*u.sig2)    # 1.5 sigma
    lb <- u.mean - eps
    ub <- u.mean + eps
  }
  
  return(list(lb=lb, ub=ub))
}


# square root matrix
sqrtm <- function(A) {
  decomp <- svd(A)
  decomp$d[decomp$d < 0] <- 0
  rootA  <- decomp$u %*% diag(sqrt(decomp$d)) %*% t(decomp$u)
  return(rootA)
}

# conditional prediction
predict <- function(y, x, eval, type="lm", tau=NULL) {
  if (is.vector(eval)) {
    eval <- matrix(eval, 1)
  }
  
  if (type == "lm") {
     beta <- .lm.fit(x, y)$coeff
  } else if (type == "qreg") {
     if (is.null(tau)) {
       tau <- c(0.05, 0.95)
     }
     beta <- rrq(y~x-1, tau=tau)$coeff
  }
  pred <- eval %*% beta
  return(pred)
}

# get true conditional w0
getw0 <- function(y.co.0, y.co.1, T0, J, model, w0.u, eq, lb) {
  B <- 5000
  A <- matrix(NA, B, J)
  for (i in 1:B) {
    data   <- dgp.cond(y.co.0, y.co.1, T0, J, model, err=TRUE, w0.u, w0.c=w0.u)
    X <- data$y.co.0
    y <- data$y.tr.0
    A[i,] <- -2 * t(y) %*% X
  } 
  
  a <- colMeans(A)
  Gram <- crossprod(y.co.0)
  
  obj <- quadfun(Q=Gram, a=a)
  
  roweq       <- t(c(rep(1, J), rep(0, 0)))
  norm.constr <- lincon(A = roweq, dir = c("=="), val = c(eq), name = "Norm-1")
  
  lb.opt  <- lbcon(val = c(rep(lb, J), rep(-Inf, 0)),
                   id = c(1:(J+0)))
  
  # roweq   <- c(rep(1, J), rep(0, 0))
  # rowineq <- cbind(diag(J), matrix(0, J, 0))
  # Mat.cst <- rbind(roweq, rowineq)
  # lc      <- lincon(A=Mat.cst, dir=c("==", rep(">=", J)), 
  #                   val=c(eq, rep(lb, J)), name=1:nrow(Mat.cst))
  
  co <- cop(f=obj, max=F, lc=norm.constr, lb=lb.opt)
  
  # optimization
  x0        <- rep(0, (J+0))
  names(x0) <- co$id
  
  result <- suppressMessages(solvecop(co, solver="slsqp", quiet=T, X=x0))
  return(result$x)
}

# Check if within a range
check <- function(x, l, r) {
  return(all((colQuantiles(x, probs=0.1) > l) & 
             (colQuantiles(x, probs=0.9) < r)))
}

####################################################################################
# The following functions are modified from Chernozhukov, Wuthrich, and Zhu (2020)
sc <- function(y1,yJ){
  E     <- matrix(1,1,ncol(yJ))
  F     <- 1
  G     <- diag(ncol(yJ))
  H     <- matrix(0,ncol(yJ),1)
  w.hat <- cbind(lsei(A=yJ, B=y1, E=E, F=F, G=G, H=H, type=2)$X)
  u.hat <- y1-yJ %*% w.hat
  return(list(u.hat=u.hat,w.hat=w.hat))
}

# Moving block permutation
# y1, yJ: T0+1 periods
# simplified: T1=q=1
moving.block.q <- function(y1,yJ,T0) {
  T01 <- T0+1 
  u.hat <- sc(y1,yJ)$u.hat
  
  #sub.size  <- 1
  #u.hat     <- c(u.hat, u.hat)
  #S.vec     <- matrix(NA,T01,1)

  S.vec <- abs(u.hat)
  
  # if (q==2){
  #   for (s in 1:(T01)){
  #     S.vec[s,1]  <- (1/sqrt(T1))*sqrt(sum((u.hat[s:(s+sub.size-1)])^2))
  #   }
  # }
  
  ind <- T0+1
  p   <- mean(S.vec>=S.vec[ind])
  return(p)
}

# vectorize
# y.tr.0: T0 periods; y.tr.1: T1 periods; data.x: (T0+T1) by J;
moving.block.pvec <- function(y.tr.0, y.tr.1, T0, T1, y.co.0, data.x) {
  return(sapply(1:T1, function(i) moving.block.q(y1=c(y.tr.0, y.tr.1[i]), yJ=rbind(y.co.0, data.x[T0+i,]), T0)))
}

# grid: predicted values (differ from the original version)
moving.block.ci <- function(y.tr.0, data.x, T0, alpha, grid) {
  T1     <- nrow(data.x) - T0
  y.co.0 <- data.x[1:T0,]
  #p.grid  <- NULL
  # T1 by ngrid
  p.grid <- sapply(1:length(grid), function(i) moving.block.pvec(y.tr.0, rep(grid[i], T1), T0, T1, y.co.0, data.x))
  if (T1==1) {
    p.grid <- matrix(p.grid,1)
  }
  
  # for (g in grid) {
  #   y10       <- c(y.tr.0, g)
  #   sapply(1:T1, function(i) moving.block.q(y10, rbind(y.co.0, data.x[T0+i,]), T0))
  #   p.grid    <- cbind(p.grid, moving.block.q(y10 ,data.x, T0))
  # }
  ci <- apply(p.grid, 1, function(x) range(grid[x > alpha]))
  #ci <- grid[(p.grid > alpha)] 
  return(list(lb=ci[1,], ub=ci[2,]))
}