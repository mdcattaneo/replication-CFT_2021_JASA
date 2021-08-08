# Replication file: empirical illustration
# Cattaneo, Feng and Titiunik (2021)
# Supporting functions
# Notes: this file contains a subset of functions defined for simulation purpose (with several small changes)
# Date: Aug 7, 2021


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
sc.pi <- function(A, B, C=NULL, x.T, eq=1, lb=0, method.u, alpha, M=200, model, u.order=1, 
                  constant=F, rho.max=1, ini=0, vce="HC0") {
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
  sc.pred <- sc.est(A=A, B=B, C=C.full, x.T=x.T.full, eq=eq, lb=lb)
  w.star  <- w <- sc.pred$w
  fit.po  <- sc.pred$fit.po 
  res     <- sc.pred$res
  #print(w)
  
  # inference
  # truncate
  Z <- cbind(B, C.full)
  ind.B <- w[1:d.B] > 0.01
  B2 <- B[,ind.B,drop=F]*B[,ind.B,drop=F]
  if (model==3) {
    # cointegration
    #rho <- sqrt(mean((res-mean(res))^2))*log(T0)/T0 * sqrt(max(colSums(B2))/T0^2)/(min(colSums(B2))/T0^2)
    rho <- sqrt(mean((res-mean(res))^2))*log(T0)/T0 /sqrt(min(colSums(B2))/T0^2)
    #rho <- sqrt(mean((res-mean(res))^2))*log(T0)/T0/(min(apply(Z[,1:d.B], 2, sd))/sqrt(T0))
  } else if (model==1) {
    # iid
    #rho <- sqrt(mean((res-mean(res))^2)*log(T0)/T0) * sqrt(max(colSums(B2)/T0))/(min(colSums(B2)/T0))
    rho <- sqrt(mean((res-mean(res))^2)*log(T0)/T0) /sqrt(min(colSums(B*B)/T0))
    #rho <- sqrt(mean((res-mean(res))^2)*log(T0)/T0)/min(apply(Z[,1:d.B], 2, sd))
  } else if (model==2) {
    # AR
    #rho <- sqrt(mean((res-mean(res))^2)*log(T0)/T0) * sqrt(max(colSums(B2)/T0))/(min(colSums(B2)/T0))
    rho <- sqrt(mean((res-mean(res))^2)*log(T0)/T0) /sqrt(min(colSums(B*B)/T0))
    #rho <- sqrt(mean((res-mean(res))^2)*log(T0)/T0)/min(apply(Z[,1:d.B], 2, sd))
  }
  # upper bound
  #print(rho); print(max(w))
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
      x.diff <- apply(rbind(B, x.T[,1:d.B,drop=F]), 2, diff)
      des.0  <- cbind(x.diff[1:(T0-1),,drop=F], C[1:(T0-1),,drop=F])[,index.sub,drop=F]
      des.1  <- cbind(x.diff[-(1:(T0-1)),,drop=F], x.T[,-(1:d.B),drop=F])[,index.sub,drop=F]
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
  # save mean and variance of u, only for sensitivity analysis
  u.1 <- u.2 <- NA
  
  if (method.u=="gaussian" | method.u=="all") {
    pi.u   <- sc.pi.u(res=res, x=des.0, eval=des.1, method="gaussian", alpha=alpha/4)
    lb.u.1 <- pi.u$lb
    ub.u.1 <- pi.u$ub
    sc.l.1 <- sc.l.0 + pi.u$lb
    sc.r.1 <- sc.r.0 + pi.u$ub
    len.1  <- sc.r.1 - sc.l.1
    u.1 <- pi.u$u.1
    u.2 <- pi.u$u.2
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
              len.2=len.2,           len.3=len.3,  len.4=len.4,
              u.1=u.1,               u.2=u.2))
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
  u.1 <- u.2 <- NA
  
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
    
    # save mean and variance of u, only for sensitivity analysis
    u.1 <- u.mean
    u.2 <- u.sig2
    
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
  
  return(list(lb=lb, ub=ub, u.1=u.1, u.2=u.2))
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

