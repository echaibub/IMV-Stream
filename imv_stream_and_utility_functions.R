library(corpcor)
library(MASS)
library(Tsphere)
library(gplots)


#######################################################
## Multivariate Stream functions
#######################################################


lmgamma <- function(q, x) {
  j <- 1:q
  (q*(q - 1)/4) * log(pi) + sum(lgamma(x + (1 - j)/2))
}



GetMultStreamModelScores3 <- function(Y, U, d, V, lambdas, a.tau, b.tau) {
  n <- nrow(Y)
  q <- ncol(Y)
  nn <- length(d)
  le <- length(lambdas)
  log.scores <- rep(NA, le)
  scores <- rep(NA, le)
  Betas <- vector(mode = "list", length = le)
  Taus <- vector(mode = "list", length = le)
  aux1 <- lmgamma(q, a.tau + n/2) - lmgamma(q, a.tau) - n * q * log(2 * pi * b.tau)/2
  UtY <- crossprod(U, Y)
  YtU <- t(UtY)
  aux0 <- crossprod(Y)/(2 * b.tau)
  for (i in 1:le) { 
    aux2 <- -0.5 * q * sum(log(1 + (d^2)/lambdas[i]))
    aux3 <- diag(1/(1 + lambdas[i]/d^2))
    aux3 <- aux3 %*% UtY
    theta <- aux3/d
    Betas[[i]] <- V %*% theta
    aux3 <- YtU %*% aux3
    aux3 <- aux3/(2 * b.tau)
    aux3 <- diag(q) + aux0 - aux3
    Taus[[i]] <- ((2 * a.tau + n)/(2 * b.tau)) * solve(aux3) 
    aux3 <- -(a.tau + n/2) * as.numeric(determinant(aux3)$modulus)
    log.scores[i] <- aux1 + aux2 + aux3
    scores[i] <- exp(log.scores[i])
  }
  list(log.scores = log.scores, scores = scores, Betas = Betas, Taus = Taus)
}



MultStream3 <- function(Y, U, d, V, lambda.grid, a.tau, b.tau) {
  out <- GetMultStreamModelScores3(Y, U, d, V, lambda.grid, a.tau, b.tau)
  LogScores <- out$log.scores
  Betas <- out$Betas
  Taus <- out$Taus
  ModPost <- ModelProbs(x = LogScores)
  Beta <- Betas[[1]] * ModPost[1]
  Tau <- Taus[[1]] * ModPost[1]
  for (i in 2:length(ModPost)) {
    Beta <- Beta + Betas[[i]] * ModPost[i]
    Tau <- Tau + Taus[[i]] * ModPost[i]
  }
  list(Beta = Beta, Tau = Tau, model.probs = ModPost, log.scores = LogScores, 
       Betas = Betas, Taus = Taus)
}


DropNARows <- function(Y, X) {
  q <- ncol(Y)
  na.counts <- apply(is.na(Y), 1, sum)
  to.keep <- which(na.counts != q)
  list(Y = Y[to.keep,], X = X[to.keep,])
}


MultStreamWithImputation4 <- function(Y, X, Y.val, X.val, U, d, V, Beta0, Tau0, 
                                      max.ite = 100) {
  n <- nrow(Y)
  In <- diag(n)
  val.mse <- rep(NA, max.ite + 1)
  current.out <- NULL
  current.Y <- NULL
  ## iteration 0 #########################
  ite <- 0
  current.mse <- mean(apply((Y.val - X.val %*% Beta0)^2, 2, mean, na.rm = TRUE))
  new.Y <- 
    ACE(Y, sig = In, delt = solve(Tau0), sigi = In, delti = Tau0, M = X %*% Beta0)[[1]]
  lambda.grid <- GetRidgeGrid(d, U, V, apply(new.Y, 1, mean))
  new.out <- MultStream3(new.Y, U, d, V, lambda.grid, a.tau = 1e-3, b.tau = 1e-3)
  Beta <- new.out$Beta
  Tau <- new.out$Tau
  new.mse <- mean(apply((Y.val - X.val %*% Beta)^2, 2, mean, na.rm = TRUE))
  val.mse[1] <- new.mse
  cat("ite ", 0, "\n")
  ########################################
  mse.diff <- 1
  while ( (mse.diff > 0) & (ite < max.ite) ) {
    current.out <- new.out
    current.mse <- new.mse
    current.Y <- new.Y
    new.Y <- ACE(Y, sig = In, delt = solve(Tau), sigi = In, delti = Tau, M = X %*% Beta)[[1]]     
    lambda.grid <- GetRidgeGrid(d, U, V, apply(new.Y, 1, mean))
    new.out <- MultStream3(new.Y, U, d, V, lambda.grid, a.tau = 1e-3, b.tau = 1e-3)
    Beta <- new.out$Beta
    Tau <- new.out$Tau
    new.mse <- mean(apply((Y.val - X.val %*% Beta)^2, 2, mean, na.rm = TRUE))
    mse.diff <- current.mse - new.mse
    ite <- ite + 1
    val.mse[ite + 1] <- new.mse
    cat("ite ", ite, "\n")
  }
  list(Y.imputed = current.Y, out = current.out, ite = ite, mse.diff = mse.diff,
       val.mse = val.mse)
}


#####################################################
## univariate stream functions
#####################################################

GetStreamModelScores2 <- function(y, U, d, V, lambdas, a.tau, b.tau) {
  n <- length(y)
  le <- length(lambdas)
  log.scores <- rep(NA, le)
  scores <- rep(NA, le)
  Betas <- matrix(NA, nrow(V), le)
  taus <- rep(NA, le)
  aux1 <- lgamma(a.tau + n/2) - lgamma(a.tau) - n * log(2 * a.tau * pi)/2
  aux2 <- -n * log(b.tau/a.tau)/2
  yty <- crossprod(y)
  Uty <- crossprod(U, y)
  ytU <- t(Uty)
  for (i in 1:le) {
    aux3 <- diag(1/(1 + lambdas[i]/d^2))
    aux3 <- aux3 %*% Uty
    theta <- aux3/d
    Betas[, i] <- V %*% theta
    aux3 <- ytU %*% aux3
    aux3 <- yty - aux3
    aux3 <- 1 + aux3/(2 * b.tau)
    taus[i] <- (2 * a.tau + n)/(2 * b.tau * aux3)
    aux3 <- -(a.tau + n/2) * log(aux3)
    aux4 <- -0.5 * sum(log(1 + (d^2)/lambdas[i]))
    log.scores[i] <- aux1 + aux2 + aux3 + aux4
    scores[i] <- exp(log.scores[i])
  }
  list(scores = scores, log.scores = log.scores, Betas = Betas, taus = taus)
}


Stream2 <- function(y, U, d, V, lambda.grid, a.tau, b.tau) {
  k <- length(lambda.grid)
  out <- GetStreamModelScores2(y, U, d, V, lambda.grid, a.tau, b.tau)
  LogScores <- out$log.scores
  ModPost <- ModelProbs(x = LogScores)
  beta <- out$Betas %*% ModPost
  tau <- sum(out$taus * ModPost)
  list(beta = beta[, 1], tau = tau, model.probs = ModPost, 
       log.scores = LogScores, Betas = out$Betas, taus = out$taus)
}



GetRidgeGrid <- function(d, U, V, y, beta.max = 0.001, epsilon = 1e-6,
                         K = 100, upper.bound = TRUE) {
  r <- length(d)
  aux <- V %*% t(U) %*% y
  if (upper.bound) {
    lambda.max <- max(d * abs(aux[1:r,1]))/beta.max
  }
  else {
    lambda.max <- max(d * abs(aux[1:r,1]))/beta.max - max(d^2)
  }
  exp(seq(log(lambda.max), log(epsilon * lambda.max), length.out = K))
}




ModelProbs <- function(x) {
  le <- length(x)
  model.probs <- rep(NA, le)
  for (i in 1:le) {
    model.probs[i] <- 1/sum(exp(x - x[i]))
  }
  model.probs
}



#####################################################
## Ridge functions
#####################################################

RidgeCoeffs <- function(U, d, V, y, lambda.grid) {
  rhs <- t(U) %*% y
  k <- length(lambda.grid)
  dx <- length(d)
  div <- d^2 + rep(lambda.grid, rep(dx, k))
  a <- drop(d * rhs)/div
  dim(a) <- c(dx, k)
  V %*% a  
}



SelectRidgeBestModel <- function(B, y.val, X.val, lambda.grid) {
  n.lambda <- length(lambda.grid)
  Y.val <- matrix(rep(y.val, n.lambda), length(y.val), n.lambda)
  Pred <- X.val %*% B
  mses <- apply((Y.val - Pred)^2, 2, mean)
  best <- which.min(mses)
  list(beta = B[, best], best.lambda = lambda.grid[best])
}


#####################################################
## utility functions
#####################################################

NormalTrans <- function(x) {
  n <- sum(!is.na(x))
  r <- rank(x, na.last = "keep")
  qnorm((r - 0.5)/n)
}



SplitData <- function(nsplits, n, seed = 123) {
  set.seed(seed)
  splitid <- sample(seq(nsplits), n, replace = TRUE)
  list(test.ids = which(splitid == 1),
       val.ids = which(splitid == 2),
       train.ids = which((splitid != 1) & (splitid != 2)))
}



CreateSigma <- function(rho, p) {
  aux1 <- matrix(rep(1:p, p), p, p)
  aux2 <- matrix(rep(1:p, each = p), p, p) 
  rho^abs(aux1 - aux2)
}



GetBlockSizes <- function(p, max.block, min.block) {
  max.block <- min(c(p, max.block))
  block.sizes <- sample(min.block:max.block, ceiling(p/min.block), 
                        replace = TRUE)
  cum.sizes <- cumsum(block.sizes)
  block.sizes <- block.sizes[which(cum.sizes <= p)]
  aux <- which.min(block.sizes)
  block.sizes[aux] <- block.sizes[aux] + p - sum(block.sizes)
  block.sizes 
}



SimulateMultData4 <- function(n, Beta, sig, rho, eff = 1, alpha = 1, max.block, min.block) {
  p <- nrow(Beta)
  q <- ncol(Beta)
  block.sizes <- GetBlockSizes(p, max.block, min.block)
  nblocks <- length(block.sizes)
  Sigma <- CreateSigma(rho, block.sizes[1])
  X <- mvrnorm(n, rep(0, block.sizes[1]), Sigma)
  if (nblocks > 1) {
    for (i in 2:nblocks) {
      Sigma <- CreateSigma(rho, block.sizes[i])
      X <- cbind(X, mvrnorm(n, rep(0, block.sizes[i]), Sigma))  
    }
  }
  z <- rnorm(n)
  Y <- matrix(NA, n, q)
  for (i in 1:q) {
    Y[, i] <- alpha * (X %*% Beta[, i]) + (1 - alpha) * eff * z + sig[i] * rnorm(n)
  }
  list(Y = Y, X = X, z = z, block.sizes = block.sizes)
}



GenerateMissing <- function(X, missing.prob = 0.3) {
  n <- nrow(X)
  p <- ncol(X)
  x <- as.vector(X)
  le <- n * p
  x[sample(seq(le), round(le * missing.prob), replace = FALSE)] <- NA
  matrix(x, n, p)
}


GetR2 <- function(Y, Pred) {
  Y.bar <- apply(Y, 2, mean, na.rm = TRUE)
  Y.bar <- matrix(rep(Y.bar, nrow(Y)), nrow(Y), ncol(Y))
  aux1 <- apply((Y - Pred)^2, 2, sum)
  aux2 <- apply((Y - Y.bar)^2, 2, sum)
  R2 <- 1 - aux1/aux2
  pmax(R2, rep(0, ncol(Y)))
}



FilterByCor <- function(feat.dat, resp.dat, cor.pval.thr = NULL) {
  featureCorrs <- apply(feat.dat, 2, "cor.test", resp.dat)
  featurePVals <- sapply(featureCorrs, FUN = function(x) x$p.value)
  which(featurePVals < cor.pval.thr)
}


FilterFeatures <- function(X, var.thr, split.index) {
  train.index <- split.index$train.ids
  val.index <- split.index$val.ids  
  test.index <- split.index$test.ids  
  Xtrain <- X[train.index,]
  Xval <- X[val.index,] 
  Xtest <- X[test.index,] 
  isNAtrain <- is.na(Xtrain)
  isNAtrain <- apply(isNAtrain, 2, sum)
  to.drop.train <- which(isNAtrain != 0)
  isNAval <- is.na(Xval)
  isNAval <- apply(isNAval, 2, sum)
  to.drop.val <- which(isNAval != 0)
  isNAtest <- is.na(Xtest)
  isNAtest <- apply(isNAtest, 2, sum)
  to.drop.test <- which(isNAtest != 0)
  to.drop.na <- unique(c(to.drop.train, to.drop.val, to.drop.test))
  X.var.train <- apply(Xtrain, 2, var, na.rm = TRUE)
  to.drop.train <- which(X.var.train <= var.thr)
  X.var.val <- apply(Xval, 2, var, na.rm = TRUE)
  to.drop.val <- which(X.var.val <= var.thr)  
  X.var.test <- apply(Xtest, 2, var, na.rm = TRUE)
  to.drop.test <- which(X.var.test <= var.thr)
  to.drop.var <- unique(c(to.drop.train, to.drop.val, to.drop.test))
  to.drop <- unique(c(to.drop.na, to.drop.var))
  sort(to.drop, decreasing = FALSE)    
}  



FilteringFeatures <- function(Y, X, var.thr, split.index, cor.pval = 0.001) {
  feat.to.drop <- FilterFeatures(X, var.thr, split.index)
  X <- X[, -feat.to.drop]
  Xtr <- X[split.index$train.ids, ]
  Ytr <- Y[split.index$train.ids,]
  n.drugs <- ncol(Y)
  to.keep <- vector(mode = "list", length = n.drugs)
  for (i in 1:n.drugs) {
    cat("drug = ", i, "\n")
    to.keep[[i]] <- FilterByCor(Xtr, Ytr[, i], cor.pval)
  }
  feat.to.keep <- sort(unique(unlist(to.keep)))
  list(Ytrain = Ytr, 
       Xtrain = Xtr[, feat.to.keep], 
       Yval = Y[split.index$val.ids,],
       Xval = X[split.index$val.ids, feat.to.keep],       
       Ytest = Y[split.index$test.ids,],
       Xtest = X[split.index$test.ids, feat.to.keep])
}



RemoveNAs <- function(y, X) {
  to.drop <- is.na(y)
  if(sum(to.drop) > 0) {
    out <-list(y = y[!to.drop], X = X[!to.drop,])
  }
  else {
    out <- list(y = y, X = X)
  }
}



RunUnivariateAnalysis <- function(Ytrain, Xtrain, Yval, Xval, Ytest, Xtest) {
  n.drugs <- ncol(Ytrain)
  p <- ncol(Xtrain)
  Rbetas <- matrix(NA, p, n.drugs)
  Sbetas <- matrix(NA, p, n.drugs)
  Staus <- rep(NA, n.drugs)
  Smse <- rep(NA, n.drugs)
  Rmse <- rep(NA, n.drugs)
  n.train <- rep(NA, n.drugs)
  n.val <- rep(NA, n.drugs)
  n.test <- rep(NA, n.drugs)
  for (i in 1:n.drugs) {
    cat("drug", i, "\n")
    ftrain <- RemoveNAs(y = Ytrain[, i], X = Xtrain)
    fval <- RemoveNAs(y = Yval[, i], X = Xval)
    ftest <- RemoveNAs(y = Ytest[, i], X = Xtest)
    ytrain <- scale(ftrain$y)[, 1]
    yval <- scale(fval$y)[, 1]
    ytest <- scale(ftest$y)[, 1]
    xtrain <- scale(ftrain$X)
    xval <- scale(fval$X)
    xtest <- scale(ftest$X)
    n.train[i] <- length(ytrain)
    n.val[i] <- length(yval)
    n.test[i] <- length(ytest)
    fsvd <- fast.svd(xtrain)
    lambda.grid <- GetRidgeGrid(fsvd$d, fsvd$u, fsvd$v, ytrain) 
    Betas <- RidgeCoeffs(fsvd$u, fsvd$d, fsvd$v, ytrain, lambda.grid)
    Rbetas[, i] <- SelectRidgeBestModel(Betas, yval, xval, lambda.grid)[[1]] 
    out <- Stream2(ytrain, fsvd$u, fsvd$d, fsvd$v, lambda.grid, a.tau = 1e-3, b.tau = 1e-3)   
    Sbetas[, i] <- out$beta
    Staus[i] <- out$tau 
    Smse[i] <- mean((ytest - xtest %*% out$beta)^2)
    Rmse[i] <- mean((ytest - xtest %*% Rbetas[, i])^2)
  }
  list(Smse = Smse, Sbetas = Sbetas, Staus = Staus, n.train = n.train, n.val = n.val,
       n.test = n.test, Rmse = Rmse, Rbetas = Rbetas)
}



RunSimulations1 <- function(mu.grid, sig.grid, my.seed, p, q, n, Sig, rho, eff, alpha,
                            missing.prob = 0.25, missing.seed = 12345) {
  n.mu <- length(mu.grid)
  n.sig <- length(sig.grid)
  n.sim <- n.mu * n.sig
  mu.nms <- as.character(mu.grid)
  sig.nms <- as.character(sig.grid)
  Y.cor <- matrix(NA, n.sig, n.mu)
  Pred.cor <- matrix(NA, n.sig, n.mu)
  Pred.R2 <- matrix(NA, n.sig, n.mu)
  Z.R2 <- matrix(NA, n.sig, n.mu)
  dimnames(Y.cor) <- list(sig.nms, mu.nms)
  dimnames(Pred.cor) <- list(sig.nms, mu.nms)
  dimnames(Pred.R2) <- list(sig.nms, mu.nms)
  dimnames(Z.R2) <- list(sig.nms, mu.nms)
  nite <- rep(NA, n.sim)
  BMSE.ridge <- matrix(NA, n.sim, q) ## baseline
  BMSE.stream <- matrix(NA, n.sim, q) ## baseline
  RMSE <- matrix(NA, n.sim, q) ## Ridge
  UMSE <- matrix(NA, n.sim, q) ## Stream
  MMSE <- matrix(NA, n.sim, q) ## IMV-Stream
  IMP.ridge <- matrix(NA, n.sig, n.mu)
  IMP.stream <- matrix(NA, n.sig, n.mu)
  dimnames(IMP.ridge) <- list(sig.nms, mu.nms)
  dimnames(IMP.stream) <- list(sig.nms, mu.nms)
  sim <- 1
  for (i in seq(n.sig)) {
    for (j in seq(n.mu)) {
      set.seed(my.seed)
      beta <- t(mvrnorm(q, rep(mu.grid[j], p), Sig))
      sig <- rep(sig.grid[i], q)
      dat1 <- SimulateMultData4(n, beta, sig, rho, eff = eff, alpha = alpha, 
                                max.block = 200, min.block = 20)
      Ytrain <- scale(dat1$Y)
      Xtrain <- scale(dat1$X)
      scaled.beta <- beta/apply(dat1$Y, 2, sd)
      scaled.eff <- rep(eff, q)/apply(dat1$Y, 2, sd)
      scaled.eff <- matrix(rep(scaled.eff, n), n, q)
      Pred.R2[i, j] <-  mean(GetR2(Ytrain, Xtrain %*% scaled.beta))
      Z.R2[i, j] <- mean(GetR2(Ytrain, scaled.eff * dat1$z))    
      aux <- cor(Ytrain)
      aux <- aux[lower.tri(aux)]
      Y.cor[i, j] <- mean(aux)
      aux <- cor(Ytrain, Xtrain %*% scaled.beta)
      aux <- diag(aux)
      Pred.cor[i, j] <- mean(aux)
      ###################################
      dat2 <- SimulateMultData4(n, beta, sig, rho, eff = eff, alpha = alpha, 
                                max.block = 200, min.block = 20)
      Yval <- scale(dat2$Y)
      Xval <- scale(dat2$X) 
      dat3 <- SimulateMultData4(n, beta, sig, rho, eff = eff, alpha = alpha, 
                                max.block = 200, min.block = 20)
      Ytest <- scale(dat3$Y)
      Xtest <- scale(dat3$X)
      cat("full data univariate anaylis for sim", sim, "\n")
      uni <- RunUnivariateAnalysis(Ytrain, Xtrain, Yval, Xval, Ytest, Xtest)
      BMSE.ridge[sim,] <- uni$Rmse
      BMSE.stream[sim,] <- uni$Smse
      cat("missing data univariate anaylis for sim", sim, "\n")
      set.seed(missing.seed)
      Ytrain.m <- GenerateMissing(Ytrain, missing.prob = missing.prob)
      uni.m <- RunUnivariateAnalysis(Ytrain.m, Xtrain, Yval, Xval, Ytest, Xtest)
      RMSE[sim,] <- uni.m$Rmse
      UMSE[sim,] <- uni.m$Smse
      Beta0 <- uni.m$Sbetas
      Tau0 <- diag(uni.m$Staus)
      cat("IMV-Stream anaylis for sim", sim, "\n")
      fsvd2 <- fast.svd(Xtrain)
      mout.m <- MultStreamWithImputation4(Ytrain.m, Xtrain, Yval, Xval, fsvd2$u, 
                                          fsvd2$d, fsvd2$v, Beta0, Tau0, max.ite = 1000)
      nite[sim] <- mout.m$ite
      MMSE[sim,] <- apply((Ytest - Xtest %*% mout.m$out$Beta)^2, 2, mean)
      IMP.ridge[i, j] <- mean( (RMSE[sim,] - MMSE[sim,])/RMSE[sim,] )
      IMP.stream[i, j] <- mean( (UMSE[sim,] - MMSE[sim,])/UMSE[sim,] )
      ####################   
      sim <- sim + 1
    }
  }
  list(Y.cor = Y.cor, Pred.cor = Pred.cor, Pred.R2 = Pred.R2, Z.R2 = Z.R2,
       IMP.ridge = IMP.ridge, IMP.stream = IMP.stream, RMSE = RMSE, UMSE = UMSE,
       MMSE = MMSE, nite = nite)
} 



RunSimulations2 <- function(mu.grid, alpha.grid, my.seed, p, q, n, Sig, rho, eff, sig,
                            missing.prob = 0.25, missing.seed = 12345) {
  n.mu <- length(mu.grid)
  n.alpha <- length(alpha.grid)
  n.sim <- n.mu * n.alpha
  mu.nms <- as.character(mu.grid)
  alpha.nms <- as.character(alpha.grid)
  Y.cor <- matrix(NA, n.alpha, n.mu)
  Pred.cor <- matrix(NA, n.alpha, n.mu)
  Pred.R2 <- matrix(NA, n.alpha, n.mu)
  Z.R2 <- matrix(NA, n.alpha, n.mu)
  dimnames(Y.cor) <- list(alpha.nms, mu.nms)
  dimnames(Pred.cor) <- list(alpha.nms, mu.nms)
  dimnames(Pred.R2) <- list(alpha.nms, mu.nms)
  dimnames(Z.R2) <- list(alpha.nms, mu.nms)
  nite <- rep(NA, n.sim)
  BMSE.ridge <- matrix(NA, n.sim, q) ## baseline
  BMSE.stream <- matrix(NA, n.sim, q) ## baseline
  RMSE <- matrix(NA, n.sim, q) ## Ridge
  UMSE <- matrix(NA, n.sim, q) ## Stream
  MMSE <- matrix(NA, n.sim, q) ## IMV-Stream
  IMP.ridge <- matrix(NA, n.alpha, n.mu)
  IMP.stream <- matrix(NA, n.alpha, n.mu)
  dimnames(IMP.ridge) <- list(alpha.nms, mu.nms)
  dimnames(IMP.stream) <- list(alpha.nms, mu.nms)
  sig <- rep(sig, q)
  sim <- 1
  for (i in seq(n.alpha)) {
    for (j in seq(n.mu)) {
      set.seed(my.seed)
      beta <- t(mvrnorm(q, rep(mu.grid[j], p), Sig))
      dat1 <- SimulateMultData4(n, beta, sig, rho, eff = eff, alpha = alpha.grid[i], 
                                max.block = 200, min.block = 20)
      Ytrain <- scale(dat1$Y)
      Xtrain <- scale(dat1$X)
      scaled.beta <- beta/apply(dat1$Y, 2, sd)
      scaled.eff <- rep(eff, q)/apply(dat1$Y, 2, sd)
      scaled.eff <- matrix(rep(scaled.eff, n), n, q)
      Pred.R2[i, j] <-  mean(GetR2(Ytrain, alpha.grid[i] * Xtrain %*% scaled.beta))
      Z.R2[i, j] <- mean(GetR2(Ytrain, (1 - alpha.grid[i]) * scaled.eff * dat1$z))    
      aux <- cor(Ytrain)
      aux <- aux[lower.tri(aux)]
      Y.cor[i, j] <- mean(aux)
      aux <- cor(Ytrain, alpha.grid[i] * Xtrain %*% scaled.beta)
      aux <- diag(aux)
      Pred.cor[i, j] <- mean(aux)
      ###################################
      dat2 <- SimulateMultData4(n, beta, sig, rho, eff = eff, alpha = alpha.grid[i], 
                                max.block = 200, min.block = 20)
      Yval <- scale(dat2$Y)
      Xval <- scale(dat2$X) 
      dat3 <- SimulateMultData4(n, beta, sig, rho, eff = eff, alpha = alpha.grid[i], 
                                max.block = 200, min.block = 20)
      Ytest <- scale(dat3$Y)
      Xtest <- scale(dat3$X)
      cat("full data univariate anaylis for sim", sim, "\n")
      uni <- RunUnivariateAnalysis(Ytrain, Xtrain, Yval, Xval, Ytest, Xtest)
      BMSE.ridge[sim,] <- uni$Rmse
      BMSE.stream[sim,] <- uni$Smse
      cat("missing data univariate anaylis for sim", sim, "\n")
      set.seed(missing.seed)
      Ytrain.m <- GenerateMissing(Ytrain, missing.prob = missing.prob)
      uni.m <- RunUnivariateAnalysis(Ytrain.m, Xtrain, Yval, Xval, Ytest, Xtest)
      RMSE[sim,] <- uni.m$Rmse
      UMSE[sim,] <- uni.m$Smse
      Beta0 <- uni.m$Sbetas
      Tau0 <- diag(uni.m$Staus)
      cat("IMV-Stream anaylis for sim", sim, "\n")
      fsvd2 <- fast.svd(Xtrain)
      mout.m <- MultStreamWithImputation4(Ytrain.m, Xtrain, Yval, Xval, fsvd2$u, 
                                          fsvd2$d, fsvd2$v, Beta0, Tau0, max.ite = 1000)
      nite[sim] <- mout.m$ite
      MMSE[sim,] <- apply((Ytest - Xtest %*% mout.m$out$Beta)^2, 2, mean)
      IMP.ridge[i, j] <- mean( (RMSE[sim,] - MMSE[sim,])/RMSE[sim,] )
      IMP.stream[i, j] <- mean( (UMSE[sim,] - MMSE[sim,])/UMSE[sim,] )
      ####################   
      sim <- sim + 1
    }
  }
  list(Y.cor = Y.cor, Pred.cor = Pred.cor, Pred.R2 = Pred.R2, Z.R2 = Z.R2,
       IMP.ridge = IMP.ridge, IMP.stream = IMP.stream, RMSE = RMSE, UMSE = UMSE,
       MMSE = MMSE, nite = nite)
} 


#####################################################################
## utility functions for the DOSE simulation
#####################################################################

TransformDesign <- function(x, par.ranges, is.discrete, lin.trans) {
  X <- x
  for (i in 1:ncol(x)) {
    par.range <- par.ranges[[i]]
    par.min <- min(par.range)
    par.max <- max(par.range)
    if (lin.trans[i]) {
      aux <- x[, i] * (par.max - par.min) + par.min
    }
    else {
      par.mid <- 1
      ind1 <- x[, i] <= 0.5
      ind2 <- x[, i] > 0.5
      aux <- rep(NA, length(x[, i]))
      aux[ind1] <- 2 * x[ind1, i] * (par.mid - par.min) + par.min
      aux[ind2] <- 
        2 * x[ind2, i] * (par.max - par.mid) + 2 * par.mid - par.max 
    }
    if (is.discrete[[i]]) {
      X[, i] <- round(aux)
    }
    else {
      X[, i] <- aux
    }  
  }
  X
}



RunSimulationsDOSE <- function(D, missing.prob = 0.25, data.generation.seeds, 
                               missing.seed = 12345) {
  n.sim <- nrow(D)
  out <- matrix(NA, n.sim, 16)
  colnames(out) <- c(colnames(D), "Y.cor", "Pred.cor", "Pred.R2", "Z.R2", 
                     "IMP.ridge", "IMP.stream", "nite")
  BMSE.ridge <- vector(mode = "list", length = n.sim) ## baseline
  BMSE.stream <- vector(mode = "list", length = n.sim) ## baseline
  RMSE <- vector(mode = "list", length = n.sim) ## Ridge
  UMSE <- vector(mode = "list", length = n.sim) ## Stream
  MMSE <- vector(mode = "list", length = n.sim) ## IMV-Stream
  for (i in seq(n.sim)) {
    cat("sim = ", i, "\n")
    out[i, 1:9] <- D[i,]
    n <- D[i, "n"]
    p <- D[i, "p"]
    q <- D[i, "q"]
    rho.X <- D[i, "rho.X"]
    rho.beta <- D[i, "rho.beta"]
    mu.beta <- D[i, "mu.beta"]
    Sig.beta <- CreateSigma(rho.beta, p)
    sigma2 <- D[i, "sigma2"]
    alpha <- D[i, "alpha"]
    gamm <- D[i, "gamma"]
    set.seed(data.generation.seeds[i])
    beta <- t(mvrnorm(q, rep(mu.beta, p), Sig.beta))
    sig <- rep(sigma2, q)
    dat1 <- SimulateMultData4(n, beta, sig, rho.X, eff = gamm, alpha = alpha, 
                              max.block = 200, min.block = 20)
    Ytrain <- scale(dat1$Y)
    Xtrain <- scale(dat1$X)
    scaled.beta <- beta/apply(dat1$Y, 2, sd)
    scaled.eff <- rep(gamm, q)/apply(dat1$Y, 2, sd)
    scaled.eff <- matrix(rep(scaled.eff, n), n, q)
    out[i, 12] <-  mean(GetR2(Ytrain, alpha * Xtrain %*% scaled.beta))
    out[i, 13] <- mean(GetR2(Ytrain, (1 - alpha) * scaled.eff * dat1$z))    
    aux <- cor(Ytrain)
    aux <- aux[lower.tri(aux)]
    out[i, 10] <- mean(aux)
    aux <- cor(Ytrain, alpha * Xtrain %*% scaled.beta)
    aux <- diag(aux)
    out[i, 11] <- mean(aux) 
    #################################
    dat2 <- SimulateMultData4(n, beta, sig, rho.X, eff = gamm, alpha = alpha, 
                              max.block = 200, min.block = 20)
    Yval <- scale(dat2$Y)
    Xval <- scale(dat2$X) 
    dat3 <- SimulateMultData4(n, beta, sig, rho.X, eff = gamm, alpha = alpha, 
                              max.block = 200, min.block = 20)
    Ytest <- scale(dat3$Y)
    Xtest <- scale(dat3$X)
    cat("full data univariate anaylis for sim", i, "\n")
    uni <- RunUnivariateAnalysis(Ytrain, Xtrain, Yval, Xval, Ytest, Xtest)
    BMSE.ridge[[i]] <- uni$Rmse
    BMSE.stream[[i]] <- uni$Smse
    cat("missing data univariate anaylis for sim", i, "\n")
    set.seed(missing.seed)
    Ytrain.m <- GenerateMissing(Ytrain, missing.prob = missing.prob)
    uni.m <- RunUnivariateAnalysis(Ytrain.m, Xtrain, Yval, Xval, Ytest, Xtest)
    RMSE[[i]] <- uni.m$Rmse
    UMSE[[i]] <- uni.m$Smse
    Beta0 <- uni.m$Sbetas
    Tau0 <- diag(uni.m$Staus)
    cat("IMV-Stream anaylis for sim", i, "\n")
    fsvd2 <- fast.svd(Xtrain)
    mout.m <- MultStreamWithImputation4(Ytrain.m, Xtrain, Yval, Xval, fsvd2$u, 
                                        fsvd2$d, fsvd2$v, Beta0, Tau0, max.ite = 1000)
    out[i, 16] <- mout.m$ite
    MMSE[[i]] <- apply((Ytest - Xtest %*% mout.m$out$Beta)^2, 2, mean)
    out[i, 14] <- mean( (RMSE[[i]] - MMSE[[i]])/RMSE[[i]] )
    out[i, 15] <- mean( (UMSE[[i]] - MMSE[[i]])/UMSE[[i]] )
  }
  list(out = out, RMSE = RMSE, UMSE = UMSE, MMSE = MMSE)
}  



################################################################
## code for figures
################################################################

GetDeltaMSE1 <- function(x1, x2) {
  n <- nrow(x1)
  delta.mse <- rep(NA, n)
  for (i in seq(n)) {
    delta.mse[i] <- mean(x1[i,] - x2[i,])
  }
  delta.mse
}


GetDeltaMSE2 <- function(x) {
  n <- length(x$RMSE)
  delta.mse <- rep(NA, n)
  for (i in seq(n)) {
    delta.mse[i] <- mean(x$RMSE[[i]] - x$MMSE[[i]])
  }
  delta.mse
}


MyLegend <- function(imp, x1, x2, y0, ydist, lower.bound, upper.bound, col.low = "red", 
                     col.mid = "green", col.high = "blue", diff = 1e-4, cex1 = 5, cex2 = 2) {
  n <- length(imp)
  legendcolors <- GetImpColor(imp, lower.bound, upper.bound, col.low, col.mid, col.high, diff)
  y <- y0
  for (i in seq(n)) {
    points(x1, y, bg = legendcolors[i], pch = 22, cex = cex1)
    text(x2, y, format(round(imp[i], 2), nsmall = 2), cex = cex2) 
    y <- y - ydist
  } 
}



GetImpColor <- function(imp, lower.bound, upper.bound, col.low = "red", 
                        col.mid = "green", col.high = "blue", diff = 1e-4) {
  col.bins <- seq(lower.bound, upper.bound, by = diff)
  ncol <- length(col.bins)
  mycol <- colorpanel(n = ncol, low = col.low, mid = col.mid, high = col.high)
  n <- length(imp)
  pos <- rep(NA, n)
  for (i in seq(n)) {
    pos[i] <- which.min(abs(col.bins - imp[i]))
  }
  mycol[pos]
}


GetMatrixForImage <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  xx <- matrix(NA, nc, nr)
  rownames(xx) <- colnames(x)
  colnames(xx) <- rev(rownames(x))
  for (i in seq(nc)) {
    xx[i, ] <- rev(x[, i])
  }
  xx
}


PlotCorImage2 <- function(resp, lin.pred, z, ncolor = 100, mycol = redblue(ncolor), 
                          mycex = 1, myline = 0.2, nskip = 1, main = NULL) {
  q <- ncol(resp)
  c1 <- cor(lin.pred)
  c2 <- cor(resp)
  l1 <- lower.tri(c1)
  u2 <- upper.tri(c2)
  m <- matrix(1, q, q)
  m[l1] <- c1[l1]
  m[u2] <- c2[u2]
  c3 <- cor(resp, z)[, 1]
  c4 <- diag(cor(resp, lin.pred))
  nas1 <- nas2 <- matrix(NA, nskip, q)
  m <- rbind(c3, nas1, m, nas2, c4)
  m <- GetMatrixForImage(m)
  mybreaks <- seq(-1, 1, length.out = ncolor + 1)
  image(m, col = mycol, xaxt = "n", yaxt = "n", breaks = mybreaks, bty = "n", 
        main = main)
  mtext("response", side = 3, cex = mycex, line = myline)
  mtext("response", side = 4, cex = mycex, line = myline)
  mtext(c("z", "predictor"), side = 2, cex = mycex, line = myline,
        at = c(1, 0.5))
  mtext("predictor", side = 1, cex = mycex, line = myline)
}


PlotMSE2 <- function(x, pchs = 1:3, cols = 1:3, cex = 1, line = 0,
                     drug.cex = 1, legend.pos = "bottomright", xlab = "",
                     cex.lab = 1, cex.axis = 1, cex.legend = 1, xlim = NULL,
                     bottom.mar = 5, left.mar = 4, top.mar = 4, right.mar = 2) {
  mod.nms <- rownames(x) 
  drug.nms <- colnames(x)
  ndrugs <- length(drug.nms)
  bars <- seq(1, 0, length.out = ndrugs)
  par(mar = c(bottom.mar, left.mar, top.mar, right.mar) + 0.1, mgp = c(2, 1, 0))
  plot(x[1,], 1:ndrugs, type = "n", xlim = xlim, 
       yaxt = "n", ylab = "", xlab = xlab, cex.lab = cex.lab, 
       cex.axis = cex.axis, ylim = c(0 - 1/ndrugs, 1 + 1/ndrugs))
  for (i in 1:ndrugs) {
    abline(h = bars[i], lty = 3, col = "black")
  }
  mtext(drug.nms, side = 4, at = bars, las = 2, cex = drug.cex, line = line)
  points(x[1,], bars, pch = pchs[1], col = cols[1], cex = cex)
  points(x[2,], bars, pch = pchs[2], col = cols[2], cex = cex)
  points(x[3,], bars, pch = pchs[3], col = cols[3], cex = cex)
  legend(legend.pos, legend = mod.nms, text.col = cols, cex = cex.legend, 
         bty = "n")
  par(mar = c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))
}


MissingDataPlot <- function(x, cex.row = 1, cex.col = 1, bottom.mar = 5, 
                            left.mar = 4, top.mar = 4, right.mar = 2, 
                            line = 0) {
  GetMatrixForImage <- function(x) {
    nr <- nrow(x)
    nc <- ncol(x)
    xx <- matrix(NA, nc, nr)
    rownames(xx) <- colnames(x)
    colnames(xx) <- rev(rownames(x))
    for (i in seq(nc)) {
      xx[i, ] <- rev(x[, i])
    }
    xx
  }
  par(mar = c(bottom.mar, left.mar, top.mar, right.mar) + 0.1)
  x <- GetMatrixForImage(x)
  y <- as.vector(x)
  yy <- y
  yy[!is.na(y)] <- 1
  yy[is.na(y)] <- 0
  xx <- matrix(yy, nrow(x), ncol(x))
  dimnames(xx) <- dimnames(x)
  image(xx, col = redblue(10), xaxt = "n", yaxt = "n")
  nr <- nrow(xx)
  nc <- ncol(xx)
  mtext(rownames(xx), side = 1, at = seq(0, 1, length.out = nrow(xx)), las = 2, 
        cex = cex.col, line = line)
  mtext(colnames(xx), side = 4, at = seq(0, 1, length.out = ncol(xx)), las = 2, 
        cex = cex.row, line = line)
  par(mar = c(5, 4, 4, 2) + 0.1)
}

