#=====================================================================================
SETUP = function(X, Y, lambda1, lambda2){
  sm = list()
  sm$X = X
  sm$Y = Y
  sm$dimension = ncol(X[[1]])
  sm$N = lapply(Y, length)
  sm$num_study = length(X)
  sm$lambda1 = lambda1
  sm$lambda2 = lambda2
  ## Fit an model to generate initial coefficients
  linsvc <- foreach(m = 1:sm$num_study) %do% {
    svm(X[[m]][,-1], factor(Y[[m]]), kernel = "linear", class.weights = NULL, fitted = FALSE)
  }
  sm$b <- foreach(m = 1:sm$num_study) %do% {
    c(linsvc[[m]]$rho,
      apply(as.vector(linsvc[[m]]$coefs) * linsvc[[m]]$SV, 2, sum))
  }
  sm$Z  <- foreach(m = 1:sm$num_study) %do% {
    Y[[m]] * X[[m]]
  }
  sm$W  <- foreach(m = 1:sm$num_study) %do% {
  diag(as.vector(1/abs(Y[[m]] - X[[m]] %*% sm$b[[m]])))
  }
  sm$One <- foreach(m = 1:sm$num_study) %do% {
  rep(1, sm$N[[m]])
  }
  s_betac_all <- sum(foreach(j = 2:sm$dimension,.combine=c) %do% {
  sqrt(sum(foreach(mm = 1:sm$num_study, .combine=c) %do% {sm$b[[mm]][j]}^2))
  })
  sm$Q_lambda <- do.call(sum, foreach(m = 1:sm$num_study) %do% {
        -1/(2*sm$N[[m]]) * sm$One[[m]] %*% sm$Z[[m]] %*% sm$b[[m]] +
        1/(4*sm$N[[m]]) * (t(Y[[m]]) %*% sm$W[[m]] %*% Y[[m]] -
        2 * t(sm$b[[m]]) %*% t(X[[m]]) %*% sm$W[[m]] %*% Y[[m]] +
        t(sm$b[[m]]) %*% t(X[[m]]) %*% sm$W[[m]] %*% X[[m]] %*% sm$b[[m]]) +
        sm$lambda1 * s_betac_all +
        sm$lambda2 * sum(abs(sm$b[[m]][-1]))
        })
  # Miscellaneous constants
  sm$scale = 1e-5
  sm$max_iteration = 50
  sm$max_iterations_tot = 5
  sm$max_iterations = 20
  sm$epsilon = 1e-5
  sm$epsilon_total = 1e-5
  sm$c = 1e-8
  return(sm)
}
#=====================================================================================
iSVM = function(sm, is.constant=TRUE){
  for(iter in 1:sm$max_iterations_tot)
  {
    ## Objective function value for coordinate descent algorithm
    old_Q_lambda = sm$Q_lambda
    ## Objective function value for Newton-Raphson algorithm
    q = sm$Q_lambda
    # Iteratively update
    tilde_coefficients_old = sm$b
    tilde_coefficients = tilde_coefficients_old
  for (m in 1:sm$num_study)
  {
    for (j in 1:sm$dimension)
    {
      cat(paste("\n", iter, "th Iteration :", "Estimating", j, "th variable \n"))
      cat(paste("Newton-Raphson iteration "))
      for(newton_iter in 1:sm$max_iteration)
      {
        cat(paste(newton_iter))
        cat(paste("est:",tilde_coefficients[[m]][j]),"\n")
        q_old = q
        sm$W <- foreach(mm = 1:sm$num_study) %do% {
          diag(as.vector(1/abs(sm$Y[[mm]] - sm$X[[mm]] %*% tilde_coefficients_old[[mm]])))
        }
        tilde_coefficients_old = tilde_coefficients
        ## L2-penalty
        group_beta <- foreach(mm = 1:sm$num_study, .combine=c) %do% {tilde_coefficients[[mm]][j]}
        s_betac = sqrt(sum(group_beta^2) + sm$c)
        if(j == 1 & is.constant){ ## We do not penalize for constant.
          d <- (- 1/(2*sm$N[[m]]) * ( t(sm$X[[m]]) %*% sm$W[[m]]  %*% (sm$Y[[m]] - sm$X[[m]] %*% tilde_coefficients[[m]]) +
                                       t(sm$Z[[m]]) %*% sm$One[[m]] ))[j]
        } else {
          d <- (- 1/(2*sm$N[[m]]) * ( t(sm$X[[m]]) %*% sm$W[[m]]  %*% (sm$Y[[m]] - sm$X[[m]] %*% tilde_coefficients[[m]]) +
                                        t(sm$Z[[m]]) %*% sm$One[[m]] ))[j] +
            sm$lambda1 * tilde_coefficients[[m]][j] / s_betac +
            sm$lambda2 * sign(tilde_coefficients[[m]][j])
        }
        ## Second derivative function
        #tmp_power <- (tilde_coefficients[[m]][j] * sign(exp(tilde_coefficients[[m]][j]))) / sm$scale
        tmp_power <- (tilde_coefficients[[m]][j]) / sm$scale
        if(exp(tmp_power)==Inf){
          cdf = 1
        } else{
          # cdf = 2*( 1 / (1 + exp(-tmp_power / sm$scale))) - 1
          cdf = 2*(exp(tmp_power) / (1+exp(tmp_power))) - 1
        }
        w =  sm$scale * cdf * (1 - cdf)
        if(j == 1 & is.constant){ ## We do not penalize for constant.
        d2 =  (1/(2*sm$N[[m]]) * (t(sm$X[[m]]) %*% sm$W[[m]] %*% sm$X[[m]]))[j, j]
        } else {
        d2 =  (1/(2*sm$N[[m]]) * (t(sm$X[[m]]) %*% sm$W[[m]] %*% sm$X[[m]]))[j, j] +
              sm$lambda1 * (1/sqrt(s_betac) - tilde_coefficients[[m]][j]*tilde_coefficients[[m]][j] / (s_betac)^(1.5)) +
              sm$lambda2 * w
        }
        diff = d / d2
        # tilde_coefficients[j] = tilde_coefficients_old[j] - diff
        for (iters in 1:sm$max_iterations)
        {
          cat(paste("diff:",diff, "estimation :", tilde_coefficients[[m]][j], " \n"))
          tilde_coefficients[[m]][j] = tilde_coefficients_old[[m]][j] - diff
          s_betac_all <- sum(foreach(jj = 2:sm$dimension,.combine=c) %do% {
            sqrt(sum(foreach(mm = 1:sm$num_study, .combine=c) %do% {tilde_coefficients[[mm]][jj]}^2))
          })
          q = do.call(sum, foreach(mm = 1:sm$num_study) %do% {
              -1/(2*sm$N[[mm]]) * sm$One[[mm]] %*% sm$Z[[mm]] %*% tilde_coefficients[[mm]] +
              1/(4*sm$N[[mm]]) * (t(sm$Y[[mm]]) %*% sm$W[[mm]] %*% sm$Y[[mm]] -
              2 * t(tilde_coefficients[[mm]]) %*% t(sm$X[[mm]]) %*% sm$W[[mm]] %*% sm$Y[[mm]] +
              t(tilde_coefficients[[mm]]) %*% t(sm$X[[mm]]) %*% sm$W[[mm]] %*% sm$X[[mm]] %*% tilde_coefficients[[mm]]) +
              sm$lambda1 * s_betac_all +
              sm$lambda2 * sum(abs(tilde_coefficients[[mm]][-1]))
          })
          # compute Q_lambda
          # q = Q + Q_penalty(tilde_coefficients)
          if(q <= q_old) break
          diff = 0.5 * diff
        }
        if(sum(foreach(mm = 1:sm$num_study, .combine=c) %do% {
          abs(tilde_coefficients[[mm]][j] - tilde_coefficients_old[[mm]][j]) }) < sm$epsilon) break
      }
    }
  }
  sm$b <- tilde_coefficients
  Q = do.call(sum, foreach(mm = 1:sm$num_study) %do% {
    -1/(2*sm$N[[mm]]) * sm$One[[mm]] %*% sm$Z[[mm]] %*% tilde_coefficients[[mm]] +
      1/(4*sm$N[[mm]]) * (t(sm$Y[[mm]]) %*% sm$W[[mm]] %*% sm$Y[[mm]] -
                            2 * t(tilde_coefficients[[mm]]) %*% t(sm$X[[mm]]) %*% sm$W[[mm]] %*% sm$Y[[mm]] +
                            t(tilde_coefficients[[mm]]) %*% t(sm$X[[mm]]) %*% sm$W[[mm]] %*% sm$X[[mm]] %*% tilde_coefficients[[mm]]) +
      sm$lambda1 * s_betac_all +
      sm$lambda2 * sum(abs(tilde_coefficients[[mm]]))
  })
  sm$Q_lambda = Q
  ## sm$Q_lambda = Q + Q_penalty(sm$b)
  ## Total value of the objective function
  if(abs(sm$Q_lambda - old_Q_lambda) < sm$epsilon_total)
    break
  }
  return(sm)
}
#=====================================================================================
## The Value of penalty terms with coefficients
Q_penalty <- function(sm, Coef){
  tmp_group_sq_sum <-c()
  for(jj in 1:length(sm$beta_group)){
    tmp_group_sq_sum <- c(tmp_group_sq_sum, sum(Coef[sm$beta_group[[jj]]$idx]^2))
  }
  penalty1 = sum(sqrt(tmp_group_sq_sum))  ## Group lasso
  penalty2 = sum(abs(Coef))     ## L1 lasso
  sm$lambda1 * penalty1 + sm$lambda2 * penalty2
}
#=====================================================================================
f_K_fold <- function(Nobs,K){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d)
    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}
#=====================================================================================
simu.data <- function(K = 10, N = 50, G = 10000, g = round(G * 0.1), p = rep(1/K, K), mu = runif(G, 0.1, 1) * sample(c(-1, 1), G, replace = TRUE), sigma = 1) {
  k <- c(sort(sample(1:K, g, replace = TRUE), decreasing = TRUE), rep(0, G - g))
  truth <- matrix(FALSE, G, K)
  for(i in 1:G) truth[i, sample(1:K, k[i])] <- TRUE
  result <- list()
  for(i in 1:K) {
    data0 <- matrix(rnorm(G * N, 0, sigma), G, N)
    mu.i <- mu
    mu.i[!truth[, i]] <- 0
    data1 <- matrix(rnorm(G * N, rep(mu.i, N), sigma), G, N)
    result[[i]] <- cbind(data0, data1)
  }
  return(list(data = result, truth = truth))
}
#=====================================================================================
sample.correlated.data <- function(p, n, rho, df, sigma) {
  s <- matrix(rho, p, p)
  s <- s + diag(rep(1 - rho, p))
  cov.mat <- riwish(df, s)
  cor.mat <- cov2cor(cov.mat)
  eig <- eigen(cor.mat)
  L <- eig[[2]] %*% diag(sqrt(eig[[1]])) %*% t(eig[[2]])
  x <- matrix(rnorm(p * n, 0, sigma), p, n)
  xx <- L %*% x
  return(xx)
}
#=====================================================================================
simu.data.correlated <- function(K = 10, N = 50, G = 10000, g = round(G * 0.1), p = rep(1/K, K),
                                 mu = matrix(runif(G * K, 0.5, 1) * sample(c(-1, 1), G * K, replace = TRUE), G, K), sigma = 1,
                                 rho = 0.5, clust.size = 20, n.clust = 200, rho.prior = 0.5, df.prior = clust.size * 3) {
  k <- c(sort(sample(1:K, g, replace = TRUE), decreasing = TRUE), rep(0, G - g))
  truth <- matrix(FALSE, G, K)
  for(i in 1:G) truth[i, sample(1:K, k[i])] <- TRUE
  ttt <- rep(0, G)
  ttt[1:(n.clust * clust.size)] <- rep(1:n.clust, clust.size)
  clust <- sample(ttt)
  result <- list()
  for(i in 1:K) {
    print(i)
    data0 <- matrix(rnorm(2 * G * N, 0, sigma), G, N * 2)
    for(j in 1:n.clust) {
      data0[clust == j, ] <- sample.correlated.data(clust.size, N * 2, rho.prior, df.prior, sigma)
    }
    mu.i <- mu[, i]
    mu.i[!truth[, i]] <- 0
    data0[, (1:N) + N] <- data0[, (1:N) + N] + mu.i
    result[[i]] <- data0
  }
  return(list(data = result, truth = truth))
}
#=====================================================================================
metalasso <- function(X.all, Y.all, obs, lam1, method, maxit, tol){
  ## starting and ending index of each dataset
  start.idx   <- cumsum(obs) + 1 - obs    # starting index of each dataset
  end.idx     <- cumsum(obs)              # ending index of each dataset
  M           <- length(start.idx)        # number of datasets
  p           <- ncol(X.all)              # number of covariates
  N           <- sum(obs)                 # total number of obserations
  gamma       <- matrix(NA, p, maxit + 1) # iterations of gamma
  gamma[, 1]  <- rep(1, p)
  theta       <- matrix(NA, M * p, maxit) # iterations of theta
  X.tha       <- matrix(NA, N, p)         # colMultiply(X.all, theta)
  beta.hat    <- vector("list", M)        # iterations of beta.hat
  coef        <- vector("list", M)        # final estimate of coefficients
  itr         <- 1
  m.diff      <- NULL                     # marginal error
  for (m in 1:M){
    beta.hat[[m]] <- matrix(NA, p, maxit)
  }
  if (method == "glmnet") {
    {
      while(!(itr > maxit)){
        ## Iterate as: theta --> gamma --> theta

        for (m in 1:M){
          ## In each dataset, fit Y.all ~ colMultiply(X.all, gamma*rho)
          theta.fit <- glmnet(t(t(X.all[start.idx[m]:end.idx[m], ]) *
                                  gamma[, itr]),
                              Y.all[start.idx[m]:end.idx[m]],
                              alpha = 1,
                              family = "binomial",
                              lambda = lam1,
                              standardize = FALSE)
          theta[((m - 1) * p + 1):(m * p), itr] <- as.vector(theta.fit$beta)
          ## adjust X.all by colMultiply(X.all, theta) for further usage
          X.tha[start.idx[m]:end.idx[m], ] <- t(t(X.all[start.idx[m]:end.idx[m], ]) *
                                                  as.vector(theta.fit$beta))
          beta.hat[[m]][, itr] <- theta[((m - 1) * p + 1):(m * p), itr] * gamma[, itr]
          ## calculate iteration difference
          if (itr == 1) {
            m.diff[m] <- max(abs(beta.hat[[m]][, itr]))
          }
          else {
            m.diff[m] <- max(abs(beta.hat[[m]][, itr] - beta.hat[[m]][, itr - 1]))
          }
        }
        if(max(m.diff) < tol) break         # break iterations if diff < tol
        itr <- itr + 1
        gamma.fit <- glmnet(X.tha, Y.all, alpha = 1, family = "binomial",
                            lambda = 1,
                            weights = rep(1 / obs, obs),
                            standardize = FALSE)
        gamma[, itr] <- as.vector(gamma.fit$beta)
      }
    }
  }
  if (method == "penalized") {
    while(!(itr > maxit)){
      ## Iterate as: theta --> gamma --> theta
      for (m in 1:M){
        ## In each dataset, fit Y.all ~ colMultiply(X.all, gamma*rho)
        theta.fit <- penalized(Y.all[start.idx[m]:end.idx[m]],
                               t(t(X.all[start.idx[m]:end.idx[m], ]) *
                                   gamma[, itr]),
                               positive = TRUE,
                               unpenalized = ~0, model = "logistic",
                               lambda1 = lam1,
                               standardize = FALSE, trace = FALSE)
        theta[((m - 1) * p + 1):(m * p), itr] <- theta.fit@penalized
        ## adjust X.all by colMultiply(X.all, theta) for further usage
        X.tha[start.idx[m]:end.idx[m], ] <- t(t(X.all[start.idx[m]:end.idx[m], ]) *
                                                theta.fit@penalized)
        beta.hat[[m]][, itr] <- theta[((m - 1) * p + 1):(m * p), itr] * gamma[, itr]
        ## calculate iteration difference
        if (itr == 1) {
          m.diff[m] <- max(abs(beta.hat[[m]][, itr]))
        }
        else {
          m.diff[m] <- max(abs(beta.hat[[m]][, itr] - beta.hat[[m]][, itr - 1]))
        }
      }
      if(max(m.diff) < tol) break         # break iterations if diff < tol
      itr <- itr + 1
      gamma.fit <- penalized(Y.all, X.tha,
                             unpenalized = ~0, model = "logistic",
                             lambda1 = 0.005,
                             standardize = FALSE, trace = FALSE)
      gamma[, itr] <- gamma.fit@penalized
    }
  }
  ## determine if convergence is achieved
  if (itr == 1) {
    iteration <- itr
    converge  <- FALSE
  }
  else {
    if (itr > maxit) {
      iteration <- itr - 1
      converge  <- FALSE
    }
    else {
      iteration <- itr
      converge  <- TRUE
    }
  }
  for (m in 1:M){
    coef[[m]] <- beta.hat[[m]][, iteration]
  }
  return(list(coe         = coef,
              gamma       = gamma[, iteration],
              iteration   = iteration,
              converge    = converge,
              lam         = lam1,
              diff        = max(m.diff)
  ))
}
#=====================================================================================
