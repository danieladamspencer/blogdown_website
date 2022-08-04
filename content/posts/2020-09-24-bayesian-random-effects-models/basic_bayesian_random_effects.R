# This is based off of the model W_ijk = X_ik + E_ijk
latent_random_effects_data <- function(n,J,K, mu, tau2, sigma2) {
  E <- array(rnorm(n*J*K, sd = sqrt(sigma2)), dim = c(n,J,K))
  X <- sapply(seq(K), function(k) {
    X_k <- rnorm(n, mean = mu[k], sd = sqrt(tau2[k]))
    return(X_k)
  })
  means <- X %o% array(1,dim = J)
  W <- aperm(means,perm = c(1,3,2)) + E
  return(W)
}

sample_X <- function(W,tau2,sigma2, prior_means) {
  N <- dim(W)[1]
  J <- dim(W)[2]
  K <- dim(W)[3]
  post_var <- 1 / (1/tau2 + rep(J/sigma2,K))
  post_mean <- matrix(prior_means / tau2,N,K,byrow = T) + apply(W,c(1,3),sum)/sigma2
  post_mean <- apply(post_mean,1,`*`,y = post_var)
  X_sample <- apply(post_mean,2,function(mu_xi) {
    rnorm(K,mean = mu_xi, sd = sqrt(post_var))
  })
  return(t(X_sample))
}

sample_tau2 <- function(X, prior_means, alpha_tau2, beta_tau2, 
                          variance_priors = "IG") {
  N <- dim(X)[1]
  K <- dim(X)[2]
  if(variance_priors == "Jeffreys") alpha_tau2 <- beta_tau2 <- 0
  post_shape <- alpha_tau2 + N/2
  post_rate <- beta_tau2 + apply((apply(X,2, `-`, y = prior_means))^2,2,sum)/2
  tau2_sample <- 1/rgamma(K,shape = post_shape, rate = post_rate)
  return(tau2_sample)
}

sample_sigma2 <- function(W, X, alpha_e, beta_e, variance_priors = "IG") {
  n <- dim(W)[1]
  J <- dim(W)[2]
  K <- dim(W)[3]
  means <- X %o% array(1,dim = J)
  sum_sq_err <- sum((W - aperm(means, perm = c(1,3,2)))^2)
  if(variance_priors == "Jeffreys") alpha_e <- beta_e <- 0
  sigma2_sample <- 1/rgamma(1,alpha_e + n*J*K/2, beta_e + sum_sq_err/2)
  return(sigma2_sample)
}

IG_param_values <- function(m,v) {
  return(matrix(c(m^2 / v + 2, m*(m^2 / v + 1)), nrow = length(m), ncol = 2, byrow = F))
}

latent_random_effects_mcmc <- function(W, n_iter = 1e4, n_burn = 5e2, 
                                       variance_priors = "Jeffreys",
                                       hyperparameters = NULL) {
  if(!variance_priors %in% c("Jeffreys","IG","eb_IG")) 
    stop("The variance_priors argument should be set to 'Jeffreys', 'IG', or 'eb_IG'.")
  # Set hyperparameters
  N <- dim(W)[1]
  J <- dim(W)[2]
  K <- dim(W)[3]
  prior_means = rep(0,K)
  alpha_tau2 = rep(2,K)
  beta_tau2  = rep(2,K)
  alpha_e = 2
  beta_e = 2
  tau2_prior_vars <- rep(5,K)
  sigma2_prior_var <- 5
  if (!is.null(hyperparameters)) {
    list2env(hyperparameters, envir = environment())
  }
  # Allocate parameter storage
  results <- list(
    X = array(NA,dim = c(N,K,n_iter)),
    tau2 = matrix(NA,K,n_iter),
    sigma2 = vector('numeric',n_iter)
  )
  # Set initials
  J <- dim(W)[2]
  X <- apply(W,c(1,3),mean)
  tau2 <- apply(X,2,var)
  sigma2 <- var(W - aperm(X%o%rep(1,J),perm = c(1,3,2)))
  if(variance_priors == "eb_IG") {
    tau2_hyperpars <- IG_param_values(tau2,tau2_prior_vars)
    alpha_tau2 <- tau2_hyperpars[,1]
    beta_tau2 <- tau2_hyperpars[,2]
    sigma2_hyperpars <- IG_param_values(sigma2, sigma2_prior_var)
    alpha_e <- sigma2_hyperpars[1]
    beta_e <- sigma2_hyperpars[2]
  }
  # Run MCMC
  cat("Starting MCMC\n")
  start_time <- proc.time()[3]
  for(s in seq(n_iter)) {
    X <- sample_X(W = W,tau2 = tau2,sigma2 = sigma2,
                  prior_means = prior_means)
    tau2 <- sample_tau2(alpha_tau2 = alpha_tau2,
                            beta_tau2 = beta_tau2,X = X,
                            prior_means = prior_means)
    sigma2 <- sample_sigma2(W = W, X = X, alpha_e = alpha_e, beta_e = beta_e)
    results$X[,,s] <- X
    results$tau2[,s] <- tau2
    results$sigma2[s] <- sigma2
  }
  total_time <- proc.time()[3] - start_time
  cat("MCMC finished,", n_iter,"samples in",total_time,"seconds\n")
  # Discard burn-in
  results$X <- results$X[,,-seq(n_burn)]
  results$tau2 <- results$tau2[,-seq(n_burn)]
  results$sigma2 <- results$sigma2[-seq(n_burn)]
  return(results)
}
