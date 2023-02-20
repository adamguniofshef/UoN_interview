# function for running forward model, i.e., estimating relationship between
# the output variables and the input variables

# in this case looking at errors-in-variables multi-output Gaussian process,
# estimating hyperparameters using MCMC

# firstly, need to define functions for covariance kernel

calcSE3_ARD_square <- function(X1,X2,l=rep(1,ncol(X1)),tau=1) {
  if(is.vector(X1)&is.vector(X2)){
    Sigma <- tau^2 * exp(-0.5*(sum((X1-X2)^2/l^2)))
  }else{
    Sigma <- matrix(rep(0, nrow(X1)*nrow(X2)), nrow=nrow(X1), ncol=nrow(X2))
    Sigma[nrow(X1),nrow(X2)] <- tau^2
    if(nrow(Sigma)>1){
      for (i in 1:(nrow(Sigma)-1)) {
        Sigma[i,i] <- tau^2
        for (j in (i+1):ncol(Sigma)) {
          Sigma[i,j] <- tau^2 * exp(-0.5*(sum((X1[i,]-X2[j,])^2/l^2)))
          Sigma[j,i] <- Sigma[i,j]
        }
      }
    }
  }
  return(Sigma)
}

calcSE3_ARD_nonsquare <- function(X1,X2,l=rep(1,ncol(X1)),tau=1) {
  if(is.vector(X1)){
    Sigma <- rep(NA,nrow(X2))
    for(i in 1:nrow(X2)){
      Sigma[i] <- tau^2 * exp(-0.5*(sum((X1-X2[i,])^2/l^2)))
    }
  }else{
    if(is.vector(X2)){
      Sigma <- rep(NA,nrow(X1))
      for(i in 1:nrow(X1)){
        Sigma[i] <- tau^2 * exp(-0.5*(sum((X1[i,]-X2)^2/l^2)))
      }
    }else{
      Sigma <- matrix(rep(0, nrow(X1)*nrow(X2)), nrow=nrow(X1), ncol=nrow(X2))
      for (i in 1:nrow(Sigma)) {
        for (j in 1:ncol(Sigma)) {
          Sigma[i,j] <- tau^2 * exp(-0.5*(sum((X1[i,]-X2[j,])^2/l^2)))
        }
      }
      
    }
  }
  return(Sigma)
}




EIVMOGP_forward <- function(Y1, Y2, X_joint = NA, X_marginal = NA, a_X = 1e-03, b_X = 1e-03, a_Xt = 1e-03, b_Xt = 1e-03, a= 1e-03, b = 1e-03,
                          a_Y = 1e-03, b_Y = 1e-03, n_m = 7, n_iter = 200000, n_thin = 10,
                          simulation = F, Yt1, Yt2, Xt, informed = F, seed = 999,
                          scaling = T, n_chains = 2,
                          l_shape = rep(2, dim.X), l_mean = rep(0.25, dim.X), sigma_shape = 1,
                          sigma_mean = 0.1, tau_shape = 3, tau_mean = 0.4, n_adapt = 1000, n_burn = 25000,
                          mixed_LOO = F, marginal_resp_obs = T, dim.X = 2, dim.X.joint = 0, dim.X.marginal = 0,
                          initialisation = T, mu_alpha = rep(0.5, 2), tau_alpha = rep(10.8241, 2),
                          a_epsilon = 1e-03, b_epsilon = 1e-03, S_delta, S_Xt, mu_X = rep(0.5, dim.X),
                          Y_range){
  # Y1: vector of observed data for first output
  # Y2: vector of observed data for second output
  # X_joint: matrix of observed data for inputs whose replicates are measured under same conditions
  # X_marginal: matrix of observed data for inputs whose replicates are measured under different conditions
  # a_X: shape parameter for gamma distribution assumed for within-groups precision of inputs
  # b_X: rate parameter for gamma distribution assumed for within-groups precision of inputs
  # a_Xt: shape parameter for gamma distribution assumed for between-groups precision of inputs
  # b_Xt: rate parameter for gamma distribution assumed for between-groups precision of inputs
  # a_Y: shape parameter for gamma distribution assumed for within-groups precision of outputs 
  # b_Y: rate parameter for gamma distribution assumed for within-groups precision of outputs
  # n_m: number of groups in the data
  # n_iter: number of samples drawn from posterior
  # n_thin: number corresponding to storing every n_thin^{th} sample from n_iter
  # simulation: whether model is fitted using simulated data
  # Yt1: true values of first output (if known, i.e., if simulation = T)
  # Yt2: true values of second output (if known, i.e., if simulation = T)
  # Xt: true values of inputs (if known, i.e., if simulation = T)
  # informed: whether informed priors are used for estimating parameters
  # seed: seed for RNG (reproducibility)
  # scaling: whether to scale the data
  # n_chains: number of parallel chains for which the model is fitted
  # l_shape: shape parameter for gamma distribution assumed for length-scale parameter in kernel
  # l_mean: mean of gamma distribution assumed for length-scale parameter in kernel
  # sigma_shape: shape parameter for gamma distribution assumed for noise sd
  # sigma_mean: mean of gamma distribution assumed for noise sd
  # tau_shape: shape parameter for gamma distribution assumed for signal sd
  # tau_mean: mean of gamma distribution assumed for signal sd
  # n_adapt: number of iterations for adaptation of MCMC algo
  # n_burn: number of iterations for burn-in phase
  # mixed_LOO: whether the mixed LOO-CV should be used to assess model fit
  # marginal_resp_obs: whether replicates for two outputs are measured under same conditions
  # dim.X: number of inputs
  # dim.X.joint: number of inputs whose replicates are measured under same conditions 
  # dim.X.marginal: number of inputs whose replicates are measured under different conditions
  # initialisation: whether or not to provide starting values for each parameter in each chain of MCMC
  # mu_alpha: mean of GP prior mean
  # tau_alpha: precision of GP prior mean
  # a_epsilon: shape of gamma dist assumed for model error PRECISION (depending on preference of parameterisation)
  # b_epsilon: rate of gamma dist assumed for model error PRECISION (depending on preference of parameterisation)
  # S_delta: scale matrix for informed prior for within-groups precision of inputs
  # S_Xt: scale matrix for informed prior for between-groups precision of inputs
  # mu_X: mean of true values for input
  # Y_range: range of possible values for output variables (for informed priors)
  
  require(rjags)
  require(car)
  require(lattice)    # load required packages
  require(loo)
  require(mvtnorm)
  n_r <- rep(NA,2)
  n_r[1] <- length(na.omit(Y1)) / n_m         # number of replicate measurements for first output
  n_r[2] <- length(na.omit(Y2)) / n_m         # same for second output
  if(marginal_resp_obs){
    dim.Y.marginal <- 2
    dim.Y.joint <- 0
  }else{
    dim.Y.marginal <- 0
    dim.Y.joint <- 2
  }
  cov_samples_dist <- "both"
  if(dim.X == 0){                    # determining whether replicate measurements for inputs should be
    cov_samples_dist <- NA           # jointly distributed, marginally distributed, or both
    n_c <- NA
  }else{
    if(dim.X.joint == 0){
      n_c <- rep(NA,dim.X.marginal)
      for(k in 1:dim.X.marginal){
        n_c[k] <- nrow(na.omit(X_marginal[,k])) / n_m
      }
      cov_samples_dist <- "marginal"
    }
    if(dim.X.marginal == 0){
      n_c <- nrow(na.omit(X_joint[,1])) / n_m
      cov_samples_dist <- "joint"
    }
    if(dim.X.joint > 0 & dim.X.marginal > 0){
      n_c <- rep(NA, 1+dim.X.marginal)
      for(k in 1:dim.X.marginal){
        n_c[k] <- nrow(na.omit(X_marginal[,k])) / n_m 
      }
      n_c[1+dim.X.marginal] <- nrow(na.omit(X_joint[,1])) / n_m 
    }
  }
  dim.Y <- 2
  nu_X <- dim.X.joint; nu_Xt <- dim.X # degrees of freedom for precision matrices
  nu_Y <- dim.Y.joint
  obs_scaling <- rep(NA,dim.X+dim.Y)      # scaling the data variables, onto the range [0,1], using maximum
  maxs <- rep(NA, dim.X+dim.Y)
  maxs[1] <- max(na.omit(Y1))
  maxs[2] <- max(na.omit(Y2))
  obs_scaling[1] <- if(floor(maxs[1])==0){0}else{nchar(floor(maxs[1]))}
  obs_scaling[2] <- if(floor(maxs[2])==0){0}else{nchar(floor(maxs[2]))}
  if(marginal_resp_obs){
    Y__joint <- NA
    Y__marginal <- array(NA, dim = c(n_m, max(n_r), dim.Y.marginal))
    Y__marginal[,1:n_r[1],1] <- matrix(na.omit(Y1), n_m, n_r[1], byrow = T) / (10^obs_scaling[1])
    Y__marginal[,1:n_r[2],2] <- matrix(na.omit(Y2), n_m, n_r[2], byrow = T) / (10^obs_scaling[2])
  }else{
    Y__marginal <- NA
    Y__joint <- array(NA, dim = c(n_m, max(n_r), dim.Y.joint))
    Y__joint[,,1] <- matrix(na.omit(Y1), n_m, max(n_r), byrow = T) / (10^obs_scaling[1])
    Y__joint[,,2] <- matrix(na.omit(Y2), n_m, max(n_r), byrow = T) / (10^obs_scaling[2])
  }
  if(dim.X.joint == 0){
    X__joint <- NA
    mu_X_joint <- NA
  }else{
    X__joint <- array(NA, dim = c(n_m, n_c[1+dim.X.marginal], dim.X.joint))
    mu_X_joint <- rep(NA, dim.X.joint)
    for(l in 1:dim.X.joint){
      
      maxs[l+2+dim.X.marginal] <- max(na.omit(as.data.frame(X_joint)[,l]))
      obs_scaling[l+2+dim.X.marginal] <- if(floor(maxs[l+2+dim.X.marginal])==0){0}else{nchar(floor(maxs[l+2+dim.X.marginal]))}
      X__joint[,,l] <- matrix(na.omit(as.data.frame(X_joint)[,l]), n_m, n_c[1+dim.X.marginal], byrow = T)/(10^obs_scaling[l+2+dim.X.marginal])
      mu_X_joint[l] <- mean(na.omit(as.data.frame(X_joint)[,l]))/(10^obs_scaling[l+2+dim.X.marginal])
      
    }
  }
  if(dim.X.marginal == 0){
    X__marginal <- NA
    mu_X_marginal <- NA
  }else{
    if(dim.X.joint == 0){
      X__marginal <- array(NA, dim = c(n_m, max(n_c), dim.X.marginal))
      mu_X_marginal <- rep(NA, dim.X.marginal)
      
    }else{
      X__marginal <- array(NA, dim = c(n_m, max(n_c[-(2+dim.X.marginal)]), dim.X.marginal))
      mu_X_marginal <- rep(NA, dim.X.marginal)
    }
    
    for(l in 1:dim.X.marginal){
      
      maxs[l+2] <- max(na.omit(as.data.frame(X_marginal)[,l]))
      obs_scaling[l+2] <- if(floor(maxs[l+2])==0){0}else{nchar(floor(maxs[l+2]))}
      X__marginal[,1:n_c[l],l] <- matrix(na.omit(as.data.frame(X_marginal)[,l]), n_m, n_c[l], byrow = T)/(10^obs_scaling[l+2])
      mu_X_marginal[l] <- mean(na.omit(as.data.frame(X_marginal)[,l]))/(10^obs_scaling[l+2])
      
    }
    
    
  }
  if(simulation){   # scaling the simulated data
    for(d_X in 1:dim.X){
      Xt[,d_X] <- Xt[,d_X] / 10^obs_scaling[d_X+2]
    }
    Yt1 <- Yt1 / 10^obs_scaling[1]
    Yt2 <- Yt2 / 10^obs_scaling[2]
    Xt[, 1] <- Xt[, 1] * 10 # NOTE: manual override of data scaling (to reduce model error)
  }
  X__joint[, , 1] <- X__joint[, , 1] * 10 # NOTE: manual override of data scaling (to reduce model error)
  if(informed == F){   # uninformed priors for model parameters
    
    a_Y <- 1e-03
    b_Y <- 1e-03
    a_Xt <- 1e-03
    b_Xt <- 1e-03
    a_X <- 1e-03
    b_X <- 1e-03
    scaleX_0 <- (1/factorX_0)*diag(dim.X.joint)
    scaleXt_0 <- (1/factorXt_0)*diag(dim.X)
    mu_X <- rep(0.5, dim.X)
    mu_alpha <- rep(0.5, 2)
    tau_alpha <- rep(10.8241, 2) # assuming that mean of each response is within [0,1]
    # with 90% probability.
    
    
  }else{
    scaleX_0 <- solve(diag(S_delta))
    scaleXt_0 <- solve(diag(S_Xt))
  }
  
  I <- diag(n_m)
  
  #scale_T <- diag(1, 2)
  #nu_T <- 2
  data_list <- list('a_X' = a_X, 'b_X' = b_X, 'a_Y' = a_Y, 'b_Y' = b_Y, 'a_Xt' = a_Xt, 'b_Xt' = b_Xt,
                    'n_m' = n_m, 'Y__marginal' = Y__marginal,
                    'X__marginal' = X__marginal, 'Y__joint' = Y__joint,
                    'mu_alpha' = mu_alpha, 'tau_alpha' = tau_alpha,
                    'a_epsilon' = a_epsilon, 'b_epsilon' = b_epsilon,
                    'X__joint' = X__joint, 'nu_Xt' = nu_Xt, 'nu_X' = nu_X, 'dim.X.joint' = dim.X.joint,
                    'dim.X.marginal' = dim.X.marginal, 'scaleX_0' = scaleX_0, 'dim.Y.joint' = dim.Y.joint, 'dim.Y.marginal' = dim.Y.marginal,
                    'n_r' = n_r, 'n_c' = n_c, 'I' = I, 'mu_X' = mu_X[1:dim.X], 'a' = a, 'b' = b, 'l_shape' = l_shape,
                    'l_mean' = l_mean, 'scaleXt_0' = scaleXt_0, 'dim.X' = dim.X, 'tau_shape' = tau_shape, 'tau_mean' = tau_mean, 'sigma_shape' = sigma_shape,
                    'sigma_mean' = sigma_mean) # assigning values for nodes in JAGS file
  
  inits_list <- vector("list",n_chains)
  
  if(initialisation){ # choosing overdispersed starting points for each MCMC chain
    
    # which parameters need initialising?
    # the GP means of alpha_1, alpha_2. 
    # the signal standard deviations sigma_k = (sigma_k1,sigma_k2)
    # the noise precisions tau_{epsilon,1}, tau_{epsilon,2}
    # the distance scaling parameters l_{1,1}, l_{1,2}, l_{2,1}, l_{2,2}
    # the precisions tau_Y_marginal[1], tau_Y_marginal[2], tau_Xt, tau_X_joint
    # the correlation rho_T, the scalars lambda_1, lambda_2
    # the true values Yt, Xt
    
    lhs_k <- 2 + 2 + 2 + 4 + 2 + 1 + 1 + 2*n_m + dim.X*n_m + 3
    require(lhs)
    set.seed(999)
    lhs_init <- randomLHS(n = n_chains, k = lhs_k)
    # matrices to initialise are tau_X_joint, tau_Xt
    
    if(dim.X.joint > 0){
      
      tau_X_joint_samp <- rWishart(100000, df = dim.X.joint,
                                   Sigma = S_delta)
      tau_X_joint_det <- rep(NA, 100000)
      
    }
    
    
    tau_Xt_samp <- rWishart(100000, df = dim.X,
                            Sigma = S_Xt)
    tau_Xt_det <- rep(NA, 100000)
    
    for(n in 1:100000){
      
      tau_X_joint_det[n] <- det(tau_X_joint_samp[, , n])
      tau_Xt_det[n] <- det(tau_Xt_samp[, , n])
      
    }
    
    tau_Y_marginal_init <- matrix(NA, dim.Y.marginal, n_chains)
    tau_X_joint_init <- array(NA, dim = c(dim.X.joint, dim.X.joint, n_chains))
    tau_Xt_init <- array(NA, dim = c(dim.X, dim.X, n_chains))
    
    Xt_init <- array(NA, dim = c(n_m, dim.X, n_chains))
    Yt_init <- array(NA, dim = c(n_m, dim.Y, n_chains))
    
    l_1_init <- matrix(NA, dim.X, n_chains)
    l_2_init <- matrix(NA, dim.X, n_chains)
    
    sigma_k_init <- matrix(NA, 2, n_chains)
    tau_epsilon_init <- matrix(NA, 2, n_chains)
    alpha_init <- matrix(NA, 2, n_chains)
    
    rho_T_init <- rep(NA, n_chains)
    lambda_init <- matrix(NA, 2, n_chains)
    
    
    for(n in 1:n_chains){
      for(d_Y in 1:dim.Y.marginal){
        
        tau_Y_marginal_init[d_Y, n] <- qgamma(p = lhs_init[n, d_Y],
                                              shape = a_Y[d_Y],
                                              rate = b_Y[d_Y])
        
      }
      for(d_Y in 1:2){
        
        tau_epsilon_init[d_Y, n] <- qgamma(p = lhs_init[n, d_Y + dim.Y.marginal],
                                           shape = a_epsilon[d_Y],
                                           rate = b_epsilon[d_Y])
        
      }
      
      q1 <- quantile(tau_X_joint_det, probs = lhs_init[n, 5])
      tau_X_joint_init[, , n] <- tau_X_joint_samp[, , which.min(abs(q1 - tau_X_joint_det))]
      q2 <- quantile(tau_Xt_det, probs = lhs_init[n, 6])
      tau_Xt_init[, , n] <- tau_Xt_samp[, , which.min(abs(q2 - tau_Xt_det))]
      
      for(d_X in 1:dim.X){
        
        l_1_init[d_X, n] <- qgamma(p = lhs_init[n, d_X + 6],
                                   shape = l_shape[d_X], rate = l_shape[d_X] / l_mean[d_X])
        
      }
      
      for(d_X in 1:dim.X){
        
        l_2_init[d_X, n] <- qgamma(p = lhs_init[n, d_X + 6 + dim.X],
                                   shape = l_shape[d_X], rate = l_shape[d_X] / l_mean[d_X])
        
      }
      
      for(d_Y in 1:dim.Y){
        
        sigma_k_init[d_Y, n] <- qgamma(p = lhs_init[n, d_Y + 6 + 2*dim.X],
                                       shape = tau_shape, rate = tau_shape / tau_mean)
        
      }
      
      for(d_Y in 1:dim.Y){
        
        alpha_init[d_Y, n] <- qnorm(p = lhs_init[n, d_Y + 6 + 2*dim.X + dim.Y],
                                    mean = mu_alpha[d_Y],
                                    sd = sqrt(1 / tau_alpha[d_Y]))
        
      }
      
      for(i in 1:n_m){
        
        for(d_X in 1:dim.X){
          
          Xt_init[i, d_X, n] <- qnorm(lhs_init[n, d_X + (i-1)*dim.X + 6 + 2*dim.X + 2*dim.Y],
                                      mean = mu_X[d_X],
                                      sd = sqrt(1 / S_Xt[d_X, d_X]))
          
        }
        
      }
      
      for(i in 1:n_m){
        
        for(d_Y in 1:dim.Y){
          
          Yt_init[i, d_Y, n] <- qnorm(lhs_init[n, d_Y + (i-1)*dim.Y + dim.X*n_m + 6 + 2*dim.X + 2*dim.Y],
                                      mean = mu_alpha[d_Y],
                                      sd = (Y_range[d_Y, 2] - Y_range[d_Y, 1]) / 3.92)
          
        }
        
      }
      
      rho_T_init[n] <- qunif(p = lhs_init[n, 1 + 2*n_m + dim.X*n_m + 6 + 2*dim.X + 2*dim.Y],
                             min = -1, max = 1)
      
      for(d_Y in 1:dim.Y){
        
        lambda_init[d_Y, n] <- qunif(p = lhs_init[n, d_Y + 1 + 2*n_m + dim.X*n_m + 6 + 2*dim.X + 2*dim.Y],
                                     min = 0, max = 5)
        
      }
      
    }
    
    
    for(c in 1:n_chains){ # storing the chains for reproducibility
      inits_list[[c]] <- list('tau_Y_marginal' = tau_Y_marginal_init[, c],
                              'tau_X_joint' = tau_X_joint_init[, , c],
                              'tau_Xt' = tau_Xt_init[, , c],
                              'tau_epsilon' = tau_epsilon_init[, c],
                              'Yt' = c(Yt_init[, , c]),
                              'Xt' = Xt_init[, , c],
                              'tau' = sigma_k_init[, c],
                              'l_1' = l_1_init[, c],
                              'l_2' = l_2_init[, c],
                              'alpha' = alpha_init[, c],
                              'rho_T' = rho_T_init[c],
                              'lambda' = lambda_init[, c],
                              .RNG.name = "base::Wichmann-Hill", .RNG.seed = 998+c)
    }
    
  }else{
    
    
    for(c in 1:n_chains){ # storing the chains for reproducibility
      inits_list[[c]] <- list(.RNG.name="base::Wichmann-Hill",.RNG.seed = 999+(c-1))
    }
    
  }
  
  

  GPmodel <- jags.model('MOGPEIVSEARD.bug', data = data_list, n.chains = n_chains, n.adapt = n_adapt,
                          inits = inits_list) # initialise the model, adaptation phase 
    init.state <- GPmodel$state(internal = FALSE)
    #}
    # return(list(init.state = init.state, X__joint = X__joint, Xt = Xt))
    #}
    update(GPmodel, n_burn) # burn-in
    posteriors <- coda.samples(GPmodel, c('tau_Y_marginal', 'tau_Y_joint', 'tau_Xt', 'Yt', 'Xt',
                                          'l_1', 'l_2', 'sigma', 'alpha', 'tau', 'tau_X_marginal',
                                          'tau_X_joint', 'rho_T', 'lambda'
                                          #,'k_l','k_tau'
    ), n.iter = n_iter,
    thin = n_thin) # sample from posterior distribution
  
  
  summary_post <- summary(posteriors)
  if(marginal_resp_obs){
    
    if(cov_samples_dist=="both"){
      tau_Xt_off_diag <- rep(NA, sum(1:(dim.X-1)))
      tau_Xt_off_diag[1:(dim.X - 1)] <- paste0("tau_Xt[", 1,",", 2:dim.X, "]")
      if(dim.X > 2){
        for(i in 2:(dim.X-1)){
          
          tau_Xt_off_diag[(1 + sum((dim.X - 1):((dim.X - 1) - (i - 2)))):((dim.X - 1)*i - sum(0:(i - 1)))] <- paste0("tau_Xt[", i,",", (i + 1):dim.X, "]")
          
        }
      }
      tau_X_joint_off_diag <- rep(NA, sum(1:(dim.X.joint - 1)))
      tau_X_joint_off_diag[1:(dim.X.joint - 1)] <- paste0("tau_X_joint[", 1,",", 2:dim.X.joint, "]")
      if(dim.X.joint > 2){
        for(i in 2:(dim.X.joint-1)){
          
          tau_X_joint_off_diag[(1 + sum((dim.X.joint - 1):((dim.X.joint - 1) - (i - 2)))):((dim.X.joint - 1)*i - sum(0:(i - 1)))] <- paste0("tau_X_joint[", i,",", (i + 1):dim.X.joint, "]")
          
        }
      }
      
      theta <- c(paste0("Xt[",rep(1:n_m,dim.X),",",rep(1:dim.X,rep(n_m,dim.X)),"]"),
                 paste0("Yt[",1:(2*n_m),"]"),
                 paste0("tau_X_joint[",1:dim.X.joint,",",1:dim.X.joint,"]"),
                 tau_X_joint_off_diag, "rho_T", "lambda[1]", "lambda[2]",
                 paste0("tau_X_marginal[",1:dim.X.marginal,"]"),paste0("tau_Y_marginal[",1:2,"]"),
                 paste0("tau_Xt[",1:dim.X,",",1:dim.X,"]"),tau_Xt_off_diag,
                 paste0("l_1[",1:dim.X,"]"), paste0("l_2[",1:dim.X,"]"),paste0("sigma[",1:2,"]"),paste0("alpha[",1:2,"]"),paste0("tau[",1:2,"]")
                 #,paste0("k_l[",rep(1:dim.X,2),",",rep(1:2,rep(dim.X,2)),"]"),paste0('k_tau[',1:2,"]")
      ) # parameter vector
    }
    if(cov_samples_dist=="marginal"){
      theta <- c(paste0("Xt[",rep(1:n_m,dim.X),",",rep(1:dim.X,rep(n_m,dim.X)),"]"),
                 paste0("Yt[",1:(2*n_m),"]"), "lambda[1]", "lambda[2]",
                 paste0("tau_X_marginal[",1:dim.X,"]"),paste0("tau_Y_marginal[",1:2,"]"),
                 paste0("tau_Xt[",1:dim.X,",",1:dim.X,"]"),"tau_Xt[1,2]","rho_T",
                 paste0("l_1[",1:dim.X,"]"), paste0("l_2[",1:dim.X,"]"),paste0("sigma[",1:2,"]"),paste0("alpha[",1:2,"]"),paste0("tau[",1:2,"]")
                 ,paste0("k_l[",rep(1:dim.X,2),",",rep(1:2,rep(dim.X,2)),"]"),paste0('k_tau[',1:2,"]")
      )
    }
    if(cov_samples_dist=="joint"){
      theta <- c(paste0("Xt[",rep(1:n_m,dim.X),",",rep(1:dim.X,rep(n_m,dim.X)),"]"),
                 paste0("Yt[",1:(2*n_m),"]"), "lambda[1]", "lambda[2]"
                 #"lambda"
                 ,"rho_T"
                 #,paste0("rho_T[", 1:2, "]")
                 ,paste0("tau_X_joint[",1:dim.X,",",1:dim.X,"]"),"tau_X_joint[1,2]",paste0("tau_Y_marginal[",1:2,"]")
                 ,paste0("tau_Xt[",1:dim.X,",",1:dim.X,"]"),"tau_Xt[1,2]"
                 ,paste0("sigma[",1:2,"]"),paste0("alpha[",1:2,"]"),
                 paste0("tau[",1:2,"]")
                 #paste0("tau[",rep(1:2, 2), ",", rep(1:2, rep(2, 2)),"]")
                 ,paste0("l_1[",1:dim.X,"]"), paste0("l_2[",1:dim.X,"]")
      )
    }
    
  }else{
    
    if(cov_samples_dist=="both"){
      theta <- c(paste0("Xt[",rep(1:n_m,dim.X),",",rep(1:dim.X,rep(n_m,dim.X)),"]"),
                 paste0("Yt[",1:(2*n_m),"]"),
                 paste0("tau_X_joint[",1:dim.X.joint,",",1:dim.X.joint,"]"),
                 "tau_X_joint[1,2]","rho_T", "lambda",
                 paste0("tau_X_marginal[",1:dim.X.marginal,"]"),paste0("tau_Y_joint[",1:2,",",1:2,"]"),
                 "tau_Y_joint[1,2]",
                 paste0("tau_Xt[",1:dim.X,",",1:dim.X,"]"),"tau_Xt[1,2]",
                 paste0("l_1[",1:dim.X,"]"), paste0("l_2[",1:dim.X,"]"),paste0("sigma[",1:2,"]"),paste0("alpha[",1:2,"]"),paste0("tau[",1:2,"]")
                 ,paste0("k_l[",rep(1:dim.X,2),",",rep(1:2,rep(dim.X,2)),"]"),paste0('k_tau[',1:2,"]")
      ) # parameter vector
    }
    if(cov_samples_dist=="marginal"){
      theta <- c(paste0("Xt[",rep(1:n_m,dim.X),",",rep(1:dim.X,rep(n_m,dim.X)),"]"),
                 paste0("Yt[",1:(2*n_m),"]"),
                 paste0("tau_X_marginal[",1:dim.X,"]"),paste0("tau_Y_joint[",1:2,",",1:2,"]"),
                 "tau_Y_joint[1,2]","rho_T", "lambda",
                 paste0("tau_Xt[",1:dim.X,",",1:dim.X,"]"),"tau_Xt[1,2]",
                 paste0("l_1[",1:dim.X,"]"), paste0("l_2[",1:dim.X,"]"),paste0("sigma[",1:2,"]"),paste0("alpha[",1:2,"]"),paste0("tau[",1:2,"]")
                 ,paste0("k_l[",rep(1:dim.X,2),",",rep(1:2,rep(dim.X,2)),"]"),paste0('k_tau[',1:2,"]")
      )
    }
    if(cov_samples_dist=="joint"){
      theta <- c(paste0("Xt[",rep(1:n_m,dim.X),",",rep(1:dim.X,rep(n_m,dim.X)),"]"),
                 paste0("Yt[",1:(2*n_m),"]")
                 ,paste0("tau_X_joint[",1:dim.X,",",1:dim.X,"]"),"tau_X_joint[1,2]",paste0("tau_Y_joint[",1:2,",",1:2,"]"),"rho_T",
                 "tau_Y_joint[1,2]"
                 ,paste0("tau_Xt[",1:dim.X,",",1:dim.X,"]"),"tau_Xt[1,2]",paste0("l[",1:dim.X,"]")
                 ,paste0("sigma[",1:2,"]"),paste0("alpha[",1:2,"]"),paste0("tau[",1:2,"]")
                 ,paste0("k_l[",rep(1:dim.X,2),",",rep(1:2,rep(dim.X,2)),"]"),paste0('k_tau[',1:2,"]")
      )
    }
    
  }
  
  
  posteriors_list <- vector("list",n_chains)
  for(c in 1:n_chains){
    posteriors_list[[c]] <- posteriors[[c]][,theta]
  }
  posteriors_list <- mcmc.list(posteriors_list)
  psrf <- gelman.diag(posteriors_list) # checking convergence to posterior
  
  # assessing model fit using mixed LOO-CV
  
  if(LOO_CV_IC){
    if(mixed_LOO){
      # looking to implement a mixed LOO-CV-IC
      # calculating the computed lppd_{loo-cv}
      # sum_{i=1}^{n_g} log( (1/S) * sum_{s=1}^{S} p(Yt_i|theta_s,D_{-i}) )
      # i.e., only one posterior distribution for the hyperparameters
      # and true values (based on all the data), but the exact LOO-CV-IC
      # is carried out, so finding the likelihood of Yt_i given the
      # rest of the data D_{-i}=(Yt_{-i},Xt_{-i})
      
      lik <- array(NA,dim=c(nrow=n_iter/n_thin,n_m,n_chains)) # storing the likelihoods
      out.loo <- rep(NA,n_chains) # storing the computed lppd_{loo-cv} for each chain
      n_fail <- rep(NA,n_chains) # tracking how many samples fail
      failed_samp <- vector("list",n_chains) # tracking which samples fail
      for(c in 1:n_chains){
        
        set.seed(seed = seed)
        alpha_post <- posteriors[[c]][, paste0("alpha[", 1:2, "]")]
        sigma_post <- posteriors[[c]][, paste0("sigma[", 1:2, "]")]
        tau_post <- posteriors[[c]][, paste0("tau[", 1:2, "]")]
        Yt_post <- posteriors[[c]][, paste0("Yt[", 1:(2*n_m), "]")]
        Xt_post <- posteriors[[c]][, paste0("Xt[", rep(1:n_m, dim.X),
                                            ",", rep(1:dim.X, rep(n_m, dim.X)),"]")]
        l_1_post <- posteriors[[c]][, paste0("l_1[", 1:dim.X, "]")]
        l_2_post <- posteriors[[c]][, paste0("l_2[", 1:dim.X, "]")]
      
        
        Xt_train <- array(NA,dim = c(n_m, dim.X * (n_m - 1), (n_iter / n_thin)))
        Xt_test <- array(NA,dim = c(n_m, dim.X, (n_iter / n_thin )))
        
        m_Xt_train <- array(NA, dim=c(2 * n_m, (n_iter / n_thin)))
        m_Xt_test <- array(NA, dim=c(2 * n_m, (n_iter / n_thin)))
        
        T_1 <- array(NA, dim = c(2, 2, (n_iter / n_thin)))
        T_2 <- array(NA, dim = c(2, 2, (n_iter / n_thin)))
        
        S_Xt_train <- array(NA, dim = c(2 * (n_m - 1), 2 * (n_m - 1), (n_iter / n_thin), n_m))
        
        K_Xt_test_Xt_train <- array(NA,dim=c(2, 2 * (n_m - 1), (n_iter / n_thin), n_m))
        
        S_Xt_test <- array(NA,dim=c(2, 2, (n_iter / n_thin), n_m))
        
        m_star_Xt <- array(NA,dim=c(2, (n_iter / n_thin), n_m))
        S_star_Xt <- array(NA,dim=c(2, 2, (n_iter / n_thin), n_m))
        
        for(s in 1:(n_iter / n_thin)){
          
          m_Xt_train[1:n_m, s] <- rep(alpha_post[s, 1], n_m)
          m_Xt_train[(n_m + 1):(2 * n_m), s] <- rep(alpha_post[s, 2], n_m)
          m_Xt_test[1:n_m, s] <- rep(alpha_post[s, 1], n_m)
          m_Xt_test[(n_m + 1):(2 * n_m), s] <- rep(alpha_post[s, 2], n_m)
          

            
          T_1[, , s] <- matrix(c(tau_post[s, 1]^2,
                                 tau_post[s, 1] * tau_post[s, 2] * rho_T_post[s],
                                 tau_post[s, 1] * tau_post[s, 2] * rho_T_post[s],
                                  tau_post[s, 2]^2), 2, 2)
          
          T_2[, , s] <- matrix(c(lambda_post[s, 1]^2 * tau_post[s, 1]^2,
                                 lambda_post[s, 1] * lambda_post[s, 2] * tau_post[s, 1] * tau_post[s ,2] * rho_T_post[s],
                                 lambda_post[s, 1] * lambda_post[s, 2] * tau_post[s, 1] * tau_post[s, 2] * rho_T_post[s],
                                 lambda_post[s, 2]^2 * tau_post[s, 2]^2), 2, 2)

            
            
          } 
          
          for(i in 1:n_m){ # want to test the GP on each of the true values of input
            
            Xt_train[i, , s] <- Xt_post[s, -(n_m * (1:dim.X) - (n_m - i))]
            Xt_test[i, , s] <- Xt_post[s, (n_m * (1:dim.X) - (n_m - i))]
            
            S_Xt_train[, , s, i] <- diag(c(rep(sigma_post[s, 1]^2, n_m - 1),
                                           rep(sigma_post[s, 2]^2, n_m - 1)),
                                         nrow = 2 * (n_m - 1)) + kronecker(T_1[, , s], calcSE3_ARD_square(X1 = matrix(Xt_train[i, , s],
                                                                                                                      n_m - 1,
                                                                                                                      dim.X),
                                                                                                          X2 = matrix(Xt_train[i, , s],
                                                                                                                      n_m - 1,
                                                                                                                      dim.X),l = l_1_post[s,])) + kronecker(T_2[, , s], calcSE3_ARD_square(X1 = matrix(Xt_train[i, , s],
                                                                                                                                                                                                       n_m - 1,
                                                                                                                                                                                                       dim.X),
                                                                                                                                                                                           X2 = matrix(Xt_train[i, , s],
                                                                                                                                                                                                       n_m - 1,
                                                                                                                                                                                                       dim.X),l = l_2_post[s,]))
            
            
            K_Xt_test_Xt_train[, , s, i] <- kronecker(T_1[, , s], calcSE3_ARD_nonsquare(X1 = t(Xt_test[i, , s]),
                                                                                        X2 = matrix(Xt_train[i, , s],
                                                                                                    n_m - 1,
                                                                                                    dim.X),
                                                                                        l = l_1_post[s, ])) + kronecker(T_2[, , s], calcSE3_ARD_nonsquare(X1 = t(Xt_test[i, , s]),
                                                                                                                                                          X2 = matrix(Xt_train[i, , s],
                                                                                                                                                                      n_m - 1,
                                                                                                                                                                      dim.X),
                                                                                                                                                          l = l_2_post[s, ]))
            
            S_Xt_test[,,s,i] <- diag(c(sigma_post[s, 1]^2,
                                       sigma_post[s, 2]^2)) + kronecker(T_1[, , s], calcSE3_ARD_square(X1 = t(Xt_test[i, , s]),
                                                                                                       X2 = t(Xt_test[i, , s]),
                                                                                                       l = l_1_post[s, ])) + kronecker(T_2[, , s], calcSE3_ARD_square(X1 = t(Xt_test[i, , s]),
                                                                                                                                                                      X2 = t(Xt_test[i, , s]),
                                                                                                                                                                      l = l_2_post[s, ]))
            
            m_star_Xt[, s, i] <- m_Xt_test[c(i, i + n_m), s] + K_Xt_test_Xt_train[, , s, i] %*% solve(S_Xt_train[, , s, i]) %*% t(t(Yt_post[s,-c(i , i + n_m)]) - m_Xt_train[-c(i, i + n_m),s])
            S_star_Xt[, , s, i] <- S_Xt_test[, , s, i] - K_Xt_test_Xt_train[, , s, i] %*% solve(S_Xt_train[, , s, i]) %*% t(K_Xt_test_Xt_train[, , s, i])
            S_star_Xt[2, 1, s, i] <- S_star_Xt[1, 2, s, i]
          } 
        }
        negative_variance <- rep(NA,(n_iter/n_thin))
        for(s in 1:(n_iter / n_thin)){ #######
          negative_variance[s] <- any(S_star_Xt[1, 1, s, ] < 0) | any(S_star_Xt[2, 2, s, ] < 0)
        }
        failed_samp[[c]] <- which(negative_variance == T)
        n_fail[c] <- length(failed_samp[[c]])
        
        if(n_fail[c] > 0){
          # removing posterior samples for which S_star_Xt<0
          S_star_Xt <- S_star_Xt[, , -failed_samp[[c]], ]
          m_star_Xt <- m_star_Xt[, -failed_samp[[c]], ]
          Yt_post <- Yt_post[-failed_samp[[c]], ]
        }
        if((n_iter / n_thin) - n_fail[c] < 30){
          # if a large majority of the samples fail (leaving less than 30 'good' samples),
          # just want to return an error of sorts
          for(i in 1:n_m){
            for(s in 1:(n_iter / n_thin)){
              lik[s, i, c] <- NA
              out.loo[[c]] <- "Error: There are less than 30 samples for which S_star_Xt variances are positive in this chain, so the estimation of the LOO-CV-IC is aborted."
            }
          }
        }else{
          for(i in 1:n_m){
            for(s in 1:((n_iter / n_thin) - n_fail[c])){
              lik[s, i, c] <- dmvnorm(x = Yt_post[s, c(i, i + n_m)],
                                      mean = m_star_Xt[, s, i],
                                      sigma = S_star_Xt[, , s, i],
                                      log = F)
            }
          }
          if(n_fail[c] > 0){ 
            # removing NA rows of lik for chains with failed samples
            # (NA rows are the final rows of lik)
            # the calculate
            # sum_{i=1}^{n_g} log( (1/S) * sum_{s=1}^{S} p(Yt_i|theta_s,D_{-i}) )
            
            # so firstly, lik_c are the p(Yt_{i,s}|theta_s,D_{-i}) for chain c
            lik_c <- lik[-(((n_iter / n_thin) - n_fail[c] + 1):(n_iter / n_thin)), , c]
          }else{
            lik_c <- lik[, , c]    
          }
          # then the post_avg_lik_c are the (1/S) * sum_{s=1}^{S} lik_c[,i,c] for chain c
          # for chain c
          post_avg_lik_c <- apply(lik_c, 2, mean)
          # and the out.loo[c] are the sum_{i=1}^{n_g} log(post_avg_lik_c) for chain c
          out.loo[c] <- -2 * sum(log(post_avg_lik_c)) # multiply by -2 to get on deviance scale
        }
      }
    }else{ # assessing model fit using approximate LOO-CV (Pareto smoothed importance sampling)
      loglik <- array(NA, dim = c(n_iter / n_thin, n_m, n_chains))
      out.loo <- vector("list", length = n_chains)
      n_fail <- rep(NA, n_chains)
      failed_samp <- vector("list", n_chains)
      for(c in 1:n_chains){
        
        set.seed(seed = seed)
        alpha_post <- posteriors[[c]][, paste0("alpha[", 1:2, "]")]
        sigma_post <- posteriors[[c]][, paste0("sigma[", 1:2, "]")]
        tau_post <- posteriors[[c]][, paste0("tau[", 1:2, "]")]
        
        Yt_post <- posteriors[[c]][, paste0("Yt[", 1:(2*n_m), "]")]
        Xt_post <- posteriors[[c]][, paste0("Xt[", rep(1:n_m, dim.X),
                                            ",", rep(1:dim.X, rep(n_m, dim.X)),"]")]
        l_1_post <- posteriors[[c]][, paste0("l_1[", 1:dim.X, "]")]
        l_2_post <- posteriors[[c]][, paste0("l_2[", 1:dim.X, "]")]
        
        if(indep_resp == F){
          
          lambda_post <- posteriors[[c]][, c("lambda[1]", "lambda[2]")]
          rho_T_post <- posteriors[[c]][, "rho_T"]
          
        }
        
        m_Xt_train <- array(NA, dim=c(2 * n_m, (n_iter / n_thin)))
        m_Xt_test <- array(NA, dim=c(2 * n_m, (n_iter / n_thin)))
        T_1 <- array(NA, dim = c(2, 2, (n_iter / n_thin)))
        T_2 <- array(NA, dim = c(2, 2, (n_iter / n_thin)))
        S_Xt_train <- array(NA, dim=c(2 * n_m, 2 * n_m, (n_iter / n_thin)))
        K_Xt_test_Xt_train <- array(NA, dim=c(2 * n_m, 2 * n_m, (n_iter / n_thin)))
        S_Xt_test <- array(NA, dim=c(2 * n_m, 2 * n_m, (n_iter / n_thin)))
        m_star_Xt <- array(NA, dim=c(2 * n_m, (n_iter / n_thin)))
        S_star_Xt <- array(NA, dim=c(2 * n_m, 2 * n_m, (n_iter / n_thin)))
        
        Xt_train <- array(NA,dim = c(n_m, dim.X * n_m, (n_iter / n_thin)))
        Xt_test <- array(NA,dim = c(n_m, dim.X * n_m, (n_iter / n_thin)))
        
        for(s in 1:(n_iter / n_thin)){
          
          m_Xt_train[1:n_m, s] <- rep(alpha_post[s, 1], n_m)
          m_Xt_train[(n_m + 1):(2 * n_m), s] <- rep(alpha_post[s, 2], n_m)
          m_Xt_test[1:n_m, s] <- rep(alpha_post[s, 1], n_m)
          m_Xt_test[(n_m + 1):(2 * n_m), s] <- rep(alpha_post[s, 2], n_m)
          
          Xt_train[,,s] <- matrix(Xt_post[s,], n_m, dim.X * n_m)
          Xt_test[,,s] <- matrix(Xt_post[s,], n_m, dim.X * n_m)
          
            
          T_1[, , s] <- matrix(c(tau_post[s, 1]^2,
                                 tau_post[s, 1] * tau_post[s ,2] * rho_T_post[s],
                                 tau_post[s, 1] * tau_post[s, 2] * rho_T_post[s],
                                 tau_post[s, 2]^2), 2, 2)
        
          T_2[, , s] <- matrix(c(lambda_post[s, 1]^2 * tau_post[s, 1]^2,
                                 lambda_post[s, 1] * lambda_post[s, 2] * tau_post[s, 1] * tau_post[s ,2] * rho_T_post[s],
                                 lambda_post[s, 1] * lambda_post[s, 2] * tau_post[s, 1] * tau_post[s, 2] * rho_T_post[s],
                                 lambda_post[s, 2]^2 * tau_post[s, 2]^2), 2, 2)
           
          
          S_Xt_train[, , s] <- diag(c(rep(sigma_post[s, 1]^2,
                                          n_m),
                                      rep(sigma_post[s, 2]^2,
                                          n_m)),
                                    nrow = 2 * n_m) + kronecker(T_1[, , s],
                                                                calcSE3_ARD_square(X1 = Xt_train[, , s],
                                                                                   X2 = Xt_train[, , s],
                                                                                   l = l_1_post[s, ])) + kronecker(T_2[, , s],
                                                                                                                   calcSE3_ARD_square(X1 = Xt_train[, , s],
                                                                                                                                      X2 = Xt_train[, , s],
                                                                                                                                      l = l_2_post[s, ]))
          
          K_Xt_test_Xt_train[, , s] <- kronecker(T_1[, , s],
                                                 calcSE3_ARD_nonsquare(X1 = Xt_test[, , s],
                                                                       X2 = Xt_train[, , s],
                                                                       l = l_1_post[s, ])) + kronecker(T_2[, , s],
                                                                                                       calcSE3_ARD_nonsquare(X1 = Xt_test[, , s],
                                                                                                                             X2 = Xt_train[, , s],
                                                                                                                             l = l_2_post[s, ])) 
          
          S_Xt_test[, , s] <- diag(c(rep(sigma_post[s, 1]^2, n_m),
                                     rep(sigma_post[s, 2]^2, n_m)),
                                   nrow = 2 * n_m) + kronecker(T_1[, , s],
                                                               calcSE3_ARD_square(X1 = Xt_test[, , s],
                                                                                  X2 = Xt_test[, , s],
                                                                                  l = l_1_post[s, ])) + kronecker(T_2[, , s],
                                                                                                                  calcSE3_ARD_square(X1 = Xt_test[, , s],
                                                                                                                                     X2 = Xt_test[, , s],
                                                                                                                                     l = l_2_post[s, ]))
          
          m_star_Xt[, s] <- m_Xt_test[, s] + K_Xt_test_Xt_train[, , s] %*% solve(S_Xt_train[, , s]) %*% (t(Yt_post[s, ]) - m_Xt_train[, s])
          S_star_Xt[, , s] <- S_Xt_test[, , s] - K_Xt_test_Xt_train[, , s] %*% solve(S_Xt_train[, , s]) %*% t(K_Xt_test_Xt_train[, , s])
          S_star_Xt[2, 1, s] <- S_star_Xt[1, 2, s] 
        }
        negative_variance <- rep(NA,(n_iter / n_thin))
        for(s in 1:(n_iter / n_thin)){ #######
          negative_variance[s] <- any(S_star_Xt[, , s] < 0)
        }
        failed_samp[[c]] <- which(negative_variance == T)
        n_fail[c] <- length(failed_samp[[c]])
        
        if(n_fail[c] > 0){
          # removing posterior samples for which S_star_Xt<0
          S_star_Xt <- S_star_Xt[, , -failed_samp[[c]]]
          m_star_Xt <- m_star_Xt[,-failed_samp[[c]]]
          Yt_post <- Yt_post[-failed_samp[[c]], ]
        }
        if((n_iter/n_thin) - n_fail[c] < 30){
          # if a large majority of the samples fail (leaving less than 30 'good' samples),
          # just want to return an error of sorts
          for(i in 1:n_m){
            for(s in 1:(n_iter / n_thin)){
              loglik[s, i, c] <- NA
              out.loo[[c]] <- "Error: There are less than 30 samples for which S_star_Xt variances are positive in this chain, so the estimation of the LOO-CV-IC is aborted."
            }
          }
        }else{
          for(i in 1:n_m){
            for(s in 1:((n_iter / n_thin) - n_fail[c])){
              loglik[s, i, c] <- dmvnorm(x = Yt_post[s, c(i, i + n_m)],
                                         mean = m_star_Xt[c(i, i + n_m), s],
                                         sigma = matrix(c(S_star_Xt[i, i, s],
                                                          S_star_Xt[i, (i + n_m), s],
                                                          S_star_Xt[(i + n_m), i, s],
                                                          S_star_Xt[(i + n_m), (i + n_m), s]), 2, 2),
                                         log = T)
            }
          }
          if(n_fail[c] > 0){ 
            # removing NA rows of loglik for chains with failed samples
            # (NA rows are the final rows of loglik)
            out.loo[[c]] <- loo(loglik[-(((n_iter / n_thin) - n_fail[c] + 1):(n_iter / n_thin)), , c])
          }else{
            out.loo[[c]] <- loo(loglik[, , c])    
          }
        }
      }
    }
  }
  
  
  return(list(posteriors=posteriors,summary=summary_post,
              #X__joint=X__joint, X__marginal = X__marginal,Y_=Y_,mu_X=mu_X,scaleXt_0=scaleXt_0,
              #data_list = data_list,
              psrf=psrf,
              ess_chain=ess_chain,
              inits_list = inits_list,
              init.state = init.state
              #,ess_total=ess_total
              , out.loo = out.loo
              , X_new=X_new, n_fail=n_fail, failed_samp=failed_samp#,m_star_Xt=m_star_Xt
              ,n_iter=n_iter,n_thin=n_thin
              #, Yt_post=Yt_post, S_star_Xt=S_star_Xt
              ,obs_scaling=obs_scaling
              #, S_Xt_train=S_Xt_train,K_Xt_test_Xt_train=K_Xt_test_Xt_train,S_Xt_test=S_Xt_test,
              #K_ARD_train=K_ARD_train
              #m_Xt_test=m_Xt_test,m_Xt_train=m_Xt_train,Xt_train=Xt_train,Xt_test=Xt_test,T_=T_
  ))
}
