# function for running the backward model, i.e., optimising the input variables
# given the estimated relationship from forward model and some desired multivariate
# output


EIVMOGP_backward <- function(posteriors, R = 2000, R_n = 10, mcmc_n = 20000, Y_star,
                               dim.X = 2, dim.Y = 2, rho_X = rep(1,dim.X), centre_X = rep(0,dim.X), n_m = 7,
                               obs_scaling, LHS = T, n_thin = 10, n_burnin = 15000, n_adapt = 500,
                               initialisation = T){
  # posteriors: posterior samples of parameters estimated in forward modelling
  # R: number of backward models to be fitted (number of posterior samples from forward model to use)
  # R_n: number of posterior samples of X* to take from each backward model
  # mcmc_n: number of samples drawn from posterior X* (in each backward model)
  # Y_star: desired output vector
  # dim.X: number of inputs
  # dim.Y: number of outputs
  # rho_X: range for uniform prior on X* (either side lengths for box prior, or radii for disk prior)
  # centre_X: centre for uniform prior on X*
  # n_m: number of groups
  # obs_scaling: scaling factors for each variable of observed data
  # LHS: whether to use a Latin hypercube sampling for choosing initial values (if initialisation == T)
  # n_thin: thinning parameter for MCMC
  # n_burnin: number of iterations for burn-in of MCMC algorithm
  # n_adapt: number of iterations for adaptation
  # initialisation: whether or not to choose the starting values of X* for the MCMC algorithm
  
  for(d_Y in 1:dim.Y){
    
    Y_star[d_Y] <- Y_star[d_Y] / 10^obs_scaling[d_Y] # scaling the desired output
    
  }
  
  for(d_X in 1:dim.X){
    rho_X[d_X] <- rho_X[d_X] / 10^obs_scaling[d_X + 2]
    centre_X[d_X] <- centre_X[d_X] / 10^obs_scaling[d_X + 2]  # scaling the uniform prior ranges and centres
  }
  
  require(rjags)
  set.seed(999)
  if(is.list(posteriors)){
    posteriors <- posteriors[[1]]
  }
  if(R<nrow(posteriors)){
    samp <- sample(nrow(posteriors), R, replace = F)
    Xt.samples.R <- as.mcmc(posteriors[samp, paste0("Xt[", 1:n_m, ",", rep(1:dim.X, rep(n_m, dim.X)), "]")])
    Yt.samples.R <- as.mcmc(posteriors[samp, paste0("Yt[", 1:(2*n_m), "]")])
    alpha.samples.R <- as.mcmc(posteriors[samp, paste0("alpha[", 1:dim.Y, "]")])
    sigma.samples.R <- as.mcmc(posteriors[samp, paste0("sigma[", 1:dim.Y, "]")])
    tau.samples.R <- as.mcmc(posteriors[samp, paste0("tau[", 1:dim.Y, "]")])
    rho_T.samples.R <- as.mcmc(posteriors[samp, "rho_T"])
    lambda.samples.R <- as.mcmc(posteriors[samp, paste0("lambda[", 1:dim.Y, "]")])
    l.samples.R <- as.mcmc(posteriors[samp,c(paste0("l_1[", 1:dim.X, "]"), paste0("l_2[", 1:dim.X, "]"))])
  }else{
    Xt.samples.R <- as.mcmc(posteriors[, paste0("Xt[", 1:n_m, ",", rep(1:dim.X, rep(n_m, dim.X)), "]")])
    Yt.samples.R <- as.mcmc(posteriors[, paste0("Yt[", 1:n_m, ",", rep(1:dim.Y, rep(n_m, dim.Y)), "]")])
    alpha.samples.R <- as.mcmc(posteriors[, paste0("alpha[", 1:dim.Y, "]")])
    sigma.samples.R <- as.mcmc(posteriors[, paste0("sigma[", 1:dim.Y, "]")])
    tau.samples.R <- as.mcmc(posteriors[, paste0("tau[", 1:dim.Y, "]")])
    rho_T.samples.R <- as.mcmc(posteriors[, "rho_T"])
    lambda.samples.R <- as.mcmc(posteriors[, paste0("lambda[", 1:dim.Y, "]")])
    l.samples.R <- as.mcmc(posteriors[samp,c(paste0("l_1[", 1:dim.X, "]"), paste0("l_2[", 1:dim.X, "]"))])
  }
  if(initialisation){ # choose randomised starting values for X* for each backward model
    if(LHS){
      require(lhs)
      set.seed(999)
      lhs_Xstar <- randomLHS(n = 2*R, k = dim.X)
      X_star_initial <- matrix(NA, nrow = 2 * R, ncol=dim.X)
      for(p in 1:dim.X){
        X_star_initial[,p] <- qunif(lhs_Xstar[,p], min = centre_X[p]-rho_X[p], max = centre_X[p]+rho_X[p])
      }
    }
  }else{
    
    X_star_initial <- matrix(NA, nrow = 2 * R, ncol=dim.X)
    for(n in 1:(2*R)){
      X_star_initial[n, ] <- centre_X
    }
    
  }
  
  
  I <- diag(n_m)
  X_star <- paste0("X_star[",1:dim.X,"]")
  models <- vector("list",length = R)          # storing of samples and MCMC diagnostics
  posteriors <- vector("list",length = R)
  posteriors_all <- vector("list", length = R)
  psrf <- vector("list", length = R)
  ess_total <- vector("list", length = R)
  ess_chain <- vector("list", length = R)
  
  for(r in 1:R){
    print(paste0("Backward model iteration ", r, " out of ", R))
    Xt_draw <- matrix(Xt.samples.R[r, ], n_m, dim.X)
    Yt_draw <- as.vector(Yt.samples.R[r, ])
    alpha_draw <- as.vector(alpha.samples.R[r, ])
    sigma_draw <- as.vector(sigma.samples.R[r, ])
    tau_draw <- as.vector(tau.samples.R[r, ])
    lambda_draw <- as.vector(lambda.samples.R[r, ])
    rho_T_draw <- rho_T.samples.R[r]
    l_draw <- matrix(l.samples.R[r, ], dim.Y, dim.X, byrow = T)
    data.list.back <- list('Y_star' = Y_star, 'Xt' = Xt_draw, 'Yt' = Yt_draw, 'alpha' = alpha_draw, 'n_m' = n_m,
                           'sigma' = sigma_draw, 'l' = l_draw, 'rho_X' = rho_X, 'centre_X' = centre_X, 'I' = I,
                           'tau' = tau_draw, 'dim.X' = dim.X, 'rho_T' = rho_T_draw, 'lambda' = lambda_draw,
                           'dim.Y' = dim.Y)
    models[[r]] <- try(jags.model('MOGPEIVSEARD_back.bug', data = data.list.back, n.chains = 2, n.adapt = n_adapt,
                                  inits = list(list('X_star' = X_star_initial[(2*r-1), ],
                                                    .RNG.name="base::Wichmann-Hill",.RNG.seed = 999),
                                               list('X_star' = X_star_initial[2*r, ],
                                                    .RNG.name="base::Wichmann-Hill",.RNG.seed = 1000))
    ))
    try(update(models[[r]], n_burnin))
    posteriors[[r]] <- try(coda.samples(models[[r]], 'X_star',
                                        #c('X_star','m_star','S_star'),
                                        mcmc_n, thin = n_thin))
    psrf[[r]] <- try(gelman.diag(mcmc.list(posteriors[[r]][[1]],
                                           posteriors[[r]][[2]])))
    ess_total[[r]] <- try(effectiveSize(mcmc.list(posteriors[[r]][[1]],
                                                  posteriors[[r]][[2]])))
    ess_chain[[r]] <- try(lapply(posteriors[[r]][, X_star], effectiveSize))
    samp.r <- sample(2 * (mcmc_n / n_thin), R_n, replace = F)
    posteriors_all[[r]] <- try(as.matrix(rbind(posteriors[[r]][[1]],
                                               posteriors[[r]][[2]])[samp.r,]))
    
  }
  
  failedsamp <- which(lapply(posteriors_all, is.matrix) == F) # tracking if models are not able to be fitted
  legitsamp <- which(lapply(posteriors_all, is.matrix) == T)
  if(length(failedsamp)==0){
    postXstar <- matrix(NA, R*R_n, dim.X)
    for(r in 1:R){
      
      # each element in list of posteriors_all should be a matrix with dim.X + dim.Y + dim.Y^dim.Y columns
      # since there are X_star has dim.X dimensions, m_star has dim.Y dimensions, and S_star has dim.Y^dim.Y dimensions
      # Let's assume dim.Y is 2, so we have dim.X + 6 columns in each element of posteriors_all list
      # the number of rows of each element in posteriors_all is 10(=R_n), as we take a random sample of R_n from each MCMC run
      # we want to take the first dim.X columns of each element of posteriors_all
      postXstar[(R_n*r-(R_n-1)):(R_n*r), ] <- posteriors_all[[r]]
      #postXstar[(R_n*r-(R_n-1)):(R_n*r),d] <- unlist(posteriors_all)[((dim.X+2)*R_n*r-((dim.X+1)*R_n-1)+(d-1)*R_n):((dim.X+2)*R_n*r-dim.X*R_n+(d-1)*R_n)]
      
    }
    
  }else{
    postXstar <- matrix(NA,(R-length(failedsamp))*R_n, dim.X)
    for(r in 1:(R-length(failedsamp))){
      
      postXstar[(R_n*r-(R_n-1)):(R_n*r), ] <- posteriors_all[[legitsamp[r]]]
      #postXstar[(R_n*i-(R_n-1)):(R_n*i),d] <- unlist(posteriors_all[-failedsamp])[((dim.X+2)*R_n*i-((dim.X+1)*R_n-1)+(d-1)*R_n):((dim.X+2)*R_n*i-dim.X*R_n+(d-1)*R_n)]
      
    }
    
  }
  return(list(posteriors_all=posteriors_all, psrf = psrf, ess_chain = ess_chain, ess_total = ess_total,
              samp=samp,models=models,
              posteriors=posteriors,postXstar=postXstar,failedsamp=failedsamp,
              n_failedsamp=length(failedsamp), X_star_initial = X_star_initial))
  
}
