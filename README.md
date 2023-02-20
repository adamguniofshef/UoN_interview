# UoN_interview

# In order to run code, JAGS must be installed on your computer. You can install the latest version here: https://sourceforge.net/projects/mcmc-jags/

# EIVMOGP_forward.R -- this script details the function EIVMOGP_forward, which fits the EIV MOGP forward model (estimating the true values for each group, and the hyperparameters of the GP), with two further functions defined that are used within EIVMOGP_forward

# EIVMOGP_backward.R -- this script details the function EIVMOGP_backward, which fits the EIV MOGP backward model (algorithm for optimising the input variables given the posterior distribution of the hyperparameters and some multivariate desired output)

# EIVMOGP_simulation.R -- this script runs the simulation which is discussed in the presentation

# MOGPEIVSEARD.bug -- this script is a .bug file, which defines the model that is accessed using rjags. This file corresponds to the forward model.

# MOGPEIVSEARD_back.bug -- this script is a .bug file, which defines the model that is accessed using rjags. This file corresponds to the backward model.
