# code to produce simulation example, two outputs, two inputs, estimating their
# joint relationship using errors-in-variables multi-output GP, then optimising
# input variables given the estimated relationship and some desired output
+
# now use the posterior distribution for the hyperparameters and true values to fit the backward model

# need some desired output
Ystar <- c(18, 237)

sim_backward <- EIVMOGP_backward(posteriors = sim_forward$posteriors,
                                 obs_scaling = c(2, 3, 1, 1), Y_star = Ystar, n_m = 13,
                                 centre_X = c(7, 5), rho_X = c(6, 4),
                                 R = 2000, n_adapt = 500, n_burnin = 15000, mcmc_n = 20000,
                                 R_n = 10, n_thin = 10, initialisation = T)


# plot the joint posterior density of X* in this case, and compare with possible points in
# input space that could produce either desired output

# want to identify these possible points. will firstly values of the outputs throughout
# the input space
x1 <- seq(1, 13, length.out=1000)
x2 <- seq(1, 9, length.out=1000)

fun_2in_2out_GP <- function(x1, x2){
  0.3*(80 + 80*sin(x1) + exp(0.55*x2))
}

fun_2in_2out_GP_sim2 <- function(x1, x2){
  3*((x1-10)^2 - 0.9*(x2-5)^3 + 60)
}

y1 <- outer(x1, x2, fun_2in_2out_GP)
y2 <- outer(x1, x2, fun_2in_2out_GP_sim2)

# now identify those outputs that are within 1% of the desired output for each output variable,
# and store the corresponding input points
Xindex_Ystar1 <- which(0.99 * Ystar[1] < y1 & y1 < 1.01 * Ystar[1], arr.ind = T)
X_Ystar1 <- matrix(NA, nrow(Xindex_Ystar1), 2)
for(n in 1:nrow(Xindex_Ystar1)){
  
  X_Ystar1[n, 1] <- x1[Xindex_Ystar1[n, 1]] / 10
  X_Ystar1[n, 2] <- x2[Xindex_Ystar1[n, 2]] / 10
  
}
X_Ystar1 <- data.frame(Xstar1 = X_Ystar1[, 1],
                                 Xstar2 = X_Ystar1[, 2])

Xindex_Ystar2 <- which(0.99 * Ystar[2] < y2 & y2 < 1.01 * Ystar[2], arr.ind = T)
X <- matrix(NA, nrow(Xindex_Ystar2), 2)
for(n in 1:nrow(Xindex_Ystar2)){
  
  X[n, 1] <- x1[Xindex_Ystar2[n, 1]] / 10
  X[n, 2] <- x2[Xindex_Ystar2[n, 2]] / 10
  
}
X_Ystar2 <- data.frame(Xstar1 = X_Ystar2[, 1],
                                 Xstar2 = X_Ystar2[, 2])

Xindex_Ystar <- which(0.99 * Ystar[1] < y1 & y1 < 1.01 * Ystar[1] & 0.99 * Ystar[2] < y2 & y2 < 1.01 * Ystar[2], arr.ind = T)
X_Ystar <- matrix(NA, nrow(Xindex_Ystar), 2)
for(n in 1:nrow(Xindex_Ystar)){
  
  X_Ystar[n, 1] <- x1[Xindex_Ystar[n, 1]] / 10
  X_Ystar[n, 2] <- x2[Xindex_Ystar[n, 2]] / 10
  
}

X_Ystar <- data.frame(Xstar1 = X_Ystar[, 1],
                                Xstar2 = X_Ystar[, 2])

require(ggplot2)
ggplot(data = rbind(data.frame(X_Ystar1, group="Ystar1"),
                    data.frame(X_Ystar2, group="Ystar2")),
       aes(x = Xstar1, y = Xstar2))+
  geom_point(aes(color = group), alpha = 0.5) +
  geom_point(data = X_Ystar, aes(x = Xstar1, y = Xstar2),
             color = "white", alpha = 0.5) +
  geom_point(data = data.frame(X1 = Xt_eivgp_2cov[, 1] / 10,
                               X2 = Xt_eivgp_2cov[, 2] / 10),
             aes(x = X1, y = X2), color = "black")


# now reproduce plot, with posterior density of X* superimposed

joint_X_star_density <- kde2d(sim_backward$postXstar[, 1],
                              sim_backward$postXstar[, 2],  n = 1001) # joint posterior density of X*
Xstar1_mode_index <- which(joint_X_star_density$z == max(joint_X_star_density$z), arr.ind = T)[1, 1]
Xstar2_mode_index <- which(joint_X_star_density$z == max(joint_X_star_density$z), arr.ind = T)[1, 2]
Xstar1_mode <- joint_X_star_density$x[Xstar1_mode_index]
Xstar2_mode <- joint_X_star_density$y[Xstar2_mode_index]
#[1] 0.5229141 0.6711764 is the global mode

Xstar_post <- data.frame(Xstar1 = sim_backward$postXstar[, 1],
                         Xstar2 = sim_backward$postXstar[, 2])

Xstar_mode2_index <- which(joint_X_star_density$z[-(1:500), ] == max(joint_X_star_density$z[-(1:500), ]), arr.ind = T)
Xstar_mode2 <- c(joint_X_star_density$x[Xstar_mode2_index[, 1] + 500], joint_X_star_density$y[Xstar_mode2_index[, 2]])
#[1] 0.9738012 0.2279954 is a local mode

ggplot(data = rbind(data.frame(X_Ystar1, group="Ystar1"),
                    data.frame(X_Ystar2, group="Ystar2")),
       aes(x = Xstar1, y = Xstar2))+
  geom_point(aes(color = group), alpha = 0.5)+
  stat_density_2d(data = Xstar_post, 
                  aes(fill = ..level..), geom = "polygon", n = 400, alpha = 0.25) +
  geom_point(aes(x = Xstar1_mode, y = Xstar2_mode), colour = "pink")+
  geom_point(aes(x = Xstar_mode2[1], y = Xstar_mode2[2]), colour = "yellow")