#####################################################################################
#                                                                                   #
#                                                                                   #
#                                       EXAMPLE                                     #
#                                                                                   #
#                                                                                   #
#####################################################################################

# parameters for the simulation
set.seed(123)
alpha <- 5
beta <- -0.23
lambda <- 0.5
k <- 2.5
p <- 0.5

# sample of real cases
n <- 500
t.real <- rep(1:10, n)[1:n]
y <- rpois(n = n,lambda)
z <- rnbinom(n = n, size = k, mu = exp(alpha + beta*t.real))
e <- (runif(n = n) < p)
x.real <- (1-e)*y + e*z 

# sample from an experiment of non-shooters
n <- 50
x.exp0 <- rpois(n = n,lambda)

# sample from an experiment of shooters
n <- 50
t.exp1 <- rep(1:10, n)[1:n]
x.exp1 <- rnbinom(n = n, size = k, mu = exp(alpha + beta*t.exp1))

# estimation based on all samples. 
est.parm <- estimate_model(x=x.real, t=t.real,x_exp1 = x.exp1, t_exp1 = t.exp1)

# est.parm contains the estimates of the parameters and two estimates
# of variances - using observed infromation and using numerical differentiation
# to obstined the estimates use:
est.parm$result

# for the variance based on Loius method:
est.parm$variance_matrix_result

# for the variance based on numerical differentation:
est.parm$inf_obs

# estimation based only on real data. 
est.parm <- estimate_model(x=x.real, t=t.real)

# confidence intervals for parameters. 
# The variance estimate and level of confidence can be changed
CI_for_theta(est.parm$result, list(est.parm$variance_matrix_result), 0.95)

# confidence intervals for the probability of being a shooter - log10 scale (two sided)
Calculate_CI_for_Probs(gunshot_values = 4,
                       time_after_shooting_values = 6,
                       is_shooter = TRUE,
                       result_matrix = est.parm$result,
                       variance_matrix_list = list(est.parm$variance_matrix_result), 
                       confidence_level = 0.95)

# confidence intervals for the likelihood ratio (one sided)
Calculate_CI_for_LR(gunshot_values = 4, 
                    time_after_shooting_values = 6, 
                    result_matrix = est.parm$result,
                    variance_matrix_list = list(est.parm$variance_matrix_result), 
                    confidence_level = 0.95)
