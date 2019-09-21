library(bbmle)
library (stats4)
library(numDeriv)
library(MASS)
library(ggplot2)
library(gridExtra)
library(grid)

set.seed(123)

#################################################################################
#                                                                               #
#                                                                               #
#                                 HELPING METHODS                               #
#                                                                               #
#                                                                               #
#################################################################################


## compute the estimated gradient vector for a spacific value
Grad <- function(x, t, p, lambda, k, alpha, beta){
  w_i <- (p*dnbinom(x, k, mu = exp(alpha + beta*t)))/
    (p*dnbinom(x, k, mu = exp(alpha + beta*t)) + (1-p)*dpois(x, lambda))
  grad_1 <- (w_i/p)-((1-w_i)/(1-p))
  grad_2 <- (1-w_i)*((x/lambda)-1)
  grad_3 <- w_i*(digamma(x+k) - digamma(k) + log(k) + 1 - log(k + exp(alpha + beta*t))
                 - k/(k + exp(alpha + beta*t)) - x/(k + exp(alpha + beta*t)))
  grad_4 <- w_i*k*((x - exp(alpha+beta*t))/(k + exp(alpha + beta*t)))
  grad_5 <- t*grad_4
  return(c(grad_1, grad_2, grad_3, grad_4, grad_5))
}

## compute the estimated negative hessian for a specific value
Hess <- function(x,t, p, lambda, k, alpha, beta){
  w_i <- (p*dnbinom(x, k, mu = exp(alpha + beta*t)))/
    (p*dnbinom(x, k, mu = exp(alpha + beta*t)) + (1-p)*dpois(x, lambda))
  row_1 <- c((w_i/p^2)+((1-w_i)/(1-p)^2),
             0,
             0,
             0,
             0)
  row_2 <- c(0,
             (1-w_i)*(x/lambda^2),
             0,
             0,
             0)
  row_3 <- c(0,
             0,
             w_i*(-trigamma(x+k) + trigamma(k) - 1/k + 
                    (2*exp(alpha + beta*t)+k-x)/((k + exp(alpha + beta*t))^2)),
             w_i*(exp(alpha + beta*t)*(exp(alpha + beta*t) - x)/(k + exp(alpha + beta*t))^2),
             w_i*t*(exp(alpha + beta*t)*(exp(alpha + beta*t) - x)/(k + exp(alpha + beta*t))^2))
  row_4 <- c(0,
             0,
             w_i*(exp(alpha + beta*t)*(exp(alpha + beta*t) - x)/(k + exp(alpha + beta*t))^2),
             w_i*k*exp(alpha + beta*t)*(x + k)/(k + exp(alpha + beta*t))^2,
             w_i*t*k*exp(alpha + beta*t)*(x + k)/(k + exp(alpha + beta*t))^2)
  row_5 <- c(0,
             0,
             w_i*t*(exp(alpha + beta*t)*(exp(alpha + beta*t) - x)/(k + exp(alpha + beta*t))^2),
             w_i*t*k*exp(alpha + beta*t)*(x + k)/(k + exp(alpha + beta*t))^2,
             w_i*t^2*k*exp(alpha + beta*t)*(x + k)/(k + exp(alpha + beta*t))^2)
  return(rbind(row_1, row_2, row_3, row_4, row_5))
}

## diag of the mix multiplicative gradient (when i=j)
Mix_Mat_diag <- function(x,t, p, lambda, k, alpha, beta){
  w_i <- (p*dnbinom(x, k, mu = exp(alpha + beta*t)))/
    (p*dnbinom(x, k, mu = exp(alpha + beta*t)) + (1-p)*dpois(x, lambda))
  col_1 <- c((w_i/p^2)+((1-w_i)/(1-p)^2),
             -((1-w_i)/(1-p))*((x/lambda)-1),
             (w_i/p)*(digamma(x+k) - digamma(k) + log(k) + 1 - log(k + exp(alpha + beta*t))
                      - k/(k + exp(alpha + beta*t)) - x/(k + exp(alpha + beta*t))),
             (w_i/p)*k*((exp(alpha + beta*t) - x)/(k + exp(alpha + beta*t))),
             (w_i*t/p)*k*((exp(alpha + beta*t) - x)/(k + exp(alpha + beta*t))))
  col_2 <- c(-((1-w_i)/(1-p))*((x/lambda)-1),
             (1-w_i)*((x/lambda)-1)^2,
             0,
             0,
             0)
  col_3 <- c((w_i/p)*(digamma(x+k) - digamma(k) + log(k) + 1 - log(k + exp(alpha + beta*t))
                      - k/(k + exp(alpha + beta*t)) - x/(k + exp(alpha + beta*t))),
             0,
             w_i*(digamma(x+k) - digamma(k) + log(k) + 1 - log(k + exp(alpha + beta*t))
                  - k/(k + exp(alpha + beta*t)) - x/(k + exp(alpha + beta*t)))^2,
             w_i*k*(digamma(x+k) - digamma(k) + log(k) + 1 - log(k + exp(alpha + beta*t))
                    - k/(k + exp(alpha + beta*t)) - x/(k + exp(alpha + beta*t)))*
               ((x - exp(alpha+beta*t))/(k + exp(alpha + beta*t))), 
             w_i*k*t*(digamma(x+k) - digamma(k) + log(k) + 1 - log(k + exp(alpha + beta*t))
                      - k/(k + exp(alpha + beta*t)) - x/(k + exp(alpha + beta*t)))*
               ((x - exp(alpha+beta*t))/(k + exp(alpha + beta*t))))
  col_4 <- c((w_i/p)*k*((exp(alpha + beta*t) - x)/(k + exp(alpha + beta*t))),
             0,
             w_i*k*(digamma(x+k) - digamma(k) + log(k) + 1 - log(k + exp(alpha + beta*t))
                    - k/(k + exp(alpha + beta*t)) - x/(k + exp(alpha + beta*t)))*
               ((x - exp(alpha+beta*t))/(k + exp(alpha + beta*t))),
             w_i*k^2*((x - exp(alpha+beta*t))/(k + exp(alpha + beta*t)))^2,
             w_i*t*k^2*((x - exp(alpha+beta*t))/(k + exp(alpha + beta*t)))^2)
  col_5 <- c((w_i*t/p)*k*((exp(alpha + beta*t) - x)/(k + exp(alpha + beta*t))),
             0,
             w_i*k*t*(digamma(x+k) - digamma(k) + log(k) + 1 - log(k + exp(alpha + beta*t))
                      - k/(k + exp(alpha + beta*t)) - x/(k + exp(alpha + beta*t)))*
               ((x - exp(alpha+beta*t))/(k + exp(alpha + beta*t))),
             w_i*t*k^2*((x - exp(alpha+beta*t))/(k + exp(alpha + beta*t)))^2,
             w_i*t^2*k^2*((x - exp(alpha+beta*t))/(k + exp(alpha + beta*t)))^2)
  return (rbind(col_1, col_2, col_3, col_4, col_5))
}

## compute the estimated variance with the information matrix
Compute_Variance <- function(information_matrix){
  return(solve(information_matrix))
}

## compute the likelihood ratio between shooters (numerator) and non-shooters (denominator)
LR <- function(x, t, lambda, alpha, beta, k){
  dnbinom(x, k, mu = exp(alpha + beta*t))/dpois(x, lambda)
}

## definite the log-likelihood function of the mixtue of poisson and neg-binomial
## the function take (k,p,lambda,alpha,beta) as input parameters
Log.Like <- function (k,p,lambda,alpha,beta){
  -sum(log((1-p)*dpois(x,lambda) + p*dnbinom(x,k,mu = exp(alpha + beta*t))))
}

## definite the log-likelihood function of the mixture of poisson and neg-binomial
## the function take (k,p,lambda,alpha,beta) as input parameters
Log.Like.1 <- function (theta){
  sum(log((1-theta[1])*dpois(x,theta[2]) + theta[1]*dnbinom(x,theta[3],mu = exp(theta[4] + theta[5]*t))))
}


###################################################################################################
#                                                                                                 #
#                                                                                                 #
#                         FUNCTIONS TO CALCULATE METHODS PROPOSED IN THE ARTICLE                  #  
#                                                                                                 #
#                                                                                                 #
###################################################################################################



## Calculate LR varinace for a given x (single observation) and a given t through proposed method by Mandel
## Parameters: x: Unique obs corresponding to a number of GSR found on  a suspect
##             t: Unique obs corresponding to the time after shooting that the corresponded x was measured
##             k: k param of the model
##             p: p param of the model
##             lambda: lambda param of the model
##             alpha: alpha param of the model
##             beta: beta param of the model
##             var_matrix: Variance matrix or Estimated variance matrix
## return: The variance of the calculate LR for x and t
LR_var_simulation <- function(x, t, k, p, lambda, alpha, beta, var_matrix){
  ## simulate 1000 random vectors from normal dist with expectation of our parameters vector and var matrix computed 
  theta_random <- mvrnorm(n = 1000, mu = c(p, lambda, k, alpha, beta), Sigma = var_matrix)
  ## calculate the LR of each vector simulated
  LR_values <- rep(0, 1000)
  for (i in 1:dim(theta_random)[1]){
    LR_values[i] <- log10(LR(x = x, t = t, lambda = theta_random[i, 2], alpha = theta_random[i, 4],
                             beta = theta_random[i, 5], k = theta_random[i, 3]))
  }
  ## calculated the var of the obtaind LR computed
  var_LR <- var(LR_values, na.rm = T)
  return(var_LR)
}

## Calculate varinace of probabilities for a given x (single observation) and a given t through proposed method by Mandel
## Parameters: x: Unique obs corresponding to a number of GSR found on  a suspect
##             t: Unique obs corresponding to the time after shooting that the corresponded x was measured
##             is_shooter: Boolean. If TRUE, variance for shooter probabilities is computed else for non-shooter
##             k: k param of the model
##             p: p param of the model
##             lambda: lambda param of the model
##             alpha: alpha param of the model
##             beta: beta param of the model
##             var_matrix: Variance matrix or Estimated variance matrix
## return: The variance of the calculate LR for x and t
Probs_var_simulation <- function(x, t, is_shooter,k, p, lambda, alpha, beta, var_matrix){
  ## simulate 1000 random vectors from normal dist with expectation of our parameters vector and var matrix computed 
  theta_random <- mvrnorm(n = 1000, mu = c(p, lambda, k, alpha, beta), Sigma = var_matrix)
  ## calculate the probs of each vector simulated
  probs_values <- rep(0, 1000)
  for (i in 1:dim(theta_random)[1]){
    if (is_shooter){
      probs_values[i] <- log10(dnbinom(x, theta_random[i, 3], mu = exp(theta_random[i, 4] + theta_random[i, 5]*t)))
    } else{
      probs_values[i] <- log10(dpois(x, theta_random[i, 2]))
    } 
  }
  ## calculated the var of the obtaind LR computed
  var_probs <- var(probs_values, na.rm = T)
  return(var_probs)
}

## Calculate two sided confidence interval for the vector of estimated parameters
## Parameters: result_matrix: A matrix that represent the estimated parameters
##             list_of_obs_mat: A list of Variance matrices corresponding to each row of result_matrix
##             confidence_level: the confidence level corresponding to 1-alpha
## return: list of coresponding CI's corresponding to each parameter of the model
CI_for_theta <- function(result_matrix, list_of_obs_mat, significance_level){
  CI_for_k <- matrix(0, ncol = 2, nrow = dim(result_matrix)[1]) 
  CI_for_p <- matrix(0, ncol = 2, nrow = dim(result_matrix)[1])
  CI_for_lambda <- matrix(0, ncol = 2, nrow = dim(result_matrix)[1])
  CI_for_alpha <- matrix(0, ncol = 2, nrow = dim(result_matrix)[1])
  CI_for_beta <- matrix(0, ncol = 2, nrow = dim(result_matrix)[1])
  colnames(CI_for_k) <- c('Lower Bound', 'Upper Bound')
  colnames(CI_for_p) <- c('Lower Bound', 'Upper Bound')
  colnames(CI_for_lambda) <- c('Lower Bound', 'Upper Bound') 
  colnames(CI_for_alpha) <- c('Lower Bound', 'Upper Bound')
  colnames(CI_for_beta) <- c('Lower Bound', 'Upper Bound')
  for (i in 1:dim(result_matrix)[1]){
    CI_for_k[i, 1] <- result_matrix[i, 1] - qnorm(p = 1 - significance_level/2)*sqrt(list_of_obs_mat[[i]][3,3])
    CI_for_p[i, 1] <- result_matrix[i, 2] - qnorm(p = 1 - significance_level/2)*sqrt(list_of_obs_mat[[i]][1,1])
    CI_for_lambda[i, 1] <- result_matrix[i, 3] - qnorm(p = 1 - significance_level/2)*sqrt(list_of_obs_mat[[i]][2,2])
    CI_for_alpha[i, 1] <- result_matrix[i, 4] - qnorm(p = 1 - significance_level/2)*sqrt(list_of_obs_mat[[i]][4,4])
    CI_for_beta[i, 1] <- result_matrix[i, 5] - qnorm(p = 1 - significance_level/2)*sqrt(list_of_obs_mat[[i]][5,5])
    CI_for_k[i, 2] <- result_matrix[i, 1] + qnorm(p = 1 - significance_level/2)*sqrt(list_of_obs_mat[[i]][3,3])
    CI_for_p[i, 2] <- result_matrix[i, 2] + qnorm(p = 1 - significance_level/2)*sqrt(list_of_obs_mat[[i]][1,1])
    CI_for_lambda[i, 2] <- result_matrix[i, 3] + qnorm(p = 1 - significance_level/2)*sqrt(list_of_obs_mat[[i]][2,2])
    CI_for_alpha[i, 2] <- result_matrix[i, 4] + qnorm(p = 1 - significance_level/2)*sqrt(list_of_obs_mat[[i]][4,4])
    CI_for_beta[i, 2] <- result_matrix[i, 5] + qnorm(p = 1 - significance_level/2)*sqrt(list_of_obs_mat[[i]][5,5])
  }
  return(list(CI_for_p = CI_for_p, CI_for_lambda = CI_for_lambda, CI_for_k = CI_for_k,
              CI_for_alpha = CI_for_alpha, CI_for_beta = CI_for_beta))
}

## Calculate one-sided confidence interval for the log(Probs) (basis 10) given the result matrix for the parameters 
## and the estimate variance matrix
## Parameters: gunshot_values: vector of values corresponding to number of GSR for which to calculate LR
##             time_after_shooting_values: vector of values corresponding to time after shooting hours for which
##                                         to calculate LR
##             is_shooter: Boolean. If TRUE, variance for shooter probabilities is computed else for non-shooter
##             result_matrix: A matrix that represent the estimated parameters
##             variance_matrix_list: A list of Variance matrices corresponding to each row of result_matrix
##             confidence_level: the confidence level corresponding to 1-alpha
## return: list. Each attribute in the list is a matrix coresponding to a bound of the 2 way CI's,
##          when the computed CI is of log(Probs) for each GSR and time after shooting
Calculate_CI_for_Probs <- function(gunshot_values, time_after_shooting_values, is_shooter,result_matrix,
                                variance_matrix_list, confidence_level){
  res <- list()
  for (i in 1:dim(result_matrix)[1]){
    if (sum(is.na(variance_matrix_list[[i]])) == 0){
      CI_matrix_lower <- matrix(0, ncol = length(time_after_shooting_values), nrow = length(gunshot_values))
      CI_matrix_upper <- matrix(0, ncol = length(time_after_shooting_values), nrow = length(gunshot_values))
      colnames(CI_matrix_lower) <- time_after_shooting_values
      rownames(CI_matrix_lower) <- gunshot_values
      colnames(CI_matrix_upper) <- time_after_shooting_values
      rownames(CI_matrix_upper) <- gunshot_values
      for (gsr in 1:length(gunshot_values)){
        for (t in 1:length(time_after_shooting_values)){
          if (is_shooter){
            probs <- log10(dnbinom(gunshot_values[gsr], result_matrix[i, 1], mu = exp(result_matrix[i, 4] + result_matrix[i, 5]*time_after_shooting_values[t])))
          } else {
            probs <- log10(dpois(gunshot_values[gsr], result_matrix[i, 3]))
          }
          probs_var <- Probs_var_simulation(gunshot_values[gsr], time_after_shooting_values[t], is_shooter,
                                            k = result_matrix[i, 1], p = result_matrix[i, 2],
                                            lambda = result_matrix[i, 3], alpha = result_matrix[i, 4],
                                            beta = result_matrix[i, 5], var_matrix = variance_matrix_list[[i]])
          CI_matrix_lower[gsr,t] <- probs - sqrt(probs_var)*qnorm(p = confidence_level+ (1-confidence_level)/2)
          CI_matrix_upper[gsr,t] <- probs + sqrt(probs_var)*qnorm(p = confidence_level+ (1-confidence_level)/2)
        }
      } 
      res[[i]] <- list(lower_bounds = CI_matrix_lower, upper_bound = CI_matrix_upper) 
    }
  }
  return(res)
}


## Calculate one-sided confidence interval for the log(LR) (basis 10) given the result matrix for the parameters and the
## estimate variance matrix
## Parameters: gunshot_values: vector of values corresponding to number of GSR for which to calculate LR
##             time_after_shooting_values: vector of values corresponding to time after shooting hours for which
##                                         to calculate LR
##             result_matrix: A matrix that represent the estimated parameters
##             variance_matrix_list: A list of Variance matrices corresponding to each row of result_matrix
##             confidence_level: the confidence level corresponding to 1-alpha
## return: list. Each attribute in the list is a matrix coresponding to one way CI's,
##          when the computed CI is the lower bound of log(LR) for each GSR and time after shooting
Calculate_CI_for_LR <- function(gunshot_values, time_after_shooting_values, result_matrix,
                                variance_matrix_list, confidence_level){
  res <- list()
  for (i in 1:dim(result_matrix)[1]){
    if (sum(is.na(variance_matrix_list[[i]])) == 0){
      CI_matrix_lower <- matrix(0, ncol = length(time_after_shooting_values), nrow = length(gunshot_values))
      colnames(CI_matrix_lower) <- time_after_shooting_values
      rownames(CI_matrix_lower) <- gunshot_values
      for (gsr in 1:length(gunshot_values)){
        for (t in 1:length(time_after_shooting_values)){
          LR_value <- LR(gunshot_values[gsr], time_after_shooting_values[t], k = result_matrix[i, 1],
                         lambda = result_matrix[i, 3], alpha = result_matrix[i, 4], beta = result_matrix[i, 5])
          LR_var <- LR_var_simulation(gunshot_values[gsr], time_after_shooting_values[t], k = result_matrix[i, 1],
                                      p = result_matrix[i, 2], lambda = result_matrix[i, 3],
                                      alpha = result_matrix[i, 4], beta = result_matrix[i, 5],
                                      var_matrix = variance_matrix_list[[i]])
          CI_matrix_lower[gsr,t] <- log10(LR_value) - sqrt(LR_var)*qnorm(p = confidence_level)
        }
      } 
      res[[i]] <- list(lower_bounds = CI_matrix_lower) 
    }
  }
  return(res)
}

## Function that estimated the GSR parameters model
## Parameters: x: Vector corresponding to a number of GSR found on  a suspect
##             t: vector corresponding to the time after shooting that the corresponded x was measured
##             experimental_data: Boolean. TRUE if experimental data is added to the model. Set by default to FALSE.
##                                IF TRUE is_shooter, x_exp and t_exp need to be provided.
##             x_exp1: Vector corresponding to a number of GSR found on  a suspect from the shooter experiment
##             t_exp1: vector corresponding to the time after shooting that the corresponded x was measured from the shooter experiment
##             x_exp0: Vector corresponding to a number of GSR found on  a suspect from the non shooter experiment
##             start_parameters: Vector corresponding to the start parameters for optimisation.
##                               The parameters are for (k, p, lambda, alpha, beta)
##             max_iter: Number of iteration to perform. Set by default to 80.
## return: list. result for the parameters estimation and variance_matrix_result for the corresponding
##         Variance matrix by Louis method and numerical derivative (inf_obs attribute)
estimate_model <- function(x, t, x_exp1=NULL, t_exp1=NULL, x_exp0=NULL,
                           start_parameters=c(0.5, 0.5, 0.001, 2, -1), max_iter=80){
  if (length(x) != length(t)){
    stop('Length of x and t should be equal')
  }
  if (length(x_exp1) != length(t_exp1)){
    stop('Length of x_exp1 and t_exp1 should be equal')
  }
  oldw <- getOption("warn")
  options(warn = -1)
  n <- length(x)
  k.start <- start_parameters[1]
  par.start <- start_parameters[2:length(start_parameters)]
  diff.par <- 1
  diff.k <- 1
  iter <- 1
  param.matrix <- NULL
  k.vector <- NULL
  p.opt <- NULL
  param.convergence <- NULL
  try(while (((diff.par > 10^-3) || (diff.k > 10^-3)) && (iter <= max_iter)){
    if (iter == 1){
      ## estimate p with EM style
      p.new <- sum((dnbinom(x, size = k.start, mu = exp(par.start[3]+ par.start[4]*t))*par.start[1])/
                     ((dnbinom(x, size = k.start, mu = exp(par.start[3]+ par.start[4]*t))*par.start[1]) +
                        (1-par.start[1])*dpois(x, par.start[2])))/n
      if (is.null(x_exp1) == FALSE & is.null(x_exp0) == FALSE){
        Log.Like.par <- function (lambda,alpha,beta){
          -sum(log((1-p.new)*dpois(x,lambda) + p.new*dnbinom(x,k.start,mu = exp(alpha + beta*t)))) -
                 sum(log(dnbinom(x_exp1, k.start, mu = exp(alpha + beta*t_exp1)))) - sum(log(dpois(x_exp0, lambda)))
        }
      } else if (is.null(x_exp1) == FALSE & is.null(x_exp0) == TRUE) { # experimental data on shooters
          Log.Like.par <- function (lambda,alpha,beta){
            -sum(log((1-p.new)*dpois(x,lambda) + p.new*dnbinom(x,k.start,mu = exp(alpha + beta*t)))) -
              sum(log(dnbinom(x_exp1, k.start, mu = exp(alpha + beta*t_exp1))))
          }
      } else if (is.null(x_exp1) == TRUE & is.null(x_exp0) == FALSE) { # experimental data on non shooters
        Log.Like.par <- function (lambda,alpha,beta){
          -sum(log((1-p.new)*dpois(x,lambda) + p.new*dnbinom(x,k.start,mu = exp(alpha + beta*t)))) -
            sum(log(dpois(x_exp0, lambda)))
        }
      } else { # no experimental data
        Log.Like.par <- function (lambda,alpha,beta){
          -sum(log((1-p.new)*dpois(x,lambda) + p.new*dnbinom(x,k.start,mu = exp(alpha + beta*t))))
        }
      }
      
      ## estimate lambda, alpha, beta separately
      est.param <- mle2(minuslogl = Log.Like.par, start = list(lambda = 0.001, alpha = 2, beta= -1),
                        method = 'Nelder-Mead', parnames = c('lambda', 'alpha', 'beta'))
    }else{
      ## estimate p with EM style
      p.new <- sum((dnbinom(x, size = k.vector[length(k.vector)], mu = exp(cur.param[2]+cur.param[3]*t))*p.opt[length(p.opt)])/
                     ((dnbinom(x,size = k.vector[length(k.vector)], mu = exp(cur.param[2]+cur.param[3]*t))*p.opt[length(p.opt)]) +
                        (1-p.opt[length(p.opt)])*dpois(x, cur.param[1])))/n
      
      if (is.null(x_exp1) == FALSE & is.null(x_exp0) == FALSE){
        Log.Like.par <- function (lambda,alpha,beta){
          -sum(log((1-p.new)*dpois(x,lambda) + p.new*dnbinom(x,best.k,mu = exp(alpha + beta*t)))) -
            sum(log(dnbinom(x_exp1, best.k, mu = exp(alpha + beta*t_exp1)))) - sum(log(dpois(x_exp0, lambda)))
        }
      } else if (is.null(x_exp1) == FALSE & is.null(x_exp0) == TRUE) { # experimental data on shooters
        Log.Like.par <- function (lambda,alpha,beta){
          -sum(log((1-p.new)*dpois(x,lambda) + p.new*dnbinom(x,best.k,mu = exp(alpha + beta*t)))) -
            sum(log(dnbinom(x_exp1, best.k, mu = exp(alpha + beta*t_exp1))))
        }
      } else if (is.null(x_exp1) == TRUE & is.null(x_exp0) == FALSE) { # experimental data on non shooters
        Log.Like.par <- function (lambda,alpha,beta){
          -sum(log((1-p.new)*dpois(x,lambda) + p.new*dnbinom(x,best.k,mu = exp(alpha + beta*t)))) -
            sum(log(dpois(x_exp0, lambda)))
        }
      } else { # no experimental data
        Log.Like.par <- function (lambda,alpha,beta){
          -sum(log((1-p.new)*dpois(x,lambda) + p.new*dnbinom(x,best.k,mu = exp(alpha + beta*t))))
        }
      }
      ## estimate lambda, alpha, beta separately
      est.param <- mle2(minuslogl = Log.Like.par, start = list(lambda = cur.param[1],
                                                               alpha = cur.param[2],beta= cur.param[3]),
                        method = 'Nelder-Mead', parnames = c('lambda', 'alpha', 'beta'))
    }
    cur.param <- est.param@coef ## set current parameters
    param.convergence <- c(param.convergence, est.param@details$convergence) ## set convergence parameters
    ## estimate k separatly with the updated  estimated other paramters
    if (iter == 1){
      est.k <- optim(par = 0.5, fn = Log.Like, p = p.new, lambda = cur.param[1],
                     alpha = cur.param[2], beta = cur.param[3], method = 'Brent', lower = 0.1, upper = 100)
    }else{
      est.k <- optim(par = k.vector[iter-1], fn = Log.Like, p = p.new, lambda = cur.param[1],
                     alpha = cur.param[2], beta = cur.param[3], method = 'Brent', lower = 0.1, upper = 100)
    }
    p.opt <- c(p.opt, p.new)
    best.k <- est.k$par
    param.matrix <- rbind(param.matrix,cur.param) ## add the parameters to the matrix to perform stat
    k.vector <- c(k.vector, best.k)
    if (iter==1){
      diff.par <- abs(cur.param - par.start[2:4])
      diff.k <- abs(best.k - 0.5)
    }else{
      diff.par <- abs(cur.param - param.matrix[iter-1,])
      diff.k <- abs(best.k - k.vector[iter-1])
    }
    iter <- iter + 1
  }, silent = TRUE)
  result <- c(best.k, p.new, cur.param)
  result <- matrix(result, nrow = 1)
  colnames(result) <- c("k", "p", "lambda", "alpha", "beta")
  sum_hess <- 0
  sum_mix <- 0
  sum_grad <- 0
  for (m in 1:length(x)){
    sum_grad <- sum_grad + Grad(x[m], t[m], p.new, cur.param[1], best.k, cur.param[2], cur.param[3])
    sum_hess <- sum_hess + Hess(x[m], t[m], p.new, cur.param[1], best.k, cur.param[2], cur.param[3])
    sum_mix <- sum_mix + Mix_Mat_diag(x[m], t[m], p.new, cur.param[1], best.k, cur.param[2], cur.param[3])
  }
  res1 <- 0
  res2 <- 0
  cnt <- 0
  for (j in 1:(length(x)-1)){
    for (f in (j+1):length(x)){
      res1 <- res1 + Grad(x[j], t[j],  p.new, cur.param[1], best.k, cur.param[2], cur.param[3]) %*% t(Grad(x[f], t[f], p.new, cur.param[1], best.k, cur.param[2], cur.param[3]))
      res2 <- res2 + Grad(x[f], t[f],  p.new, cur.param[1], best.k, cur.param[2], cur.param[3]) %*% t(Grad(x[j], t[j], p.new, cur.param[1], best.k, cur.param[2], cur.param[3]))
      cnt <- cnt + 1
    }
  }
  inf_mat <- sum_hess - sum_mix - res1 - res2
  matrices <- Compute_Variance(inf_mat)
  inf_obs <- -solve(hessian(Log.Like.1, c(p.new, cur.param[1], best.k, cur.param[2], cur.param[3])))
  colnames(matrices) <- c( "p", "lambda", "k", "alpha", "beta")
  rownames(matrices) <- c("p", "lambda", "k", "alpha", "beta")
  options(warn = oldw)
  if (iter == max_iter){
    return('Max iteration reach without convergence')
  }
  return(list(result = result, variance_matrix_result = matrices, inf_obs = inf_obs))
}

## Calculating one sided CI for confidence_level (alpha) wich corresponding to lower bound of LR
## on LR_values_boot_list returned by Bootstrap_simulation function
## Parameters: LR_values_boot_list: A list of LR values simulated from bootstrap
##             confidence_level: alpha value to calculate lower bound CI through quantile function
## return: list of coresponding one way CI's, for each GSR a vector corresponding to the CI per time after shooting
Calculate_CI_for_LR_with_EQ <- function(LR_values_boot_list, confidence_level){
  number_of_GSR <- nrow(LR_values_boot_list[[1]])
  matrix_values <- do.call(rbind,lapply(LR_values_boot_list,matrix,ncol=ncol(LR_values_boot_list[[1]])))
  list_of_result <- list()
  for (i in 1:number_of_GSR){
    temp_mat <- NULL
    ind <- seq(from = i, to = nrow(matrix_values), by = number_of_GSR)
    for (j in ind){
      if(j <= nrow(matrix_values)){
        temp_mat <- rbind(temp_mat, matrix_values[j,])
      }
    }
    list_of_result[[i]] <- apply(temp_mat, 2, quantile, probs = confidence_level, na.rm = T)
  }
  return(list_of_result)
}



#####################################################################################
#                                                                                   #
#                                                                                   #
#                           SIMULATIONS EXAMPLE                                     #
#                                                                                   #
#                                                                                   #
#####################################################################################


n <- 500
alpha <- 5
beta <- -0.23
lambda <- 0.5
k <- 2.5
p <- 0.5

## set up sample
t <- rep(1:10, n)[1:n]
y <- rpois(n = n,lambda)
z <- rnbinom(n = n, size = k, mu = exp(alpha + beta*t))
e <- (runif(n = n) < p)
x <- (1-e)*y + e*z 

estimated_parameters <- estimate_model(x, t)
CI_theta <- CI_for_theta(estimated_parameters$result, list(estimated_parameters$variance_matrix_result), 0.95)
gunshot_values <- 0:5
time_after_shooting_values <- c(2, 4, 6, 8, 10)
CI_LR_Mandel_method <- Calculate_CI_for_LR(gunshot_values, time_after_shooting_values, estimated_parameters$result,
                                           list(estimated_parameters$variance_matrix_result), 0.95)
boot_values <- Bootstrap_simulation(1000, estimated_parameters$result, length(x), gunshot_values,
                                    time_after_shooting_values)
EQ_CI <- Calculate_CI_for_LR_with_EQ(boot_values, 0.05)

estimated_parameters_with_experimental_data <- estimate_model(x = x, t = t, x_exp0 = rep(0, 25))
