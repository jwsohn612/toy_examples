# Data Scenario
# y ~ N(mu, sigma^2) and mu ~ N(0,1)
# sigma is fixed 
# 
# Reference : Handbook of Markov Chain Monte Carlo, Brooks,et al.
# 

# Data Configuration
N <- 1000
mu <- 3; sigma <- 2
y <- rnorm(N, mu, sigma)

# Define potential energy function : target function
U <- function(y, mu, sigma){
  prior_density <- sum(dnorm(mu, mean=0 ,sd=1, log=TRUE))
  likelihood_density <- sum(dnorm(y, mean=mu,sd=sigma,log=TRUE))
  return(-(likelihood_density+prior_density))
}

# Define its derivative with respect to the target variable
grad_mu <- function(y, mu, sigma) mu - sum(y-mu)/(sigma^2)

# This will be implemented for each iteration.
HMC <- function(y, U, grad_mu, epsilon, L, current_q){
  
  q <- current_q 
  p <- rnorm(length(q), 0, 1) # m is set 1
  current_p <- p 
  
  # Leapfrog method
  for(i in 1:L){
    p <- p - epsilon/2 * grad_mu(y, q, sigma)
    q <- q + epsilon * p 
    p <- p - epsilon/2 * grad_mu(y, q, sigma)
  }
  p <- -p
  
  current_U <- U(y, current_q, sigma)
  current_K <- sum(current_p^2)/2
  
  proposed_U <- U(y, q, sigma)
  proposed_K <- sum(p^2)/2
  
  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)){
    return(q)
  }else{
    return(current_q)
  }
}

# The total iteration number is T
T <- 1000
posterior_mu <- c()
for(t in 1:T){
  q <- HMC(y, U, grad_mu, epsilon = 0.01, L=50, current_q = quantile(y, 0.5))
  posterior_mu <- append(posterior_mu, q)
}

# Describe results
ts.plot(posterior_mu)
hist(posterior_mu)



# HMC <- function(U, grad_U, epsilon, L, current_q){
#   
#   q <- current_q 
#   p <- rnorm(length(q), 0, 1) # m is set 1
#   current_p <- p 
#   
#   # Leapfrog method
#   for(i in 1:L){
#     p <- p - epsilon/2 * grad_U(q)
#     q <- q + epsilon * p 
#     p <- p - epsilon/2 * grad_U(q)
#   }
#   
#   p <- -p
#   
#   current_U <- U(current_q)
#   current_K <- sum(current_p^2)/2
#   
#   proposed_U <- U(q)
#   proposed_K <- sum(p^2)/2
#   
#   if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)){
#     return(q)
#   }else{
#     return(current_q)
#   }
# }
