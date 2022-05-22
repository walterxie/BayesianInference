# set working dir here
# setwd("~/WorkSpace/BayesianInference")

# This function calculates the probability of a 
# sequence of coin flip data, given a value of P_heads_guess
calc_prob_coin_flip_data <- function(P_heads_guess, coin_flips, log=T)
{
  # Empty list of probabilities
  probs_list = rep(NA, times=length(coin_flips))
  probs_list
  
  for (i in 1:length(coin_flips))
  {
    # Print an update
    #cat("\nAnalysing coin flip #", i, "/", length(coin_flips), sep="")
    
    # Get the current coin flip
    coin_flip = coin_flips[i]
    
    # If the coin flip is heads, give that datum
    # probability P_heads_guess.
    # If tails, give it (1-P_heads_guess)
    
    if (coin_flip == "H")
    {
      probs_list[i] = P_heads_guess
    } # End if heads
    
    if (coin_flip == "T")        
    {
      probs_list[i] = (1-P_heads_guess)
    } # End if tails
  } # End for-loop
  
  # Look at the resulting probabilities
  probs_list
  
  # We get the probability of all the data by multiplying
  # all the probabilities
  likelihood_of_data_given_P_heads_guess = prod(probs_list)
  
  if (log) 
    likelihood_of_data_given_P_heads_guess = log(likelihood_of_data_given_P_heads_guess)
  
  # Return result
  return(likelihood_of_data_given_P_heads_guess)
}


# this must be log likelihood
likelihood = function(param){
  return(calc_prob_coin_flip_data(P_heads_guess=param, coin_flips=coin_flips))
}

# prior 
prior = function(param, shape1 = 2, shape2 = 5){
  return(dbeta(param, shape1 = shape1, shape2 = shape2, log = T))
}

# make a proposal in MCMC
getProposal = function(current, n=1, w = .01) {
  proposal = current + runif(n, min = -(w/2), max = (w/2)) 
  return(proposal)
}

# accept or not the currect state 
accept = function(proposal, current) {
  if (proposal < 0 | proposal > 1) 
    return(FALSE);
  ratio = exp( likelihood(proposal) + 
                 prior(proposal) - 
                 likelihood(current) - 
                 prior(current) )
  if (ratio > runif(1)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# MCMC given a start value and how many iterations
run_metropolis_MCMC = function(startvalue, iterations){
  cat("Prior beta(2, 10)\n")
  chain = array(dim = c(iterations+1,2))
  chain[1,1] = startvalue
  # record likelihood not log-scale
  chain[1,2] = calc_prob_coin_flip_data(P_heads_guess=startvalue, coin_flips=coin_flips, log = F)
  reject = 0
  accept = 0
  i=1
  while (i <= iterations) {
    #cat("i = ", i, " : ")
    proposal = getProposal(chain[i,1], w=.2) 
    isAccepted = accept(proposal, chain[i,1])
    
    if (isAccepted) {
      #cat(" accept ", proposal, "\n")
      chain[i+1,1] = proposal
      chain[i+1,2] = calc_prob_coin_flip_data(P_heads_guess=proposal, coin_flips=coin_flips, log = F)
      accept = accept + 1
    } else {
      #cat(" reject ", proposal, "\n")
      chain[i+1,] = chain[i,]
      reject = reject + 1
    }
    i = i + 1
  }
  cat("\nAccept = ", accept, ", reject  = ", reject, "\n")
  
  return(mcmc(chain))
}

# plot prior, posterior, and likelihood together
createPlot = function(chain, chain_len, shape1 = 2, shape2 = 5, burnin=.1) {
  # rm burnin
  posterior = chain[(chain_len*burnin):chain_len, 1]
  length(posterior)
  
  llhd = chain[(chain_len*burnin):chain_len, 2]
  # rescaled by max of posterior
  llhd = llhd / max(llhd) * 0.9
    
  library("tidyverse")
  df1 = tibble(x = posterior, dist = rep("posterior", length(posterior))) 
  df2 = tibble(x = llhd, dist = rep("likelihood", length(llhd)))
  
  library("ggplot2")
  p <- ggplot(df1, aes(x=x)) + 
    geom_density(aes(colour = "posterior")) + 
    stat_function(fun = function(x) dbeta(x, shape1 = shape1, shape2 = shape2), 
                  aes(colour = "prior")) +
    geom_density(aes(x=x, colour = "likelihood"), data = df2) +
    xlim(0, 1) + scale_colour_manual(values = c("orange", "green", "blue")) + 
    xlab("") + ylab("") + theme_classic() 
  return(p)
}


library(coda)
startvalue = c(0.5) # P_heads_guess
chain_len = 10000

###### case 1: 10 flips ###### 
coin_flips = c('H','T','H','T','H','H','T','H','H','H')
length(coin_flips)

chain = run_metropolis_MCMC(startvalue, chain_len)
summary(chain[(chain_len*.1):chain_len,])

p <- createPlot(chain, chain_len)
p

ggsave("b25d10l.pdf", p, width = 5, height = 5)


###### case 2: 50 flips ###### 

coin_flips = c('H','T','H','T','H','H','T','H','H','H','T','H','H','T','T','T','T','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','T','T','T','H','T','T','T','H','T','T','T','H','H','H','T','T','H')
length(coin_flips)
chain = run_metropolis_MCMC(startvalue, chain_len)
summary(chain[(chain_len*.1):chain_len,])

p <- createPlot(chain, chain_len)
p

ggsave("b25d50l.pdf", p, width = 5, height = 5)


###### case 3: 100 flips ###### 

coin_flips = c('H','T','H','T','H','H','T','H','H','H','T','H','H','T','T','T','T','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','T','T','T','H','T','T','T','H','T','T','T','H','H','H','T','T','H','H','H','T','H','H','H','T','T','H','H','H','H','H','H','H','T','T','H','H','H','H','T','T','H','H','H','T','T','H','H','H','H','H','H','T','T','T','H','H','H','H','H','H','T','H','T','H','H','T','T')
length(coin_flips)
chain = run_metropolis_MCMC(startvalue, chain_len)
summary(chain[(chain_len*.1):chain_len,])

p <- createPlot(chain, chain_len)
p

ggsave("b25d100l.pdf", p, width = 5, height = 5)



