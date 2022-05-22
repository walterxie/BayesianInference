
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
prior = function(param, shape1 = 1, shape2 = 1){
  return(dbeta(param, shape1 = shape1, shape2 = shape2, log = T))
}

# make a proposal in MCMC
getProposal = function(current, n=1, w = .01) {
  proposal = current + runif(n, min = -(w/2), max = (w/2)) 
  return(proposal)
}

# accept or not the currect state 
accept = function(proposal, current, prior_shape) {
  if (proposal < 0 | proposal > 1) 
    return(FALSE);
  ratio = exp( likelihood(proposal) + 
                 prior(proposal, shape1 = prior_shape, shape2 = prior_shape) - 
                 likelihood(current) - 
                 prior(current, shape1 = prior_shape, shape2 = prior_shape) )
  if (ratio > runif(1)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# MCMC given a start value and how many iterations
run_metropolis_MCMC = function(startvalue, iterations, prior_shape){
  cat("Prior beta(", prior_shape, ", ", prior_shape, ")\n")
  chain = array(dim = c(iterations+1,1))
  chain[1,] = startvalue
  reject = 0
  accept = 0
  i=1
  while (i <= iterations) {
    #cat("i = ", i, " : ")
    proposal = getProposal(chain[i,], w=.2) 
    isAccepted = accept(proposal, chain[i,], prior_shape)
    
    if (isAccepted) {
      #cat(" accept ", proposal, "\n")
      chain[i+1,] = proposal
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

# plot prior, posterior together
createPlot = function(chain, chain_len, prior_shape, burnin=.1) {
  # rm burnin
  posterior = chain[(chain_len*burnin):chain_len, 1]
  length(posterior)
  
  library("tidyverse")
  df = tibble(x = posterior, dist = rep("posterior", length(posterior)))
  
  library("ggplot2")
  p <- ggplot(df, aes(x=x)) + 
    geom_density(aes(colour = "posterior")) + 
    stat_function(fun = function(x) dbeta(x, prior_shape, prior_shape), 
                  aes(colour = "prior")) +
    xlim(0, 1) + scale_colour_manual(values = c("red", "green")) + 
    xlab("") + ylab("") + theme_classic() + 
    # Add mean line
    geom_vline(aes(xintercept=mean(posterior)), color="blue", linetype="dashed", size=.5) +
    annotate(geom="text", x=mean(posterior), y=0.2, label=round(mean(posterior),3), 
             color="blue", alpha=.7)
  return(p)
}



####### data 1 #######

# Here are 100 coin flips:
coin_flips = c('H','T','H','T','H','H','T','H','H','H','T','H','H','T','T','T','T','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','T','T','T','H','T','T','T','H','T','T','T','H','H','H','T','T','H','H','H','T','H','H','H','T','T','H','H','H','H','H','H','H','T','T','H','H','H','H','T','T','H','H','H','T','T','H','H','H','H','H','H','T','T','T','H','H','H','H','H','H','T','H','T','H','H','T','T')

# Look at the data
coin_flips
length(coin_flips[coin_flips=='H'])

# some likelihoods :
calc_prob_coin_flip_data(P_heads_guess=0.5, coin_flips=coin_flips)
calc_prob_coin_flip_data(P_heads_guess=0.6, coin_flips=coin_flips)
calc_prob_coin_flip_data(P_heads_guess=0.7, coin_flips=coin_flips)


####### case 1 #######

library(coda)

startvalue = c(0.5) # P_heads_guess
chain_len = 10000

# data 100 flip, 65 head
# prior beta(1,1)
chain = run_metropolis_MCMC(startvalue, chain_len, prior_shape=1)

summ = summary(chain[(chain_len*.1):chain_len,])
plot(chain)

p <- createPlot(chain, chain_len, prior_shape = 1)
p
ggsave("b11d100.pdf", p, width = 5, height = 5)

####### case 2 #######

# data 100 flip, 65 head
# prior beta(5,5)
chain = run_metropolis_MCMC(startvalue, chain_len, prior_shape=5)

summ = summary(chain[(chain_len*.1):chain_len,])
plot(chain)

p <- createPlot(chain, chain_len, prior_shape = 5)
p
ggsave("b55d100.pdf", p, width = 5, height = 5)

####### case 3 #######

coin_flips = c('H','H','H','H')
               
# data 4 flip, 4 head
# prior beta(1,1)
chain = run_metropolis_MCMC(startvalue, chain_len, prior_shape=1)

summ = summary(chain[(chain_len*.1):chain_len,])
plot(chain)

p <- createPlot(chain, chain_len, prior_shape = 1)
p
ggsave("b11d4.pdf", p, width = 5, height = 5)


####### case 3 #######

# data 4 flip, 4 head
# prior beta(5,5)
chain = run_metropolis_MCMC(startvalue, chain_len, prior_shape=5)

summ = summary(chain[(chain_len*.1):chain_len,])
plot(chain)

p <- createPlot(chain, chain_len, prior_shape = 5)
p
ggsave("b55d4.pdf", p, width = 5, height = 5)


