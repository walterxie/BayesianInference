# S I R
# http://alexeidrummond.org/bayesian_phylo_lectures/lecture13/SIR_demo.js

Tinf <- function(S,I,R) {
  return(pub.beta * S * I)
}

Trem <- function(S,I,R) {
  return(pub.mu*I)
}

# T is final time
getStochastic <- function(T) {
  S = c(pub.N-1) 
  I=c(1)
  R=c(0)
  t=(0)
  
  i = 1 # R index starts from 1
  while (t[i] <= T) {
    Ti = Tinf(S[i], I[i], R[i])
    Tr = Trem(S[i], I[i], R[i])
    T0 = Ti + Tr
    
    if (T0 == 0)
      break
    
    dt = -log(runif(1))/T0
    if (t[i] + dt > T)
      break
    
    t = c(t, t[i] + dt)
    
    if (runif(1)*T0 < Ti) {
      # Infection
      S = c(S, S[i]-1)
      I = c(I, I[i]+1)
      R = c(R, R[i])
    } else {
      # Recovery
      S = c(S, S[i])
      I = c(I, I[i]-1)
      R = c(R, R[i]+1)
    }
    
    i = i + 1
  }
  
  S = c(S, S[i])
  I = c(I, I[i])
  R = c(R, R[i])
  t = c(t, T)
  R0 = getR0(S,I,R,t)
  
  return(list(S=S,I=I,R=R,t=t, R0=R0))
}

getDeterministic <- function(T, nSteps) {
  
  dt = T/(nSteps - 1)
  
  S = c(pub.N-1) 
  I=c(1)
  R=c(0)
  t=(0)
  
  # R index starts from 1
  for (i in 2:nSteps) {
    #Semi-implicit algorithm (uses numerical estimate of dx/dt at t+0.5dt)
    Sp = c(S[i-1])
    Ip = c(I[i-1])
    Rp = c(R[i-1]) 
    
    for (iter in 1:3) {
      Spp = S[i-1] - Tinf(Sp,Ip,Rp)*0.5*dt
      Ipp = I[i-1] + (Tinf(Sp,Ip,Rp) - Trem(Sp,Ip,Rp))*0.5*dt
      Rpp = R[i-1] + Trem(Sp,Ip,Rp)*0.5*dt
      
      Sp = Spp
      Ip = Ipp
      Rp = Rpp
    }
    
    S = c(S, Sp*2 - S[i-1])
    I = c(I, Ip*2 - I[i-1])
    R = c(R, Rp*2 - R[i-1])
    t = c(t, t[i-1] + dt)
  }
  
  R0 = getR0(S,I,R,t)
  return(list(S=S,I=I,R=R,t=t,R0=R0))
}

getR0 <- function(S,I,R,t) {
  R0 = c()
  for (j in 1:length(t)) {
    R0 = c(R0, Tinf(S[j],I[j],R[j])/Trem(S[j],I[j],R[j]) )
  }
  return(R0)
}

createPlot <- function(result=list()) {
  library(tidyverse)
  print(names(sto))
  df <- data.frame(sto) %>% as_tibble %>% gather("name", "y", -t)
  
  
  library(ggplot2)
  p <- ggplot(df, aes(x=t, y=y, group=name, color=name)) +
    geom_line() + xlab("Time") + ylab("Population size") + theme_classic() 
  return(p)
}



pub.N = 100
pub.beta = 0.01
pub.mu = 0.1

sto <- getStochastic(50)
sto

p <- createPlot(sto)
p
ggsave("SIRstochastic.pdf", p, width = 5, height = 5)

det <- getDeterministic(50, nSteps=50)
det

p <- createPlot(det)
p
ggsave("SIRdeterministic.pdf", p, width = 5, height = 5)
