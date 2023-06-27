# SBM functions source file

# G is binary adjacency matrix for co-occurrence network
# Q is binary matrix for taxonomic tree information
# k is the fixed number of communities
# f is a positive constant or zero for the Markov random field prior setting

Bayesian_SBM_MRF <- function(G, Q, k, f){
  
  # number of taxa
  p <- ncol(G)
  
  # initialize random block membership z
  z <- sample(1:k, size = p, replace = TRUE)
  
  # initial settings
  T <- 1000   # iterations
  B <- T/2   # burnins
  a_omega <- 1  # beta prior shape 1
  b_omega <- 1  # beta prior shape 2
  f <- f  # MRF param
  counter <- 0
  
  # storage for parameters
  omega_store <- array(data = NA, dim = c(T,k,k))
  z_store <- array(data = NA, dim = c(T,p))
  
  # algorithm
  for(t in 1:T){
    # Monitor the Process
    if(t*100/T >= counter){
      print(paste0(counter, "% has been completed with k = ",k,"."))
      counter <- counter + 10
    }
    
    # update omegas
    for(i in 1:k){
      for(j in 1:k){
        # k and k' are the same group
        if(i == j){
          N <- choose(length(which(z == i)),2)
          M <- sum(G[z == i, z == j])/2
          w <- rbeta(1, a_omega + M, b_omega + N - M)
          omega_store[t,i,j] <- w 
        }
        # k and k' are not the same group
        if(i < j){
          N <- length(which(z == i))*length(which(z == j))
          M <- sum(G[z == i, z == j])
          w <- rbeta(1, a_omega + M, b_omega + N - M)
          omega_store[t,i,j] <- w
          omega_store[t,j,i] <- w
        }
      }
    }
    
    #update z (discrete uniform prior)
    eta <- rep(1/k, k)
    
    # MRF prior
    for (i in 1:p) {
      logprob_temp <- rep(0,k)
      for (K in 1:k) {

        logprob_temp[K] <- log(eta[K]) + f*sum(z[which(Q[i,] == 1)] == K)
        
        # Likelihood part
        for (ii in 1:p) {
          if (ii == i) {
            next
          } else {
            if (K < z[ii]) {
              logprob_temp[K] <- logprob_temp[K] + log(omega_store[t, K, z[ii]])*G[i, ii] + log(1 - omega_store[t, K, z[ii]])*(1 - G[i, ii])
            } else {
              logprob_temp[K] <- logprob_temp[K] + log(omega_store[t, z[ii], K])*G[i, ii] + log(1 - omega_store[t, z[ii], K])*(1 - G[i, ii])
            }
          }
        }
      }
      
      # normalize density
      logprob_temp <- logprob_temp - max(logprob_temp)
      prob_temp <- exp(logprob_temp)
      probs <- prob_temp/sum(prob_temp)
      # sample community label
      z[i] <- sample(1:k, 1, prob = probs)
    }
    
    # store z
    z_store[t,] <- z
    
  }
  
  results <- list(z_store, omega_store)
  return(results)
}
