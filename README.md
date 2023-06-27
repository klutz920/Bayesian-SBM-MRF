# Bayesian-SBM-MRF
 
## Tutorial

1. Save/load the function code in R or download it from GitHub.
2. Process the abundance data as described in our paper to estimate the co-occurrence network G. Our network in the real data analysis consisted of species level taxa. The MCLR transformation is available in the SPRING package in R.
3. Matrix Q should be from one other taxonomic tree level such as genus.
4. Run the code for one value of K and f at a time. The function returns all 1000 iterations for parameters Omega and z. You should discard the first half as burn in. Then, use our formulation from the paper to compute the posterior density, MAP, BIC, etc. To compute ARI for simulation data, use the ARI function in the aricode package in R.
5. Sample code for a specified G and Q with 5 communities and our recommended MRF prior setting.

\# Run function

set.seed(1)

my_results = Bayesian_SBM_MRF(G = G, 
                              Q = Q, 
                              k = 5, 
                              f = 1)

\# results

z = my_results[[1]]

omegas = my_results[[2]]

\# discard first half as burn in

z = z[-c(1:500),]

omegas = omegas[-c(1:500),,]

\# then compute posterior density to identify MAP, compute BIC
