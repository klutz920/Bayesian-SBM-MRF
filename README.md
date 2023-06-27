# Bayesian-SBM-MRF
 
## Tutorial

1. Save/load the function code in R as a source function or download it from GitHub.
2. Process the abundance data as described in our paper to estimate the co-occurrence network G. Our network in the real data analysis consisted of species. The MCLR transformation is available in the SPRING package in R.
3. Matrix Q should be from one other taxonomic tree level such as genus.
4. You can run the code for one value of K and f at a time. The function returns all 1000 iterations for parameters Omega and z. You should discard the first half as burn in. Then, use our formulation from the paper to compute the posterior density, MAP, BIC, etc. To compute ARI for simulation data, use the ARI function in the aricode package in R.
5. Sample codeonce you hav
