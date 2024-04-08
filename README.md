# Bayesian-SBM-MRF
## Data File

In the data file, you will find the following:

1.  real_urinary_microbiome.csv: The real urinary microbiome dataset from our manuscript "A generalized Bayesian stochastic block model for microbiome community detection" with $p=99$ taxa and $n=75$ samples (postmenopausal women).
2.  simulation.rdata: the simulated networks with $K=3,6,9$ communities.  This was our original simulation data for the manuscript.
3.  simulation2.rdata: contains additional simulated networks with $K=12$ communities.  These were added later when we revised our paper for submission to Statistics In Medicine, which included an extension of the original simulation study.  

## Tutorial

1. Save/load the function code in R or download it from GitHub.
2. Process the abundance data as described in our paper to estimate the co-occurrence network G. Our network in the real data analysis consisted of species level taxa. The MCLR transformation is available in the SPRING package in R.
3. Matrix Q should be from one other taxonomic tree level such as genus.
4. Run the code for one value of K and f at a time. The function returns all 1000 iterations for parameters Omega and z. You should discard the first half as burn in. Then, use our formulation from the paper to compute the posterior density, MAP, BIC, etc. To compute ARI for simulation data, use the ARI function in the aricode package in R.
5. Sample code for a specified G and Q with 5 communities and our recommended MRF prior setting.

## Run function

set.seed(1)

my_results = Bayesian_SBM_MRF(G = G, 
                              Q = Q, 
                              k = 5, 
                              f = 1)

## results from all 1000 iterations

z = my_results[[1]]

omegas = my_results[[2]]

## discard first half as burn in

z = z[-c(1:500),]

omegas = omegas[-c(1:500),,]

