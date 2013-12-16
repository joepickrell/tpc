Test for population continuity

R code for the TPC from Brandt et al. (2013) 

Requirements:

R library MCMCpack

R library mvtnorm

Usage:

R --vanilla --args [index1] [index2] [prior on c] [fpropose] [tpropose] [input file] [output file] < MCMC_new.R

where:

[index1] is the line number in the input file containing the data from the earlier population in time (excluding the header)

[index2] is the line number in the input file containing the data from the later population in time (exclusing the header)

[prior on c] is the prior mean on the drift parameter t/N, where t is the number of generations separating the populations and N is the effective population size

[fpropose] is the tuning parameter for the update of the vector of haplotype frequencies f

[tpropose] is the variance on the update of the vector of thetas

[input file] the input file, format described below

[output file] the output file 



Input file:

Example in sim166.in. Each line contains the population id, the total number of haplotypes, and the counts of each individual haplotype. The first line in the header.



Output file:

Stores the output of all parameters every 10 iterations. At the end of the file is the posterior predictive P-value.



Example:

cd R/

R --vanilla --args  1 2 0.04 3000 0.00001 ../sim166.in sim166.out < MCMC_new.R

