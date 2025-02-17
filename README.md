# SARIS algorithm

This the code used to reproduce the simulations of the preprint "Guédon, T., Baey, C., & Kuhn, E. (2024). Estimation of ratios of normalizing constants using stochastic approximation: the SARIS algorithm. arXiv preprint arXiv:2408.13022."

## Description

Two scripts :
- utils.R : contains the auxiliary functions used for the simulations, the main ones are : 
  - metropolis(alpha, sample_size, **args) : metropolis-hasting sampler. return a sample of size sample_size, from the distribution with density proportional to the function alpha
  - bissection(f, a, b, **args) : return the approximated root of the funciton f
  - BridgeSARIS(sampletot, f0, f1, r0, **args) : compute the bridge like SARIS estimator based on a sample sampletot$\sim 0.5(p_0+p_1)$ $$r_{k+1} = r_k + \frac{f_0(Z{k+1}) - r_k f_1(Z{k+1})}{f_0(Z{k+1}) + r_k f_1(Z{k+1})}, Z_{k+1}\sim 0.5(p_0+p_1)$$
  - SARIS_RM(f0, f1, alpha, **args) : compute the SARIS-EXT estimator using the density proportional to alpha(.,r) : $$r_{k+1} = r_k + \frac{f_0(Z{k+1}) - r_k f_1(Z{k+1})}{alpha(r_k, Z_{k+1})}, Z_{k+1}\sim alpha(r_k,.)$$
- run_simu_SARIS : scripts to execute in order to obtain the different figures of the article:
  - Required Packages:
    - dplyr (1.1.4)
    - ggplot2 (3.5.0)
    - reshape2 (1.4.4)
    - gridExtra (2.3)
  - Usage:
    - This script shall be executed linearly
