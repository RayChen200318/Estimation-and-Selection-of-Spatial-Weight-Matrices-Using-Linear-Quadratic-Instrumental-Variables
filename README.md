# Estimation-and-Selection-of-Spatial-Weight-Matrices-Using-Linear-Quadratic-Instrumental-Variables
Building upon recent developments in spatial econometric models that address the misspecification of spatial weight matrices through adaptive LASSO techniques, my research aims to enhance the selection of instrumental variables within this framework. Specifically, I introduce linear quadratic forms of instrumental variables to improve the estimation properties of the model. This approach seeks to refine the estimation process by leveraging more informative instruments, thereby achieving better convergence rates and oracle properties. By incorporating these advanced IV selection methods, my work contributes to a more robust and efficient unified framework for the estimation and selection of spatial weight matrices in spatial econometric analysis.

However, I encountered issues with my self-developed optimization algorithm, where the penalized objective function failed to shrink variables to zero. As a result, I am systematically studying optimization techniques and relevant concepts related to LASSO to address this problem.

Sparse Matrix Data Generating Process
=====
This script generates a large-scale panel dataset with sparse matrix structures for econometric analysis.


BCD_lasso stage
=====
This script implements a block coordinate descent algorithm for LASSO-based optimization. The process involves iterating over the elements of matrix `A`, estimating variables `beta`, `delta`, and performing convergence checks. 


BCD_adaptive lasso stage
=====
This script implements an Adaptive LASSO-based block coordinate descent algorithm. It iteratively optimizes the elements of matrix `A`, `delta`, and `beta` with adaptive weights.
