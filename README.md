# Estimation-and-Selection-of-Spatial-Weight-Matrices-Using-Linear-Quadratic-Instrumental-Variables
Building upon recent developments in spatial econometric models that address the misspecification of spatial weight matrices through adaptive LASSO techniques, my research aims to enhance the selection of instrumental variables within this framework. Specifically, I introduce linear quadratic forms of instrumental variables to improve the estimation properties of the model. This approach seeks to refine the estimation process by leveraging more informative instruments, thereby achieving better convergence rates and oracle properties. By incorporating these advanced IV selection methods, my work contributes to a more robust and efficient unified framework for the estimation and selection of spatial weight matrices in spatial econometric analysis.

However, I encountered issues with my self-developed optimization algorithm, where the penalized objective function failed to shrink variables to zero. As a result, I am systematically studying optimization techniques and relevant concepts related to LASSO to address this problem.

Sparse Matrix Data Generating Process
=====
# Sparse Matrix Data Generating Process

This script generates a large-scale panel dataset with sparse matrix structures for econometric analysis. Below is a summary of the key steps involved:

## 1. Initial Parameters and Matrix Generation
- **Adjacency Matrix (`A_star`)**: A sparse `N×N` matrix is created, with 5% of off-diagonal elements set to 0.5, normalized by row.
- **Sparse Matrix Sequence (`prespecified_spm_sequence`)**: A list of `M` binary `N×N` sparse matrices is generated, representing binary relationships between nodes over time.
  
## 2. Error and Covariate Generation
- **Error Sequence (`error_sequence`)**: Disturbance terms are generated using a multivariate normal distribution with a random covariance matrix structure.
- **Covariate Sequence (`covariate_sequence`)**: Covariates are generated with a time-varying structure, adjusted by the disturbance terms.

## 3. Instrumental Variables and Demeaning
- **Instrumental Variables (`IV_sequence`)**: Covariates and the sparse matrices are combined to generate time-varying instrumental variables.
- **Demeaning Process**: The average of IVs is subtracted to obtain demeaned IVs, which are used for further transformations.

## 4. Transformation and Full Matrix Estimation
- **Dependent and Independent Variable Transformations**: The demeaned IVs are used to transform the dependent variables (`dependence_tilde_sequence`) and covariates (`covariate_tilde_sequence`).
- **Estimation Matrices (`B_full`, `K_full`, `H_full`, `V_full`)**: These matrices are constructed for econometric estimation, using sparse matrix operations to handle large data efficiently.

## 5. Sparse Matrix Optimization
- The script heavily relies on sparse matrix operations via the `Matrix` package to reduce memory usage and improve computational performance, making it suitable for large-scale panel data analysis.

BCD_lasso stage
=====
# Block Coordinate Descent Algorithm (LASSO Stage)

This script implements a block coordinate descent algorithm for LASSO-based optimization. The process involves iterating over the elements of matrix `A`, estimating variables `beta`, `delta`, and performing convergence checks. Below is a brief overview of key functions and their roles:

## 1. **`constraint_functions`**
- **Purpose**: Defines the inequality constraints for the optimization problem.
- **Inputs**: 
  - `x`: Vector of variables to optimize.
- **Output**: A matrix representing the constraints on the matrix `A` and the weights `delta`.
- **Role**: Ensures that the rows of matrix `A` are properly normalized and that `delta` satisfies the constraint.

## 2. **`objective_functions`**
- **Purpose**: Defines the objective function to minimize, including the LASSO penalty.
- **Inputs**: 
  - `x`: Vector of variables to optimize (elements of `A`).
- **Output**: The value of the objective function (a scalar).
- **Role**: Calculates the error between the predicted and actual values, while applying the LASSO regularization on the matrix `A`.

## 3. **`objective_functions_grad`**
- **Purpose**: Computes the gradient of the objective function.
- **Inputs**: 
  - `x`: Vector of variables to optimize (elements of `A`).
- **Output**: The gradient of the objective function with respect to `x`.
- **Role**: Helps the optimization algorithm to follow the gradient for faster convergence.

## 4. **`constraint_functions_grad`**
- **Purpose**: Computes the Jacobian matrix (gradient of the constraints).
- **Inputs**: 
  - `x`: Vector of variables to optimize.
- **Output**: Jacobian matrix of the constraint functions.
- **Role**: Provides the necessary gradients for constraint optimization.

## 5. **Optimization Loop**
- **Purpose**: Sequentially optimizes each row of matrix `A` using the `nloptr` package.
- **Inputs**: 
  - Initial values (`x0`), objective function, constraints, and their gradients.
- **Output**: Optimized solutions for the elements of `A`.
- **Role**: Implements block coordinate descent by optimizing each row of matrix `A` independently.

## 6. **Convergence Check**
- **Purpose**: Ensures that the algorithm stops when the results converge.
- **Role**: Compares the difference between subsequent solutions for `eta` and terminates the loop if the change is below a set threshold (`pred_rel`).

## 7. **Result Storage**
- **Purpose**: Stores the results of each iteration, including the optimized values of `A`, `beta`, and `delta`.
- **Output**: Final estimated matrices `A_tilde`, `beta_tilde`, and `delta_tilde`.

## 8. **Progress Bar (`progress_bar`)**
- **Purpose**: Displays the progress of the iterative process.
- **Role**: Updates after each iteration to show completion percentage.

## Final Output:
- Optimized matrices `A_tilde`, `beta_tilde`, `delta_tilde`, and the final LASSO results.


BCD_adaptive lasso stage
=====
# Adaptive LASSO Algorithm

This script implements an Adaptive LASSO-based block coordinate descent algorithm. It iteratively optimizes the elements of matrix `A`, `delta`, and `beta` with adaptive weights.
