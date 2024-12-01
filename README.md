# Uma análise numérica dos métodos do Gradiente Projetado ao longo de arcos de projeção e de direções viáveis
This repository contains the Julia codes referring to the course conclusion work entitled "Uma análise numérica dos métodos do Gradiente Projetado ao longo de arcos de projeção e de direções viáveis", authored by Emanuel Mendes Queiroz.

## main.jl
This file performs a test on the Rosenbrock or Elipsoid function using the projected gradient method equipped with one of the PG1, PG2, PG3 or PG4 strategies and generates a table with the data obtained during the optimization process.

### Code overview
The script is organized into the following sections:

1. **Imports and Inclusions:** The script imports the necessary packages (`LinearALgebra`, `DataFrames` and `Random`) and includes the essential Julia files: `projections.jl`, `GPA1.jl`, `GPA2.jl`, `GPA3.jl` e `GPA4.jl`.

2. **Rotated-Elipsoid function:**
   - **$f$ (Function):** The Rotated-Elipsoid function
   - **$\nabla f$ (Function):** The gradient of the Rotated-Elipsoid function

2. **Rosenbrock function:**
   - **$f$ (Function):** The Rosenbrock function
   - **$\nabla f$ (Function):** The gradient of the Rosenbrock function

4. **Parameters:**
Definition of all parameters to be used in the optimization process.

6. **Optimization Cycle:**
The script performs the gradient method optimization for the chosen case and show the informations at the end of process.

## GPA1.jl
This file contains the implementation of the Projected Gradient Method equipped with Armijo's linesearch along feasible directions.

### Function GPA1 - Monotone Armijo linesearch along feasible directions

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Point where the linear search will applied
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\sigma$ (FLoat64):** The parameter of Armijo linesearch
- **min_step (FLoat64):** The minimum stepsize parameter acceptable
- **$\gamma$ _start (FLoat64):** Initial stepsize
- **$zk$ (Vector):** $P_C(x_k - \beta_k\nabla f(x_k))$

#### Outputs

- **$\boldsymbol{\gamma}$ (Float64):** The stepsize parameter in iterarion k
- **ierror (Int64):** The error of linesearch (0 - no error, 1 - step length too small)
- **evalf (Int64):** Number of function evaluations in the current iteration

### Function method1 - Gradient Projected Method with monotone Armijo linesearch along feasible directions

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Initial point
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\varepsilon$ (FLoat64):** The convergence tolerance parameter
- **max_iter (Int64):** Maximum number of iterations allowed 
- **PG1 (Function):** The linesearch function
- **projection (Function):** The projection function

#### Outputs

- **$\boldsymbol{x}$ (Vector):** Minimizer point
- **$\boldsymbol{f(x)}$ (Float64):** Objective function evaluated in the minimizer
- **info (DataFrame):** Data frame with information about the minimization process
- **et (Float64):** Time spent solving the problem
- **ierror (Int64):** Error returned (0 - no error, 1 - step length too small, 2 - Maximum of iterations was achieved)
- **seqx (Matrix{Vector{Float64}}):** Matrix with the terms of the sequence $\{x_k\}$
- **sum(evalf_\gamma) (Int64):** Total number of function evaluations

## GPA2.jl
This file contains the implementation of the Projected Gradient Method equipped with nonmonotone Max type linesearch along feasible directions.

### Function GPA2 - Nonmonotone Max type linesearch along feasible directions.

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Point where the linear search will applied
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\sigma$ (FLoat64):** The parameter of Armijo linesearch
- **$M (Int64):** Non-monotonicity parameter
- **min_step (FLoat64):** The minimum stepsize parameter acceptable
- **$\gamma$ _start (FLoat64):** Initial stepsize
- **$zk$ (Vector):** $P_C(x_k - \beta_k\nabla f(x_k))$
- **fvals:** Vector with the values ​​of the objective function evaluated in previous iterations

#### Outputs

- **$\boldsymbol{\gamma}$ (Float64):** The stepsize parameter in iterarion k
- **ierror (Int64):** The error of linesearch (0 - no error, 1 - step length too small)
- **evalf (Int64):** Number of function evaluations in the current iteration
- **fmax (Float64):** Maximum of the last $M$ iterations

### Function method2 - Gradient Projected Method with monotone Armijo linesearch along feasible directions

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Initial point
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\varepsilon$ (FLoat64):** The convergence tolerance parameter
- **$M (Int64):** Non-monotonicity parameter
- **max_iter (Int64):** Maximum number of iterations allowed 
- **PG2 (Function):** The linesearch function
- **projection (Function):** The projection function

#### Outputs

- **$\boldsymbol{x}$ (Vector):** Minimizer point
- **$\boldsymbol{f(x)}$ (Float64):** Objective function evaluated in the minimizer
- **info (DataFrame):** Data frame with information about the minimization process
- **et (Float64):** Time spent solving the problem
- **ierror (Int64):** Error returned (0 - no error, 1 - step length too small, 2 - Maximum of iterations was achieved)
- **seqx (Matrix{Vector{Float64}}):** Matrix with the terms of the sequence $\{x_k\}$
- **sum(evalf_\gamma) (Int64):** Total number of function evaluations

## GPA3.jl
This file contains the implementation of the Projected Gradient Method equipped with nonmonotone Average type linesearch along feasible directions.

### Function GPA3 - Nonmonotone Average type linesearch along feasible directions

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Point where the linear search will applied
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\sigma$ (FLoat64):** The parameter of Armijo linesearch
- **min_step (FLoat64):** The minimum stepsize parameter acceptable
- **$\gamma$ _start (FLoat64):** Initial stepsize
- **$zk$ (Vector):** $P_C(x_k - \beta_k\nabla f(x_k))$
- **f_average:** Weighted average of previous iterations

#### Outputs

- **$\boldsymbol{\gamma}$ (Float64):** The stepsize parameter in iterarion k
- **ierror (Int64):** The error of linesearch (0 - no error, 1 - step length too small)
- **evalf (Int64):** Number of function evaluations in the current iteration

### Function method3 - Gradient Projected Method with nonmonotone Average type linesearch along feasible directions

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Initial point
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\varepsilon$ (FLoat64):** The convergence tolerance parameter
- **$\boldsymbol{\eta}$ (Float64):** Non-monotonicity parameter
- **max_iter (Int64):** Maximum number of iterations allowed 
- **GPA3 (Function):** The linesearch function
- **projection (Function):** The projection function

#### Outputs

- **$\boldsymbol{x}$ (Vector):** Minimizer point
- **$\boldsymbol{f(x)}$ (Float64):** Objective function evaluated in the minimizer
- **info (DataFrame):** Data frame with information about the minimization process
- **et (Float64):** Time spent solving the problem
- **ierror (Int64):** Error returned (0 - no error, 1 - step length too small, 2 - Maximum of iterations was achieved)
- **seqx )Matrix{Vector{Float64}}):** Matrix with the terms of the sequence $\{x_k\}$
- **sum(evalf_\gamma) (Int64):** Total number of function evaluations

## GPA4.jl
This file contains the implementation of the Projected Gradient Method equipped with Armijo linesearch along the projection arc.

### Function GPA2 - Armijo linesearch along the projection arc

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Point where the linear search will be carried out
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\sigma$ (FLoat64):** The parameter of Armijo linesearch
- **min_step (FLoat64):** The minimum stepsize parameter acceptable
- **$\beta$ _start (FLoat64):** Initial stepsize

#### Outputs

- **$\boldsymbol{\beta}$ (Float64):** The stepsize parameter in iterarion k
- **zkj (Vector):** $P_C(x_k - \beta_k\nabla f(x_k))$
- **ierror (Int64):** The error of linesearch (0 - no error, 1 - step length too small)
- **evalf (Int64):** Number of function evaluations

### Function method4 - Gradient Projected Method with Armijo linesearch along the projection arc

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Initial point
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\varepsilon$ (FLoat64):** The convergence tolerance parameter
- **max_iter (Int64):** Maximum number of iterations allowed 
- **GPA2 (Function):** The linesearch function
- **projection (Function):** The projection function

#### Outputs

- **$\boldsymbol{x}$ (Vector):** Minimizer point
- **$\boldsymbol{f(x)}$ (Float64):** Objective function evaluated in the minimizer
- **info (DataFrame):** Data frame with information about the minimization process
- **et (Float64):** Time spent solving the problem
- **ierror (Int64):** Error returned (0 - no error, 1 - step length too small, 2 - Maximum of iterations was achieved)
- **seqx )Matrix{Vector{Float64}}):** Matrix with the terms of the sequence $\{x_k\}$
- **sum(evalf_\gamma) (Int64):** Total number of function evaluations

## projections.jl
This file contains projection functions for 6 different feasible sets, whose parameters can be easily modified.

### Examples of projection functions for diferent feasible sets
- **projection1:** $\Omega_1 =$$`\left\{x \in \mathbb{R}^n; |x-c| \le 10 \right\}`$, onde $c = (1,1,\ldots,1)$.
- **projection2:** $\Omega_2 =$$`\left\{x \in \mathbb{R}^n; x_j \in (-\infty,5], \ j=1,2,\ldots,n \right\}`$.
- **projection3:** $\Omega_3 =$$`\left\{x \in \mathbb{R}^n; x_j \in [-5,10], \ j=1,2,\ldots,n \right\}`$.
- **projection4:** $\Omega_4 =$$`\left\{x \in \mathbb{R}^n ; \langle a, x \rangle = 10\right\}`$, onde $a = (1,1,\ldots,1)$.
- **projection5:** $\Omega_5 =$$`\left\{x \in \mathbb{R}^n ; \langle a, x \rangle \le 10\right\}`$, onde $a = (10,10,\cdots,10)$.
- **projection6:** $\Omega_6 =$$`\left\{x \in \mathbb{R}^n ; Ax = (1,1,1)\right\}`$, onde $A \in \mathbb{R}(3,n)$ tem posto $3$.
