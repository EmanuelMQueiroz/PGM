## Main Code ##

include("projections.jl")
include("PG1.jl")
include("PG2.jl")
include("PG3.jl")
include("PG4.jl")

using LinearAlgebra, DataFrames, Random

# Problems #

# 1 - Hiper-Rotated Elipsoid

function elipsoid(x)
   n = length(x)
   outer = 0

   for k in 1:n 
       inner = 0
       for j in 1:k
           inner = inner + x[j]^2
       end
       outer = outer + inner
   end

   return outer 
end

function ∇elipsoid(x)
   n = length(x)
   grad = zeros(n) 

   for i in 1:n
       grad[i] = 2 * x[i] * (n - i + 1)
   end

   return grad
end

# 2 - Rosenbrock Function 

function rosenbrock(x::Vector)
   n = length(x)
   f = 0.0
   for i in 1:(n-1)
       f += 100.0 * (x[i+1] - x[i]^2)^2 + (1.0 - x[i])^2
   end
   return f
end

function ∇rosenbrock(x::Vector{Float64})
   n = length(x)
   g = zeros(Float64, n)
   for i in 1:n
       if i == 1
           g[i] = -400.0 * x[i] * (x[i+1] - x[i]^2) + 2.0 * (x[i] - 1.0)
       elseif i == n
           g[i] = 200.0 * (x[i] - x[i-1]^2)
       else
           g[i] = 200.0 * (x[i] - x[i-1]^2) - 400.0 * x[i] * (x[i+1] - x[i]^2) + 2.0 * (x[i] - 1.0)
       end
   end
   return g
end

## Parameters 
guess = MersenneTwister(1234567)
dimension = 5                             # Dimension
x_rand = 20 * rand(guess, dimension) .-10 # Define the guess in interval [-10, 10] in each variable
projection = projection6                  # Select the projection operator (Options: projection1, projection2, projection3, projection4, projection5 and projection6)
#(A,b) = matrix(3, dimension)             # Only for projection6
x0 = projection(x_rand)                   # Projects the initial guess into the set
σ = 1.e-4                                 # Linesearch parameter
ε = 1.e-5                                 # Convergence tolerance parameter
β_start = 1.0                             # Initial step length parameter 
β1 = 0.1                                  # Lower limit of the range of step length parameters calculated by quadratic interpolation
β2 = 0.9                                  # Upper limit of the range of step length parameters calculated by quadratic interpolation
γ_start = 1.0                             # Initial step length parameter for PG1, PG2 and PG3 strategies
min_step = 1.e-5                          # Tolerance for minimum step length
max_iter = 50000                          # Maximum number of iterations allowed
M = 4                                     # Non-monotonicity parameter for the PG2 strategy
η = 0.7                                   # Non-monotonicity parameter for the PG3 strategy
objective = rosenbrock                    # Objective function
strategy = "PG1"                          # Strategy employed (Options: PG1, PG2, PG3 and PG4)

if objective == elipsoid
   gradf = ∇elipsoid
   elseif objective == rosenbrock
   gradf = ∇rosenbrock
end

if strategy == "PG1"
   (x, f(x), info, et, ierror, seqx, evalsf, evalsproj) = method1(x0, objective, gradf, ε, max_iter, PG1, projection)
   elseif strategy == "PG2"
   (x, f(x), info, et, ierror, seqx, evalsf, evalsproj) = method2(x0, objective, gradf, ε, M, max_iter, PG2, projection)
   elseif strategy == "PG3"
   (x, f(x), info, et, ierror, seqx, evalsf, evalsproj) = method3(x0, objective, gradf, ε, η, max_iter, PG3, projection) 
   elseif strategy == "PG4"
   (x, f(x), info, et, ierror, seqx, evalsf, evalsproj) = method4(x0, objective, gradf, ε, max_iter, PG4, projection)
end

# Show the result #

ENV["LINES"] = 10000
println(info)
println("Minimum value of f: ", objective(x))
println("Minimizador = ", x)
println("Total time spent: ", et)
println("Function evaluations = ", evalsf)
println("Ierror = ", ierror)
