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

# Parameters 
guess = MersenneTwister(1234567)
x_rand = 20 * rand(guess, 10) .-10
n = length(x_rand) 
projection = projection3
x0 = projection(x_rand)
σ = 1.e-4 
ε = 1.e-5
β_start = 1.0
β1 = 0.1
β2 = 0.9
γ_start = 1.0
min_step = 1.e-5
max_iter = 50000
M = 4
η = 0.9
objective = rosenbrock
strategy = "PG1"
#strategy = "PG2"
#strategy = "PG3"
#strategy = "PG4"

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
