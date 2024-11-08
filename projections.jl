#
# Projection functions for diferent feasible sets 
#

function projection1(x)
    n = length(x)
    c = fill(10, n) 
    δ = 30 
    cx = norm(c-x)
    if cx <= δ
        return x
    else
        return c - (δ*(c-x)/cx)
    end
end

## Projection 02: C is R^n+ ##

function projection2(x)
    return min.(0, x)
end

## Projection 03: C is a hypercube of dimension n. For each dimension, the range is [-1,1] ##

function projection3(x)
    n = length(x)
    a = fill(2, n)
    b = fill(10, n)
    for j in 1:n   
        if x[j] < a[j]
           x[j] = a[j]
        elseif a[j] <= x[j] <= b[j]
           x[j] = x[j]
        elseif x[j] > b[j]
           x[j] = b[j] 
        end
    end
    return x
end

## Projection 04: C is the set of points that stisfy an equality with inner product ##

function projection4(x)
    n = length(x)
    a = fill(10, n)
    b = 10
    return x + ((b - dot(a, x))/norm(a)^2)*a
end

## Projection 05: C is the set of points that stisfy an inequality with inner product ##

function projection5(x)
    n = length(x)
    a = fill(10, n)
    b = 1
    return x + ((min(0, b - dot(a, x)))/norm(a)^2)*a
end

## Projection 06: C is the set of points that satisfy an equality with a product between a vector and a matrix with defined rank ##

function matrix(rank, n)
    rng = MersenneTwister(3214)
    matrix_0 = rand(rng, rank, n)  
    U, Σ, V = svd(matrix_0)
    Σ[rank + 1:end] .= 0      
    A = U[:, 1:rank] * Diagonal(Σ[1:rank]) * V[:, 1:rank]'
    return A
end

function projection6(x)
    n = length(x)
    rank = 3  
    B = matrix(rank, n)
    b = fill(1, size(B, 1))
    return x - B' * inv(B * B') * (B * x - b)
end
