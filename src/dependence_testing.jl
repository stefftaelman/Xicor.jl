#using Xicor
include("../src/Xicor.jl")

using Random, Distributions, StatsBase


"""
    permutation_test(xi, Y; n=1000)

Computes the p-value and standard deviation of the correlation coefficient
`xi` by calculating to coefficients for `n` random vectors against `Y`.
"""
function permutation_test(ξᵢ::Float64, Y::AbstractVector; n::Int=1000)
    ξ_perm = zeros(n)
    for i in 1:n
        X_perm = rand(length(Y))
        ξ_perm[i] = ξ(X_perm, Y)
    end
    pval = mean(ξ_perm .> ξᵢ)
    sd = sqrt(var(ξ_perm))
    return pval, sd
end


#=
TO DO
""" 
    asymptotic_test(xi, Y; n=1000)

Computes the p-value and standard deviation of the correlation coefficient
`xi` by an asymptotic test. 
"""
function asymptotic_test(ξᵢ::Float64, X::AbstractVector, Y::AbstractVector; rank=denserank)
    # intermediates depending on the tie breaks
    n = length(Y)
    v = 

    # results
    sd = sqrt(v/n)
    pval = 1 - cdf(Normal(), sqrt(n)*ξᵢ/sqrt(v))
    return pval, sd
end
=#


"""
    test_dependence(xi, X, Y; n_perm=1000, noties=false)

Computes a p-value and standard deviation quantifying the significance of the 
correlation coefficient `xi` between `X` and `Y`. If `noties` is `true`, an 
asymptotic test is used. If `noties` is `false`, a permutation test is used if 
there are less than 20 observations, Otherwise, an asymptotic test is used based
on the tie breaks of achieved by the `rank` function.
"""
function test_dependence(ξᵢ::Float64, X::AbstractVector, Y::AbstractVector; n_perm::Int=1000, rank=denserank, noties::Bool=false)
    if noties
        n = length(Y)
        pval = 1 - cdf(Normal(), sqrt(n)*ξᵢ/sqrt(2/5))
        sd = sqrt(2/(5*n))
        return pval, sd
    else#=if=# n <= 20  # perform permutation test
        pval, sd = permutation_test(ξᵢ, Y; n=n_perm)
        return pval, sd
    #=else            # perform asymptotic test
        pval, sd = asymptotic_test(ξᵢ, X, Y; rank=rank)
        return pval, sd
    =#
    end
end
