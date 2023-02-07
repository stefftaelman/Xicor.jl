using StatsBase, Distributions
include("rank.jl")

"""
A function that returns the Xi correlation for two vectors xvec and yvec.
"""
function calc_xi(xvec, yvec; simple::Bool=true, seed=missing)
    @assert typeof(seed) <: Integer || typeof(seed) <: Missing

    n = length(xvec)
    PI = rank(xvec; ties_method="random", seed=seed)
    fr = rank(yvec; ties_method="max")./n
    gr = rank(-1 .* yvec; ties_method="max")./n 
    ord = sortperm(PI)
    fr = fr[ord]
    A1 = sum(abs.(fr[1:(n - 1)] - fr[2:n]))/(2 * n)
    CU = mean(gr .* (1 .- gr))
    xi = 1 - A1/CU
    if simple
        return xi
    else
        return xi, fr, CU
    end
end


"""
This function computes the xi coefficient between two vectors x and y.
It can be used to test independence using a Monte Carlo permutation 
test or through an asymptotic approximation test.
Input:
    - xvec :    A Vector of Floats, Integers, or Bools
    - yvec :    A Vector of Floats, Integers, or Bools
    - pvalue :  A Bool indicating whether to compute the p-value
    - ties :    A Bool indicating whether to use ties in the calculation
    - method :  A String indicating the method to use for the p-value
                calculation. Options are "asymptotic" and "permutation".
    - nperm :   An Integer indicating the number of permutations to use
                for the permutation test.
    - seed :    either an Integer to seed the random number generator 
                or `missing` to use the default seed.
"""
function xicor(xvec, yvec; 
               pvalue::Bool=false, ties::Bool=true, 
               method::String="asymptotic", nperm::Integer=1000, 
               seed=missing)
    @assert typeof(seed) <: Integer || typeof(seed) <: Missing

    # TO DO: check if dataframe type and extract vector if so

    # calculate xi corr
    xi, fr, CU = calc_xi(xvec, yvec; simple=false, seed=seed)
    n = length(xvec)

    # testing
    if pvalue
        if !ties
            sd = sqrt(2/(5 * n))
            pval = 1 - cdf(Normal(), sqrt(n) * xi/sqrt(2/5))
            return xi, sd, pval
        end
        @assert method ∈ ["asymptotic", "permutation"]
        if method == "asymptotic"
            # intermediates
            sorted_fr = sort(fr)
            ind = 1:n
            ind2 = 2 * n .- 2 .* ind .+ 1
            ai = mean(ind2 .* sorted_fr .* sorted_fr)/n
            ci = mean(ind2 .* sorted_fr)/n
            cq = cumsum(sorted_fr)
            m = (cq .+ (n .- ind) .* sorted_fr)/n
            b = mean(m.^2)
            v = (ai - 2 * b + ci^2)/(CU^2)
            # results 
            sd = sqrt(v/n)
            pval = 1 - cdf(Normal(), sqrt(n) * xi/sqrt(2/5))
            return xi, sd, pval
        elseif method == "permutation"
            rp = zeros(nperm)
            for i in 1:nperm
                x1 = rand(n)
                rp[i] = calc_xi(x1, yvec)
            end
            sd = sqrt(var(rp))
            pval = mean(rp .> xi)
            return xi, sd, pval
        end
    end
    return xi
end
