module Xicor

using Distributions, Random, StatsBase

"""
    xicor(X, Y; rank=tiedrank, noties=false)

Computes the asymmetric ξ correlation coefficient sbetween two vectors `X` and `Y`.
`rank` species how the rank is computed, default uses `tiedrank`. If there
are no ties in `Y`, `noties` can be set to `true`, which speeds up computation.

`ξ` is an alias for this function.
"""
function xicor(X::AbstractVector, Y::AbstractVector; rank=tiedrank, noties=false)
    @assert length(X) == length(Y) "length `X` and `Y` missmatch"
    n = length(Y)
    xind = sortperm(X)
    r = rank(xind, by=i->Y[i])
    if !noties
        l = rank(xind, by=i->Y[i], rev=true)
    end
    ξn = 0.0
    for i in 1:n-1
        ξn += (abs(r[i+1] - r[i]))
    end
    ξn *= noties ? 3 / (n^2 - 1) : n / 2sum(li->li*(n-li), l)
    ξn = 1 - ξn
    return ξn
end

ξ = xicor


include("correlation.jl")
include("rank.jl")

export xicor, ξ

end
