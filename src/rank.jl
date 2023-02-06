using Random, StatsBase

"""
A function that returns the sample ranks using maximum ranks 
for equal values.
    - x : A Vector of Floats, Integers, or Bools
"""
function maxrank(x)
    sorted_x = sort(x)
    ranks = zeros(Int, length(x))
    for (idx,i) in enumerate(x)
        ranks[idx] = findlast(isequal(i), sorted_x)
    end
    return ranks
end


"""
A function that returns the sample ranks using minimum ranks 
for equal values.
    - x : A Vector of Floats, Integers, or Bools
"""
function minrank(x)
    sorted_x = sort(x)
    ranks = zeros(Int, length(x))
    for (idx,i) in enumerate(x)
        ranks[idx] = findfirst(isequal(i), sorted_x)
    end
    return ranks
end


"""
A function that returns the sample ranks using random ranking 
of equal values.
    - x : A Vector of Floats, Integers, or Bools
    - seed : either an Integer to seed the random number 
             generator or `missing` to use the default seed.
"""
function randrank(x; seed=missing)
    @assert typeof(seed) <: Integer || typeof(seed) <: Missing
    if !ismissing(seed)
        Random.seed!(seed)
    end
    sorted_x = sort(x)
    ranks = zeros(Int, length(x))
    for (idx,i) in enumerate(x)
        ties = findall(isequal(i), sorted_x)
        remaining_ranks = [r for r in ties if r ∉ ranks]
        ranks[idx] = rand(remaining_ranks)
    end
    return ranks
end


"""
A function that returns the sample ranks of the values in a vector.
    - x : A Vector of Floats, Integers, or Bools
    - ties_method : a String indicating how equal values (i.e. ties) should be
                    handled. 
TO DO: the original R `rank` function has an argument that let's you choose the
       placement of the missing values.
"""
function rank(x; ties_method::String="average", seed=missing)
    @assert ties_method ∈ ["average", "random", "max", "min"]
    @assert typeof(seed) <: Integer || typeof(seed) <: Missing

    # get NA/missing
    missing_values = ismissing.(x)
    # rank based on the correct ties method
    if ties_method == "average"
        ranks = tiedrank(x)
    elseif ties_method == "random"
        ranks = randrank(x; seed=seed)
    elseif ties_method == "max"
        ranks = maxrank(x)
    elseif ties_method == "min"
        ranks = minrank(x)
    end
    # add ranks of the missing values (if any) as the last ones
    if any(missing_values)
        counter = 1
        for rdx in 1:lastindex(ranks)
            if missing_values[rdx]
                ranks[rdx] = length(ranks)+counter
                counter += 1
            end
        end
    end
    return ranks
end