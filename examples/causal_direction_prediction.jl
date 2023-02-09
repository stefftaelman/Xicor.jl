include("../src/Xicor.jl")
#using Xicor

#=
We use the Tuebingen cause-effect pairs database from Mooij et al. (2017) 
to test the potential of Xi correlation as a causal direction prediction 
method. The data can be downloaded from 
https://webdav.tuebingen.mpg.de/cause-effect/pairs.zip, and consists of a
number of txt files, each containing a list of cause-effect pairs and 
description files containing information about the data sources and the 
ground truth direction of the causal relationships.
=#


"""
A function to load a pair of variables X and Y from the Tuebingen 
cause-effect pairs database.  The data can be downloaded from 
https://webdav.tuebingen.mpg.de/cause-effect/pairs.zip, and consists of a
number of txt files, each containing a list of cause-effect pairs and 
description files containing information about the data sources and the 
ground truth direction of the causal relationships.
    - pair_id: the Integer id of the pair to load
    - path_to_data_dir: the path to the directory containing the data files
"""
function load_pair(pair_id::Integer; path_to_data_dir="./data")
    pair_value_file = path_to_data_dir * "/pair" * lpad(string(pair_id), 4, '0') * ".txt"
    pair_values = hcat(
        [
            parse.(
                Float64, 
                [q for q in split(i, r"\s+") if q != ""]
                ) for i in split(chomp(readchomp(pair_value_file)), "\n")
        ]...
        )
    X = pair_values[1, :]
    Y = pair_values[2, :]
    return X, Y
end


"""
A function to get the ground truth direction of a pair of variables X and Y
based on the description file associated with the pair. The directions for 
pair 86 and 88 were hard-coded, as the description files do not explicitly
state the ground truth (arrow).
    - pair_id: the Integer id of the pair to load
    - path_to_data_dir: the path to the directory containing the data files
"""
function get_ground_truth_direction(pair_id::Integer; path_to_data_dir="./data")
    # Hard-code the ground truth direction for pairs 86 and 88
    if pair_id == 86 || pair_id == 88
        return "X->Y"
    end

    pair_descr_file = path_to_data_dir * "/pair" * lpad(string(pair_id), 4, '0') * "_des.txt"
    pair_descr = read(pair_descr_file, String)
    xy_ground_truth_text = match(
        r"[Xx]\s*(<{0,1}\s*(-\s*)+>{0,1})\s*[Yy]", 
        pair_descr
        )
    if isnothing(xy_ground_truth_text)
        yx_ground_truth_text = match(
            r"[Yy]\s*(<{0,1}\s*(-\s*)+>{0,1})\s*[Xx]", 
            pair_descr
            )
        if isnothing(yx_ground_truth_text)
            @error "No ground truth direction found for pair $(pair_id)"
        else
            return occursin('<', yx_ground_truth_text[1]) ? "X->Y" : "Y->X"
        end
    else
        return occursin('>', xy_ground_truth_text[1]) ? "X->Y" : "Y->X"
    end
end


"""
A function to predict the causal direction of a pair of variables X and Y
based on the Xi correlation.
    - X: a vector of values for variable X
    - Y: a vector of values for variable Y
    - alpha: the significance level for the Xi correlation test
    - nperm: the number of permutations to use for the Xi correlation test
"""
function predict_direction(X, Y; alpha=0.05, method="permutation", nperm=1000)
    xy_xi, xy_sd, xy_p = Xicor.xicor(
        X, Y; 
        pvalue=true, method=method, nperm=nperm
        )
    yx_xi, yx_sd, yx_p = Xicor.xicor(
        Y, X; 
        pvalue=true, method=method, nperm=nperm
        )

    if xy_p < alpha && yx_p < alpha
        if abs(xy_xi) > abs(yx_xi)
            return "X->Y"
        elseif abs(xy_xi) < abs(yx_xi)
            return "Y->X"
        elseif xy_sd > yx_sd
            return "X->Y"
        elseif xy_sd < yx_sd
            return "Y->X"
        else
            return "uncertain"
        end
    elseif xy_p < alpha
        return "X->Y"
    elseif yx_p < alpha
        return "Y->X"
    else
        return "uncertain"
    end
end


### check performance
TP = 0
FP = 0
TN = 0
FN = 0
for pair_id in 1:108
    # get data
    X, Y = load_pair(pair_id)
    dir_truth = get_ground_truth_direction(pair_id)

    # predict direction
    dir_pred = predict_direction(X, Y; method="asymptotic") # this can take a while
    if dir_pred == "uncertain"
        println("Pair $(pair_id) is uncertain")
    else
        if dir_pred == dir_truth
            if dir_pred == "X->Y"
                TP += 1
            else
                TN += 1
            end
        else
            if dir_pred == "X->Y"
                FP += 1
            else
                FN += 1
            end
        end
    end
end

println("TP: $(TP)")
println("FP: $(FP)")
println("TN: $(TN)")
println("FN: $(FN)")
println("Accuracy: $(100*(TP+TN)/(TP+FP+TN+FN))%")
println("Precision: $(100*TP/(TP+FP))%")
println("Recall: $(100*TP/(TP+FN))%")
println("F1: $(100*2TP/(2TP+FP+FN))%")