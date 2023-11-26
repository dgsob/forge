# using InteractiveUtils

#=
    This is a direct translation of python's code provided by: 
    Li, S., Bohman, B., & Jayatilaka, D. (2022). 
    Enumerating Possible Molecular Formulae in Mass Spectrometry Using a Generating Function Based Method. 
    MATCH Communications in Mathematical and in Computer Chemistry, 88(2).
    It needs some Julia-specific optimization. Notably it is type unstable.  
=#

function is_leaf(t)
    return length(t) == 1
end

function sum_tree(t)
    return is_leaf(t) ? sum(t) : t[1] + sum_tree(t[2])
end

function traverse_tree(tree)
    res = []
    
    function s(t, r)
        if is_leaf(t)
            push!(res, r[2:end] .+ t)
        else
            [s(branch, [r; t[1]]) for branch in t[2:end]]
        end
    end
    
    s(tree, [])
    return res
end

using IterTools: groupby

function multiply_gf(f1, f2, M)
    r = []
    
    for i in f1
        for j in f2
            if sum_tree(vcat(j,i)) <= M
                push!(r, vcat(j,i))
            end
        end
    end
    sort!(r, by=sum_tree)
    
    result = []
    
    for g in groupby(sum_tree, r)
        push!(result, collect(g))
    end
    
    return result
end

function produce_gfs(alphabet, M)
    return [[[j] for j in 0:i:M] for i in alphabet]
end

function formula_tree(alphabet, M)
    l = produce_gfs(alphabet, M)
    prod::Vector{Any} = [[i] for i in l[end]]
    for i in length(l)-2:-1:1
        prod = multiply_gf(prod, l[i], M)
    end
    return prod[end]
end

function generate_lst(t)
    return traverse_tree(vcat([0], t))
end

function produce_formula(mon, alphabet, lst::Vector{String})
    compomer = [div(mon[i], alphabet[i]) for i in eachindex(mon)]
    res = [lst[i] * (compomer[i] == 1 ? "" : string(compomer[i])) for i in eachindex(compomer) if compomer[i] != 0]
    return join([i == "1" ? i[1] : i for i in res])
end

function fill_res!(res, lst, compomer)
    for i in eachindex(compomer)
        append!(res, fill(lst[i], compomer[i]))
    end
end

function produce_formula(mon, alphabet, lst::Vector{Vector{Int}})
    compomer = [div(mon[i], alphabet[i]) for i in eachindex(mon)]
    res = Vector{Vector{Int}}()
    fill_res!(res, lst, compomer)
    return res
end

"""
    aufbau_generator(alphabet, lst, M) -> formulae
    Takes masses as alphabet and either symbols or valences as lst.
    Returns either a vector of strings representing formulae, e.g. ["C2H6", ...]
    or a vector of valences (vector of vectors of ints) representing formulae, 
    e.g. [[[2, 4], [2, 4], [1], [1], [1], [1], [1], [1]], ...].
"""
function aufbau_generator(alphabet, lst, M)
    return [produce_formula(i, alphabet, lst) for i in generate_lst(formula_tree(alphabet, M))]
end

# println(aufbau_generator(masses, valences, 30))

# @code_warntype aufbau_generator(masses, valences, 30)
