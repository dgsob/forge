using InteractiveUtils

#=
    This is a direct translation of python's code provided by: 
    Li, S., Bohman, B., & Jayatilaka, D. (2022). 
    Enumerating Possible Molecular Formulae in Mass Spectrometry Using a Generating Function Based Method. 
    MATCH Communications in Mathematical and in Computer Chemistry, 88(2).
    It is type unstable.  
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
            push!(res, [r[2:end]; t...])
        else
            [s(branch, [r..., t[1]]) for branch in t[2:end]]
        end
    end
    s(tree, [])
    return res
end

using IterTools: groupby

function multiply_gf(f1, f2, M)
    r = [vcat(j,i) for i in f1 for j in f2 if sum_tree(vcat(j,i)) <= M]
    sort!(r, by = sum_tree)
    return [collect(g) for g in groupby(sum_tree, r)]
end

function produce_gfs(alphabet, M)
    return [[[j] for j in 0:i:M] for i in alphabet]
end

function formula_tree(alphabet, M)
    l = produce_gfs(alphabet, M)
    prod = [[i] for i in l[end]]
    for i in reverse(1:length(l)-1)
        prod = multiply_gf(prod, l[i], M)
    end
    return prod[end]
end

function generate_lst(t)
    return traverse_tree([0; t])
end

function fill_res!(res, lst, compomer)
    for i in eachindex(compomer)
        append!(res, fill(lst[i], compomer[i]))
    end
end

function produce_MF(mon, alphabet)
    lst = ["C", "H", "Cl", "F", "N", "O", "P", "S"]
    compomer = [div(mon[i], alphabet[i]) for i in 1:length(mon)]
    res = [lst[i] * string(compomer[i]) for i in 1:length(compomer) if compomer[i] != 0]
    return join([i == "1" ? i[1] : i for i in res])
end

function enumerate_MF(alphabet, M)
    return [produce_MF(i, alphabet) for i in generate_lst(formula_tree(alphabet, M))]
end

using BenchmarkTools
# @time enumerate_MF([12,1,35,19,14,16,31,32],500)
@time println(length(enumerate_MF([12,1,35,19,14,16,31,32],775)))
