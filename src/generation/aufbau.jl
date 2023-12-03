using InteractiveUtils
using ProfileView, Profile
using StaticArrays
using BenchmarkTools

#=
    This is based on the python's code provided by: 
    Li, S., Bohman, B., & Jayatilaka, D. (2022). 
    Enumerating Possible Molecular Formulae in Mass Spectrometry Using a Generating Function Based Method. 
    MATCH Communications in Mathematical and in Computer Chemistry, 88(2).
=#

struct Branch <: Any
    leaf::Int
    subbranch::Union{Int, Branch, Vector{Branch}}
end

function Branch(l::Int, b::Branch, bb::Branch, args...)
    return Branch(l, [b, bb, args...])
end

function sum_tree(t)
    return isa(t, Int) ? t : t.leaf + sum_tree(t.subbranch)
end

function sum_tree(t::Vector{Branch})
    return t[1].leaf + sum_tree(t[1].subbranch)
end

function collect_compomers(t::Int, r, res, alphabet)
    push!(res, [r[2:end]; t] .÷ alphabet)
end

function collect_compomers(t::Branch, r, res, alphabet)
    collect_compomers(t.subbranch, [r..., t.leaf], res, alphabet)
end

function collect_compomers(t::Vector{Branch}, r, res, alphabet)
    for branch in t
        collect_compomers(branch.subbranch, [r..., branch.leaf], res, alphabet)
    end
end

function traverse_tree(tree, alphabet)
    res = Vector{Vector{Int}}()
    collect_compomers(tree, Int[], res, alphabet)
    return res
end

using IterTools: groupby

function multiply_gf(f1, l, M, i)
    f2 = l[i]
    r = collect((Branch(j,i...) for i in f1 for j in f2 if sum_tree(Branch(j,i...)) <= M))
    sort!(r, alg = QuickSort, by = sum_tree)
    if i-1 !== 0
        multiply_gf(collect(groupby(sum_tree, r)), l, M, i-1)
    else
        return collect(groupby(sum_tree, r))[end]
    end
end

function produce_gfs(alphabet, M)
    return [(j for j in 0:i:M) for i in alphabet]
end

function formula_tree(alphabet, M)
    l = produce_gfs(alphabet, M)
    return multiply_gf(l[end], l, M, lastindex(l)-1)
end

function generate_lst(t, alphabet)
    return traverse_tree(Branch(0, t...), alphabet)
end

# function fill_res!(res, lst, compomer)
#     for i in eachindex(compomer)
#         append!(res, fill(lst[i], compomer[i]))
#     end
# end

# function produce_MF(mon, alphabet, lst)
#     compomer = [div(mon[i], alphabet[i]) for i in eachindex(mon)]

#     # res = [lst[i] * string(compomer[i]) for i in eachindex(compomer) if compomer[i] != 0]
#     # return join([i == "1" ? i[1] : i for i in res])

#     # res = Vector{Vector{Int}}()
#     # fill_res!(res, lst, compomer)
#     # return res
# end

function enumerate_MF(alphabet, lst, M)
    list = Vector{Vector{Int64}}()
    list = generate_lst(formula_tree(alphabet, M), alphabet)
    return list
    # return collect(produce_MF(i, alphabet, lst) for i in list)
end

# using BenchmarkTools
# # @time enumerate_MF([12,1,35,19,14,16,31,32],500)
# function main()
#     enum = enumerate_MF([12,1,14,16,32],43)
#     # ProfileView.@profview enumerate_MF([12,1,14,16],133)
#     println(length(enum), " ", enum)
#     # ben = @benchmark enumerate_MF([12,1,14,16],133)
#     # display(ben)
#     # @btime enumerate_MF([12,1,14,16,32],900)
#     # @time enumerate_MF([12,1,35,19,14,16,31,32],43)
# end

# main()