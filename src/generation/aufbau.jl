using InteractiveUtils
using ProfileView, Profile
using StaticArrays
using BenchmarkTools
include("./filtering/basic_organic.jl")
include("./filtering/mass_precision.jl")
include("./filtering/mgraph.jl")

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

function ∑(t)
    return isa(t, Int) ? t : t.leaf + ∑(t.subbranch)
end

function ∑(t::Vector{Branch})
    return t[1].leaf + ∑(t[1].subbranch)
end

function build_repeat_seq(compomer, degrees)
    result = Vector{Vector{Int}}(undef, sum(compomer))
    index = 1
    for i in eachindex(compomer)
        if compomer[i] === 0
            continue
        end
        for j in 1:compomer[i]
            result[index] = degrees[i]
            index += 1
        end
    end
    return result
end

function sc(t::Int, r, res, alphabet, masses_precise, M_precise, ϵ, valences)
    compomer = [r[2:end]; t] .÷ alphabet
    if basic_organic_filter(compomer) && mass_precision_filter(compomer, masses_precise, M_precise, ϵ) && mgraph_filter(build_repeat_seq(compomer, valences))
        push!(res, compomer)
    end
end

function sc(t::Branch, r, res, alphabet, masses_precise, M_precise, ϵ, valences)
    sc(t.subbranch, [r..., t.leaf], res, alphabet, masses_precise, M_precise, ϵ, valences)
end

function sc(t::Vector{Branch}, r, res, alphabet, masses_precise, M_precise, ϵ, valences)
    for branch in t
        sc(branch.subbranch, [r..., branch.leaf], res, alphabet, masses_precise, M_precise, ϵ, valences)
    end
end

function traverse_tree(tree, alphabet, masses_precise, M_precise, ϵ, valences)
    res = Vector{Vector{Int}}()
    sc(tree, Int[], res, alphabet, masses_precise, M_precise, ϵ, valences)
    return res
end

using IterTools: groupby

function multiply_gf(f1, l, M, i)
    f2 = l[i]
    r = collect((Branch(j,i...) for i in f1 for j in f2 if ∑(Branch(j,i...)) <= M))
    sort!(r, alg = QuickSort, by = ∑)
    if i-1 !== 0
        multiply_gf(collect(groupby(∑, r)), l, M, i-1)
    else
        return collect(groupby(∑, r))[end]
    end
end

function produce_gfs(alphabet, M)
    return collect((j for j in 0:i:M) for i in alphabet)
end

function formula_tree(alphabet, M)
    l = produce_gfs(alphabet, M)
    return multiply_gf(l[end], l, M, lastindex(l)-1)
end

function generate_lst(t, alphabet, masses_precise, M_precise, ϵ, valences)
    return traverse_tree(Branch(0, t...), alphabet, masses_precise, M_precise, ϵ, valences)
end

function enumerate_MF(alphabet, M, masses_precise, M_precise, ϵ, valences)
    list = Vector{Vector{Int64}}()
    list = generate_lst(formula_tree(alphabet, M), alphabet, masses_precise, M_precise, ϵ, valences)
    return list
end