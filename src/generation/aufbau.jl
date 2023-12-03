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

function ∑(t)
    return isa(t, Int) ? t : t.leaf + ∑(t.subbranch)
end

function ∑(t::Vector{Branch})
    return t[1].leaf + ∑(t[1].subbranch)
end

function sc(t::Int, r, res, alphabet)
    push!(res, [r[2:end]; t] .÷ alphabet)
end

function sc(t::Branch, r, res, alphabet)
    sc(t.subbranch, [r..., t.leaf], res, alphabet)
end

function sc(t::Vector{Branch}, r, res, alphabet)
    for branch in t
        sc(branch.subbranch, [r..., branch.leaf], res, alphabet)
    end
end

function traverse_tree(tree, alphabet)
    res = Vector{Vector{Int}}()
    sc(tree, Int[], res, alphabet)
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

function generate_lst(t, alphabet)
    return traverse_tree(Branch(0, t...), alphabet)
end

function enumerate_MF(alphabet, M)
    list = Vector{Vector{Int64}}()
    list = generate_lst(formula_tree(alphabet, M), alphabet)
    return list
end