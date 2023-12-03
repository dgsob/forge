# using InteractiveUtils
# using ProfileView, Profile

# #=
#     This is a direct translation of python's code provided by: 
#     Li, S., Bohman, B., & Jayatilaka, D. (2022). 
#     Enumerating Possible Molecular Formulae in Mass Spectrometry Using a Generating Function Based Method. 
#     MATCH Communications in Mathematical and in Computer Chemistry, 88(2).
#     It is type unstable.  
# =#

# function is_leaf(t)
#     return length(t) == 1
# end

# function sum_tree(t)
#     return is_leaf(t) ? sum(t) : t[1] + sum_tree(t[2])
# end

# function traverse_tree(tree)
#     res = []
#     function s(t, r)
#         if is_leaf(t)
#             push!(res, [r[2:end]; t...])
#         else
#             [s(branch, [r..., t[1]]) for branch in t[2:end]]
#         end
#     end
#     s(tree, [])
#     return res
# end

# using IterTools: groupby

# function multiply_gf(f1, f2, M)
#     r = collect((vcat(j,i) for i in f1 for j in f2 if sum_tree(vcat(j,i)) <= M))
#     sort!(r, by = sum_tree)
#     return collect(groupby(sum_tree, r))
# end

# function multiply_gf(f1, l, M, i)
#     f2 = l[i]
#     r = collect((vcat(j,i) for i in f1 for j in f2 if sum_tree(vcat(j,i)) <= M))
#     # println("R BEFORE: ", r)
#     # p = sortperm(r, alg = QuickSort, by = sum_tree)
#     sort!(r, alg = QuickSort, by = sum_tree)
#     # println("R SORTED: ", r)
#     if i-1 !== 0
#         multiply_gf(collect(groupby(sum_tree, r)), l, M, i-1)
#     else
#         # return collect(groupby(sum_tree, r))
#         return collect(groupby(sum_tree, r))[end]
#     end
# end

# function produce_gfs(alphabet, M)
#     return [((j) for j in 0:i:M) for i in alphabet]
# end

# function formula_tree(alphabet, M)
#     l = produce_gfs(alphabet, M)
#     # prod = map(x->[x], l[end])
#     # for i in reverse(eachindex(l[1:end-1]))
#     #     prod = multiply_gf(prod, l[i], M)
#     # end

#     # prod = multiply_gf(map(x->[x], l[end]), l, M, lastindex(l)-1)
#     # return prod[end]

#     return multiply_gf(map(x->(x), l[end]), l, M, lastindex(l)-1)
# end

# function generate_lst(t)
#     return traverse_tree([0; t])
# end

# function fill_res!(res, lst, compomer)
#     for i in eachindex(compomer)
#         append!(res, fill(lst[i], compomer[i]))
#     end
# end

# function produce_MF(mon, alphabet, lst)
#     # lst = ["C", "H", "Cl", "F", "N", "O", "P", "S"]
#     # lst = ["C", "H", "N", "O", "S"]
#     # lst = [[2,4], [1], [3,5], [2], [2,4,6]]
#     compomer = [div(mon[i], alphabet[i]) for i in eachindex(mon)]
#     # res = [lst[i] * string(compomer[i]) for i in eachindex(compomer) if compomer[i] != 0]
#     # return join([i == "1" ? i[1] : i for i in res])
#     # make it return vectors of valences instead of formulae:
#     # if compomer[1] != 0 && compomer[2] != 0
#     res = Vector{Vector{Int}}()
#     fill_res!(res, lst, compomer)
#     return res
#     # else
#     #     return 
#     # end
# end

# function enumerate_MF(alphabet, lst, M)
#     # @code_warntype formula_tree(alphabet, M)
#     return [produce_MF(i, alphabet, lst) for i in generate_lst(formula_tree(alphabet, M))]
# end