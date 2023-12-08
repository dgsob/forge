f1 = Vector{Vector{Any}}[[[0, Any[0, Any[0, [0, 0]]]]], [[100800, Any[0, Any[0, [0, 0]]]]], [[201600, Any[0, Any[0, [0, 0]]]]], [[302400, Any[0, Any[0, [0, 0]]]]], [[403200, Any[0, Any[0, [0, 0]]]]], [[504000, Any[0, Any[0, [0, 0]]]]], [[604800, Any[0, Any[0, [0, 0]]]]], [[705600, Any[0, Any[0, [0, 0]]]]], [[806400, Any[0, Any[0, [0, 0]]]]], [[907200, Any[0, Any[0, [0, 0]]]]], [[1008000, Any[0, Any[0, [0, 0]]]]], [[1108800, Any[0, Any[0, [0, 0]]]]], [[1209600, Any[0, Any[0, [0, 0]]]]], [[1310400, Any[0, Any[0, [0, 0]]]]], [[0, Any[1400700, Any[0, [0, 0]]]]]]
f2 = [0, 1201100]
M = 1402700

function is_leaf(t)
    return length(t) == 1
end

function sum_tree(t)
    return is_leaf(t) ? sum(t) : t[1] + sum_tree(t[2])
end

function multiply_gf(f1, f2)
    r = collect((vcat(j,i) for i in f1 for j in f2 if sum_tree(vcat(j,i)) <= M))
    sort!(r, by = sum_tree)
    println(r)
end

multiply_gf(f1, f2)