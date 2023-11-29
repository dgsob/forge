using OffsetArrays
using InteractiveUtils

function generate_combinations(masses, elements::Vector{String}, M, ϵ)
    d = OffsetArray(fill(M + 1, 0:M), 0)
    d[0] = 0
    for i in 1:M
        for c in masses
            if i - c >= 0
                d[i] = min(d[i], d[i - c] + 1)
            end
        end
    end
    combinations = Vector{Vector{String}}()
    find(M, [], d, masses, combinations, elements, 1, ϵ)
    return combinations
end

function find(M, result, d, masses, combinations, elements, index, ϵ)
    if abs(M) <= ϵ
        push!(combinations, result)
        return
    end
    for i in index:length(masses)
        if M - masses[i] >= -ϵ
            find(M - masses[i], [result; elements[i]], d, masses, combinations, elements, i, ϵ)
        end
    end
end
# # masses = [12,1,35,19,14,16,31,32]
# masses = [12,1,14,16,32]
# # symbols = ["C", "H", "Cl", "F", "N", "O", "P", "S"]
# symbols = ["C", "H", "N", "O", "S"]
# include("../elements.jl")
# valences = valences_dict
# M = 43
# ϵ = 0

function fill_valences!(vc, v, sc)
    for i in eachindex(sc)
        for j in eachindex(sc[i])
            vc[i][j] = v[sc[i][j]]
        end
    end
end

function naive_generator(valences, masses, symbols, M, ϵ)
    combinations = generate_combinations(masses, symbols, M, ϵ)
    valence_combinations = [Vector{Vector{Int}}(undef, length(c)) for c in combinations]
    fill_valences!(valence_combinations, valences, combinations)
    return valence_combinations
end

# function main(v, ms, s, _M, _ϵ)
#     @time println(naive_generator(v, ms, s, _M, _ϵ))
# end

# main(valences, masses, symbols, M, ϵ)
