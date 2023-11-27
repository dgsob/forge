include("./elements.jl")
include("./generation/aufbau.jl")
include("./generation/naive.jl")
include("./filtering/mgraph.jl")
using BenchmarkTools

function collect_formulae!(real_strings, real_vectors, all_strings, all_vectors)
    index = 1
    for i in eachindex(all_vectors)
        if all_vectors[i] in real_vectors
            real_strings[index] = all_strings[i]
            index += 1
        end
    end

end

# function collect_formulae!(real_strings, real_vectors, symbols, order)
#     for compound in real_vectors
#         atom_count = Dict{String, Int}()
#         for atom in compound
#             symbol = symbols[atom]
#             if haskey(atom_count, symbol)
#                 atom_count[symbol] += 1
#             else
#                 atom_count[symbol] = 1
#             end
#         end
#         formula = ""
#         for symbol in sort(collect(keys(atom_count)), by = x -> findfirst(isequal(x), order))
#             formula *= symbol
#             if atom_count[symbol] > 1 # so that 1 is ommited as in chemical formulas
#                 formula *= "$(atom_count[symbol])"
#             end
#         end    
#         push!(real_strings, formula)
#     end
# end

function fill_real_vectors!(real_vectors, all_vectors)
    for formula in all_vectors
        if !mgraph_filter(formula)
            continue
        end
        push!(real_vectors, formula)
    end
end

function convert_to_ints(M, masses, ϵ)
    M_int = round(Int, M * ϵ)
    masses_int = round.(Int, masses .* ϵ)
    return M_int, masses_int
end

function main(M_precise, ϵ, symbols, valences, masses_precise)
    # Set-up
    println("Set-up starts")
    M, masses = convert_to_ints(M_precise, masses_precise, ϵ)
    println("Set-up completed")

    # Stage 1: building up all formulaecollect_formulae!(formulae_to_display)
    all_vectors, all_strings = aufbau_generator(masses, valences, symbols, M)
    println("All formulae generated, L: ", length(all_vectors))

    # Stage 2: filtering 
    real_vectors = Vector{Vector{Int64}}[]
    fill_real_vectors!(real_vectors, all_vectors)

    # Display formulae 
    real_strings = Vector{String}(undef, length(real_vectors))
    collect_formulae!(real_strings, real_vectors, all_strings, all_vectors)
    println("Only graphical formulae left, L: ", real_strings)
end




main(130.14301, 1e5, element_symbols, element_valences, element_masses)