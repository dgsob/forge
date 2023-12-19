#=
    Functions and methods mostly for converting data structures 
=#

function convert_to_ints(M, masses, ϵ)
    # TODO: Figure out a way to round the mass values to integers in a right way
    M_int = round(Int, M * ϵ)
    # M_int = floor(Int, M * ϵ)
    masses_int = round.(Int, masses .* ϵ)
    return M_int, masses_int
end

# function build_repeat_dict(compomer, degrees, symbols)
#     return Dict(zip((symbols, degrees), collect(i for i in compomer if i !== 0)))
# end

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

function fill_formulae!(formulae, compomers, symbols)
    index = 1
    for compomer in compomers
        formula = ""
        for i in eachindex(compomer)
            if compomer[i] === 0
                continue
            end
            formula *= symbols[i]
            if compomer[i] > 1
                formula *= "$(compomer[i])"
            end
        end
        formulae[index] = formula
        index += 1
    end
end

function display_formulae(realizables, symbols, size)
    formulae = Vector{String}(undef, size)
    fill_formulae!(formulae, realizables, symbols)
    println("Compomers converted to formulae: ", formulae)
    println("")
end
