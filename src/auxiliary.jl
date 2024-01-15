#=
    Functions and methods mostly for converting data structures 
=#
"""
Converts input mass to integer mass. It will round it to the nearest integer which is not always correct.
Sometimes the value corresponding to combinatios of all the mass of individual isotopes sum up to a preceise
mass of which integer value should be floor. Pass nearest=false if the results obtained by default are not satisfactory.
"""
function convert_to_ints(M, masses, ϵ, nearest=true)
    M_int = round(Int, M * ϵ)
    if !nearest
        M_int = floor(Int, M * ϵ)
    end
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
