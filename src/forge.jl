include("./elements.jl")
include("./generation/aufbau.jl")
include("./generation/naive.jl")
include("./filtering/mgraph.jl")
using BenchmarkTools

function collect_formulae!(real_strings, real_vectors, symbols, order)
    for compound in real_vectors
        atom_count = Dict{String, Int}()
        for atom in compound
            symbol = symbols[atom]
            if haskey(atom_count, symbol)
                atom_count[symbol] += 1
            else
                atom_count[symbol] = 1
            end
        end
        formula = ""
        for symbol in sort(collect(keys(atom_count)), by = x -> findfirst(isequal(x), order))
            formula *= symbol
            if atom_count[symbol] > 1 # so that 1 is ommited as in chemical formulas
                formula *= "$(atom_count[symbol])"
            end
        end    
        push!(real_strings, formula)
    end
end

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

function main(M_precise, ϵ, symbols, s_dict, v_dict, valences, masses_precise, generator = "aufbau")
    # Set-up
    println("Set-up starts")
    M, masses = convert_to_ints(M_precise, masses_precise, ϵ)
    println("Total mass:", M)
    for i in eachindex(symbols)
        println(symbols[i], " ", masses[i])
    end
    println("Set-up completed")
    println("")

    # Stage 1: building up all formulaecollect_formulae!(formulae_to_display)
    all_vectors = Vector{Vector{Int64}}[]
    if generator == "aufbau"
        all_vectors = enumerate_MF(masses, valences, M)
    elseif generator == "naive"
        all_vectors = naive_generator(v_dict, masses, symbols, M, 0)
    end
    println("All vectors generated, L: ", length(all_vectors))
    if length(all_vectors) < 4
        println("Vectors: ", all_vectors)
    end
    println("")

    # Stage 2: filtering 
    real_vectors = Vector{Vector{Int64}}[]
    fill_real_vectors!(real_vectors, all_vectors)
    println("Filtering completed, L: ", length(real_vectors))
    if length(real_vectors) < 4
        println("Filtered vectors: ", real_vectors)
    end
    println("")

    # Display formulae 
    real_strings = Vector{String}()
    collect_formulae!(real_strings, real_vectors, s_dict, symbols)
    println("Vectors converted to formulae: ", real_strings)
    println("")
end




main(103.12100, 1e5, element_symbols, symbols_dict, valences_dict, element_valences, element_masses, "naive")