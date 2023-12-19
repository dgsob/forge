include("../src/elements.jl")
include("../src/generation/aufbau.jl")
include("../src/generation/naive.jl")
include("../src/filtering/mgraph.jl")
include("../src/auxiliary.jl")
include("../src/filtering/basic_organic.jl")
using BenchmarkTools

function fill_realizables!(realizables, compomers, valences)
    for compomer in compomers
        if !basic_organic_filter(compomer)
            continue
        end
        if !mgraph_filter(build_repeat_seq(compomer, valences))
            continue
        end
        push!(realizables, compomer)
    end
end

function main(M_precise, ϵ, symbols, valences, masses_precise)
    # Set-up
    println("Set-up starts")
    M, masses = convert_to_ints(M_precise, masses_precise, ϵ)
    println("Total mass:", M)
    for i in eachindex(symbols)
        println(symbols[i], ":", masses[i])
    end
    println("Set-up completed")
    println("")

    # Stage 1: building up all formulae
    compomers = Vector{Vector{Int64}}[]
    println("Generation time: ")
    @time compomers = enumerate_MF(masses, M)
    println("Compomers generated, L: ", length(compomers))
    if length(compomers) < 50
        println("Compomers: ", compomers)
    end
    println("")
    
    # Stage 2: filtering 
    realizables = Vector{Int}[]
    println("Filtration time: ")
    @time fill_realizables!(realizables, compomers, valences)
    n = length(realizables)
    println("Filtering completed, L: ", n)
    if n < 20
        println("Filtered vectors: ", realizables)
    end
    println("")

    # Display formulae, only first 30 if larger
    if n ≤ 30
        display_formulae(realizables, symbols, n)
    else
        display_formulae(realizables[1:30], symbols, 30)
    end
end
# 103.12100 ; 46 ; 775
main(326.10016, 1e3, element_symbols, element_valences, element_masses)