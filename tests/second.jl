include("../src/elements.jl")
include("../src/generation/aufbau.jl")
include("../src/generation/filtering/mgraph.jl")
include("../src/auxiliary.jl")

include("../src/elements.jl")
include("../src/generation/aufbau.jl")
include("../src/generation/filtering/mgraph.jl")
include("../src/auxiliary.jl")
include("./old_m_versions.jl")

using BenchmarkTools

MM = 180
function arts!(realizables, compomers, valences)
    for compomer in compomers
        seq = build_repeat_seq(compomer, valences)
        # if length(seq) < 60
        #     continue
        # end
        if mgraph_filter(seq)
            push!(realizables, compomer)
            continue
        end
        measure = 1
        for i in seq
            measure *= length(i)
        end
        if measure < 0
            println("Something's wierd")
            continue
        end
        measure = measure^(1/length(seq))
        
        if measure <= 1.5
            continue
        end
        println("-----------------------------------")
        println(seq)
        println("Measure: ")
        println(measure, ", ", length(seq))
        println("M filtering timing starts...")
        @btime mgraph_filter($seq)
        println("M filtering timing completed.")
        # println("Brute force timing starts...")
        # @btime check_sequence_iter($seq)
        # println("Brute force timing completed.")
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
    compomers = enumerate_MF(masses, M, masses_precise, M_precise, ϵ, valences)
    # @btime enumerate_MF($masses, $M)
    println("Compomers generated, L: ", length(compomers))
    if length(compomers) < 5
        println("Compomers: ", compomers)
    end
    println("")
    
    # Stage 2: filtering 
    realizables_arts = Vector{Int}[]
    arts!(realizables_arts, compomers, valences)
    println("L_ARTS: ", length(realizables_arts))
end
# 103.12100 ; 46 ; 775
main(MM, 1, element_symbols, element_valences, element_masses)