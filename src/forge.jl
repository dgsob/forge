include("./elements.jl")
include("./generation/aufbau.jl")
include("./auxiliary.jl")
using BenchmarkTools

function main(M_precise, ϵ, symbols, valences, masses_precise)
    # Set-up
    println("Set-up starts")
    M_int, masses_int = convert_to_ints(M_precise, masses_precise, 1)
    println("Total mass:", M_int, ", ", M_precise, ", ", M_precise+ϵ, ", ", M_precise-ϵ)
    for i in eachindex(symbols)
        println(symbols[i], ":", masses_int[i])
    end
    println("Set-up completed")
    println("")

    # Building up all formulae + filtration
    compomers = Vector{Vector{Int64}}[]
    @time compomers = enumerate_MF(masses_int, M_int, masses_precise, M_precise, ϵ, valences)
    println("Compomers generated, L: ", length(compomers))

    n = length(compomers)
    if n ≤ 30
        println("Compomers: ", compomers)
    end
    println("")

    # Display formulae, only first 30 if larger
    if n ≤ 30
        display_formulae(compomers, symbols, n)
    else
        display_formulae(compomers[1:30], symbols, 30)
    end
end
# 103.12100 ; 46 ; 775
main(106.00885, 5e-5, element_symbols, element_valences, element_masses)