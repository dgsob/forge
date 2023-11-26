include("./elements.jl")
include("./generation/aufbau.jl")
include("./filtering/mgraph.jl")

function convert_to_ints(M, masses, ϵ)
    M_int = round(Int, M * ϵ)
    masses_int = round.(Int, masses .* ϵ)
    return M_int, masses_int
end

function main(M_precise, ϵ, symbols, valences, masses_precise)
    # Set-up
    M, masses = convert_to_ints(M_precise, masses_precise, ϵ)
    println("Set-up completed")

    # Stage 1: building up all formulae
    all_formulae = aufbau_generator(masses, valences, M)
    println("All formulae generated, L: ", length(all_formulae))

    # Stage 2: filtering 
    realizable_formulae = Vector{Vector{Int64}}[]
    for formula in all_formulae
        if !mgraph_filter(formula)
            continue
        end
        push!(realizable_formulae, formula)
    end

    # Display formulae 
    println("Only graphical formulae left, L: ", realizable_formulae)
end




main(130.14301, 1e5, element_symbols, element_valences, element_masses)