include("../src/elements.jl")
include("../src/generation/aufbau.jl")
include("../src/generation/naive.jl")
include("../src/filtering/mgraph.jl")
include("../src/auxiliary.jl")
include("../src/filtering/basic_organic.jl")
include("./old_m_versions.jl")
using BenchmarkTools

function arts!(realizables, seqs)
    for seq in seqs
        if !mgraph_filter(seq)
            continue
        end
        push!(realizables, seq)
    end
end

function awmi!(realizables, seqs)
    for seq in seqs
        if !ima(seq)
            continue
        end
        push!(realizables, seq)
    end
end

function awik!(realizables, seqs)
    for seq in seqs
        if !check_sequence_iter(seq)
            continue
        end
        push!(realizables, seq)
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
    compomers = enumerate_MF(masses, M)
    # @btime enumerate_MF($masses, $M)
    println("Compomers generated, L: ", length(compomers))
    if length(compomers) < 5
        println("Compomers: ", compomers)
    end
    println("")
    
    seqs = collect((build_repeat_seq(compomer, valences) for compomer in compomers))
    # dict = Dict(zip(seqs, compomers))
    # Stage 2: filtering 
    println("ARTS starts")
    realizables_arts = Vector{Vector{Int}}[]
    @time arts!(realizables_arts, seqs)
    # println("ARTS completed, timing starts...")
    # @btime arts!(Vector{Vector{Int}}[], $seqs)
    # println("Timing ARTS completed")

    # realizables_awmi = Vector{Vector{Int}}[]
    # @time awmi!(realizables_awmi, seqs)
    # println("AWMI completed, timing starts...")
    # awmi!(Vector{Vector{Int}}[], seqs)
    # println("Timing AWMI completed")

    # realizables_awik = Vector{Vector{Int}}[]
    # awik!(realizables_awik, seqs)
    # println("AWIK completed, timing starts...")
    # @time awik!(Vector{Vector{Int}}[], seqs)
    # println("Timing AWIK completed")

    println("L_ARTS: ", length(realizables_arts))
    # println("L_AWMI: ", length(realizables_awmi))
    # println("L_AWIK: ", length(realizables_awik))

    # set_awik = Set(realizables_awik)
    # set_arts = Set(realizables_arts)
    # set_awmi = Set(realizables_awmi)

    # uset_awik = setdiff(set_awik, set_arts)
    # uset_arts = setdiff(set_arts, set_awmi)
    # uset_awmi = setdiff(set_awmi, set_arts)

    # println("Found only by ARTS: ", uset_arts, " L: ", length(uset_arts))
    # println("Found only by AWMI: ", uset_awmi, " L: ", length(uset_awmi))
    # println("Found only by AWIK: ", uset_awik, " L: ", length(uset_awik))
end
# 103.12100 ; 46 ; 775
main(1030, 1, element_symbols, element_valences, element_masses)