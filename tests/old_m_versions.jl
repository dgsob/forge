function check_conditions(sequence)
    total_sum = sum(sequence)
    max_degree = maximum(sequence)
    num_degrees = length(sequence)
    
    if (total_sum % 2 == 0) && (total_sum >= max_degree * 2) && (total_sum >= (num_degrees - 1) * 2)
        return true
    else
        return false
    end
end

function check_sequence_iter(seq::Vector{Vector{Int64}})
    n = length(seq)
    indices = ones(Int, n)

    current_combination = Vector{Int}(undef, n)
    for i in eachindex(seq)
        current_combination[i] = seq[i][1]
    end

    while true
        if check_conditions(current_combination)
            return true
        end

        indices[end] += 1

        for i in n:-1:2
            if indices[i] > length(seq[i])
                indices[i] = 1
                indices[i - 1] += 1
            end
        end

        for i in 1:n
            if indices[i] <= length(seq[i])
                current_combination[i] = seq[i][indices[i]]
            else
                break
            end
        end

        if indices[1] > length(seq[1])
            break
        end
    end

    return false
end

function try_to_satisfy_handshake_lemma!(largest_degrees_sum, largest_sequence, sequence)
    for i in eachindex(sequence)
        if length(sequence[i]) < 2
            break
        end
        for j in eachindex(sequence[i])
            if j == lastindex(sequence[i])
                break
            end
            difference = sequence[i][end] - sequence[i][end-j] # sequence[i][end] = largest_sequence[i]
            if difference % 2 != 0 
                largest_sequence[i] = sequence[i][end-j]
                largest_degrees_sum -= difference
                return largest_degrees_sum
            end
        end
    end
    return largest_degrees_sum
end

function lookup_smaller_valences(act)
    for i in eachindex(act)
        if length(act[i]) == 1
            return false
        end
        for j in eachindex(act[i][begin:end-1])#1:length(act[i]) - 1
            if (act[i][end] - act[i][j]) == 1
                return true
            end
        end
    end
    return false
end

function count_possible_connections(met::AbstractArray{T}) where T
    num_of_the_same_edges = 0
    unique_edges_ends = Set()
    for i in eachindex(met)
        # Skip first element (for we compare with i-1) as well as those of size one (not a completed edge)
        if length(met[i]) != 2 || i == 1
            continue
        end
        # count redundant multiple edges - more obvious than the condition below
        if isequal(met[i], met[i-1]) 
            num_of_the_same_edges += 1
        # weaken the condition for graph connectivity - count unnecessary edges which are single connections
        elseif met[i][1] != 1 && met[i][1] in unique_edges_ends && met[i][2] in unique_edges_ends
            num_of_the_same_edges += 1
        end
        # update the set now  
        push!(unique_edges_ends, met[i][2])
    end

    return num_of_the_same_edges * 2
end

"""
Fires when there is too many hydrogen atoms left in the sequence.
met - multiple_edges_tracker    
"""
function if_too_many(met, num_of_open_connections, num_of_hydrogens, act)
    num_of_possible_connections = count_possible_connections(met)
    # println(num_of_possible_connections)
    num_of_leftover_hydrogens = num_of_hydrogens - num_of_open_connections
    num_of_final_connections = num_of_possible_connections - num_of_leftover_hydrogens
    if num_of_final_connections < 0
        return false
    end
    if num_of_final_connections == 0 || num_of_final_connections % 2 == 0
        return true
    else
        return lookup_smaller_valences(act)
    end
    return false
end

#= 
    TODO: This does not yet handle the case when there is actually too few hydrogens.
    In order to do so it has to look up different valences of the atoms that can make connections
    and check if they can satsify handshake lemma in case of the laargest degrees not being able to.
=#
function if_too_few(act, num_of_open_connections, num_of_hydrogens)
    num_of_open_connections = num_of_open_connections - num_of_hydrogens
    # check if open_connection can make connections with themselves
    # at this point we know the number of connections each atom is trying to make = act wthout hydrogens
    if num_of_open_connections % 2 == 0
        return true
    else
        # find an atom where the difference between largest degree and any other is odd
        # check if num_of_open_connections - this difference is larger than 0
        # if so return true 
        return lookup_smaller_valences(act)
    end
    return false
end

function build_incidence_matrix!(matrix, degree, i, met)
    if degree >=4 # if the element makes 4 or more connections set 3 as the max it can make with one other element
        bound = degree - 3
    else # if it makes fewer than 4, set it to 2 or 1
        bound = 1
    end
    # Main loop of this part where we assign edges beginings and ends to vertices by setting 1s in the matrix
    for j in eachindex(eachcol(matrix)) # eachcol(matrix) or matrix[1,:]? they seem to perform similarily 
        if degree == 0
            break
        end
        if degree > bound
            # sum calls could be replaced by custom functions built for this specific purpose and thus possibly increasing performance
            if sum(matrix[:,j]) != 2
                matrix[i,j] = 1 # we set 1 here whereever there's place 
                degree -= 1
                push!(met[j], i)
            end
        else
            if sum(matrix[:,j]) == 0 # when we reached the bound, we only 'reach out hands' for potentail future elements
                matrix[i,j] = 1
                degree -= 1
                push!(met[j], i)
            end
        end
    end
end

function ima(sequence, exclude_without_hydrogens = false)
    # ----------------------- Part 1 - Initialization -------------------------------------------
    #= The sequence should come in sorted by its elements (sets) size, with tailing hydrogens.
    Either perform sorting here or add it outside of this file =#
    sequence = sort(sequence, by = x -> (-length(x), -x[1]))
    largest_sequence = Vector{Int}(undef, length(sequence))
    largest_degrees_sum = 0
    num_of_hydrogens = 0

    for i in eachindex(sequence)
        degree = sequence[i][end]
        largest_sequence[i] = degree
        largest_degrees_sum += degree
        if degree == 1
            num_of_hydrogens += 1
        end
    end
    # not a compound but a single element which maybe happened to have the same mass as something we are looking for
    if length(largest_sequence) < 2
        return false
    end
    # Handshake Lemma
    #= 
      If we assume that valences in all elements that can have a varying number of valences 
      are either all even or odd, e.g. [3,5], [2,4] but not [2,3,4], then we can check the condition 
      for any cartesian product of the initial sequence. Based on this assumption we can just return false
      if largest_degrees_sum % 2 != 0 at this point.
    =#
    # But here is an alternative solution working around the above assumption, needs more testing though.
    if largest_degrees_sum % 2 != 0 #&& !assumption
        # ints are immutable so we reassign them, largest_sequence however is mutated by this call:
        largest_degrees_sum = try_to_satisfy_handshake_lemma!(largest_degrees_sum, largest_sequence, sequence)
    end
    # If after all the stretching it is still odd, give up and return false
    if largest_degrees_sum % 2 != 0
        return false
    end
    # Initialization Continues
    # Vertices
    rows = length(sequence) # i.e. vertices/atoms
    # Edges
    # if handshake lemma satisied, edges exist: 
    columns = Integer(largest_degrees_sum/2) # i.e. edges/connections

    matrix = zeros(Int, rows, columns)
    too_many_hydrogens = false
    too_few_hydrogens = false
    num_of_open_connections = 0
    multiple_edges_tracker = [Int[] for _ in 1:columns]
    # ----------------------- Part 2 - The Matrix and Hydrogens --------------------------------------------------
    # We build the graph only for elements that have varying valences, i.e. not hydrogens
    for i in eachindex(sequence[begin:end-num_of_hydrogens])
        degree = sequence[i][end] # get largest degree
        build_incidence_matrix!(matrix, degree, i, multiple_edges_tracker)
    end

    # for row in eachrow(matrix)
    #     println(row)
    # end
    
    # After this loop there are only hydrogens left
    matrix = nothing # no longer needed
    for j in eachindex(multiple_edges_tracker)
        if length(multiple_edges_tracker[j]) == 1
            num_of_open_connections += 1
        end
        if num_of_open_connections == num_of_hydrogens
            return true
        end
    end
    if num_of_open_connections < num_of_hydrogens
        too_many_hydrogens = true
    elseif num_of_open_connections > num_of_hydrogens
        too_few_hydrogens = true
    end

    # ----------------------- Part 3 - Finishing --------------------------------------------------------
    is_graphical = false
    # Handle the case where there are no hydrogens
    if num_of_hydrogens == 0
        if !exclude_without_hydrogens
            if num_of_open_connections % 2 == 0
                return true
            end
            return false
        else
            return false
        end
    end
    # Compound for sure contains hydrogen at this point
    if too_many_hydrogens
        is_graphical = if_too_many(multiple_edges_tracker, num_of_open_connections, num_of_hydrogens, sequence)
    elseif too_few_hydrogens
        is_graphical = if_too_few(sequence, num_of_open_connections, num_of_hydrogens)
    end
    return is_graphical
end

#= 
    This is a custom algorithm for filtering sequences based on their molecular graph realizability.
    It is similar to SENIOR rule check, but instead accounts for all of the valences for a given element.
=#

function is_monovalent(element)
    return length(element) === 1
end
function search_degree_diffsm(degrees, smallest_diff)
    degree_tracker = 0
    for i in eachindex(degrees)
        if i === lastindex(degrees)
            return degree_tracker, smallest_diff
        end
        diff = degrees[end] - degrees[end-i]
        if diff % 2 != 0 
            if diff === 1 
                return i, diff
            end
            if diff < smallest_diff
                smallest_diff = diff
                degree_tracker = i
            end
        end

    end
end

"""
    Try to satisfy handshake lemma. 
    It searches through the sequence for maxdegrees sequence, 
    of which elements sum is even, thus satisfying handshake lemma. 
    It will mutate both maxdegrees_sum and maxdegrees_sequence if successful.
"""
function find_evenm!(msum, msequence, seq)
    smallest_diff = 100 # just a large number to initialize
    i_tracker = 0
    j = 0
    for i in eachindex(seq)
        if is_monovalent(seq[i])
            if smallest_diff !== 100
                msequence[i_tracker] = seq[i_tracker][end-j]
                msum -= smallest_diff
                return msum
            end
            return msum
        end
        # if is_monovalent(seq[i])
        #     continue
        # end
        j, diff = search_degree_diffsm(seq[i], smallest_diff)
        if diff === 1
            msequence[i] = seq[end-j]
            msum -= diff
            return msum
        elseif diff !== 100
            smallest_diff = diff
            i_tracker = i
        end
    end
    return msum
end

function get_index(indices)
    if length(indices) > 1
        for i in eachindex(indices)
            if i === lastindex(indices)
                return 0
            end
            if sequence[indices[i]] != sequence[indices[i+1]]
                return indices[1] # TODO: return indices[i] instead and pop all up to indices[i]
            end
        end
    else
        return indices[1]
    end
end

function mrecu_conditions_check!(mdegree_sequence, mdegree_sum, seq)
    # Step 1: Handshake Lemma
    if mdegree_sum % 2 != 0
        mdegree_sum = find_evenm!(mdegree_sum, mdegree_sequence, seq)
    end
    if mdegree_sum % 2 != 0
        # println("Rejected")
        return false
    end

    # Step 2: Connectivity
    if mdegree_sum < (length(mdegree_sequence) - 1) * 2
        return true
    end

    # Step 3: Exclude Loops
    if mdegree_sum >= maximum(mdegree_sequence) * 2
        return true
    else
        indices = findall(x -> x == maximum(mdegree_sequence), mdegree_sequence)
        index = get_index(indices)
        if index === 0
            return false
        end
        mdegree_sum -= maximum(mdegree_sequence)
        new_sequence = deepcopy(seq)
        seq = nothing
        pop!(new_sequence[index])
        if isempty(new_sequence[index])
            return false
        end
        old_max_replacement = new_sequence[index][end]
        mdegree_sequence[index] = old_max_replacement
        mdegree_sum += old_max_replacement
        println("Got to 3")
        return mrecu_conditions_check!(mdegree_sequence, mdegree_sum, new_sequence)
    end
end

using BenchmarkTools
"""
    Returns true on graphically valid sequence of valences. 
"""
function mod_mgraph_filter(sequence)
    sort!(sequence, by = x -> (-last(x), -length(x)))

    # maxdegrees_sequence = map(degrees -> last(degrees), sequence)
    maxdegrees_sequence = collect(last(degrees) for degrees in sequence)
    maxdegrees_sum = sum(maxdegrees_sequence)
    
    return mrecu_conditions_check!(maxdegrees_sequence, maxdegrees_sum, sequence)

end

