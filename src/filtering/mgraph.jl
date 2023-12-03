#= 
    This is a custom algorithm for filtering sequences based on their molecular graph realizability.
    It is similar to SENIOR rule check, but instead accounts for all of the valences for a given element.
=#

function is_hydrogen(element)
    return length(element) === 1 && element[1] === 1
end
function is_monovalent(element)
    return length(element) === 1
end
function search_degree_diffs(degrees, smallest_diff)
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

function update_maxvals!(msum, msequence, seq, i, j, sdiff)

end

"""
    Try to satisfy handshake lemma. 
    It searches through the sequence for maxdegrees sequence, 
    of which elements sum is even, thus satisfying handshake lemma. 
    It will mutate both maxdegrees_sum and maxdegrees_sequence if successful.
"""
function find_even!(msum, msequence, seq)
    smallest_diff = 100 # just a large number to initialize
    i_tracker = 0
    j = 0
    for i in eachindex(seq)
        if is_hydrogen(seq[i])
            if smallest_diff !== 100
                msequence[i_tracker] = seq[i_tracker][end-j]
                msum -= smallest_diff
                return msum
            end
            return msum
        end
        if is_monovalent(seq[i])
            continue
        end
        j, diff = search_degree_diffs(seq[i], smallest_diff)
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

function recu_conditions_check!(mdegree_sequence, mdegree_sum, seq)
    # Step 1: Handshake Lemma
    # println("BEFORE: ", mdegree_sum)
    if mdegree_sum % 2 != 0
        mdegree_sum = find_even!(mdegree_sum, mdegree_sequence, seq)
    end
    # println("AFTER: ", mdegree_sum)
    if mdegree_sum % 2 != 0
        # println("Rejected on 1")
        return false
    end

    # Step 2: Connectivity
    if mdegree_sum < (length(mdegree_sequence) - 1) * 2
        # println("Rejected on 2")
        return false
    end

    # Step 3: Exclude Loops
    if mdegree_sum >= maximum(mdegree_sequence) * 2
        # println("Accepted")
        return true
    else
        indices = findall(x -> x == maximum(mdegree_sequence), mdegree_sequence)
        index = get_index(indices)
        if index === 0
            return false
        end
        mdegree_sum -= maximum(mdegree_sequence)
        new_sequence = deepcopy(seq)
        sequence = nothing
        pop!(new_sequence[index])
        if isempty(new_sequence[index])
            return false
        end
        old_max_replacement = new_sequence[index][end]
        mdegree_sequence[index] = old_max_replacement
        mdegree_sum += old_max_replacement
        # println("Got to 3")
        return recu_conditions_check!(mdegree_sequence, mdegree_sum, new_sequence)
    end
end

using BenchmarkTools
"""
    Returns true on graphically valid sequence of valences. 
"""
function mgraph_filter(sequence)
    sort!(sequence, by = x -> (-last(x), -length(x)))

    # maxdegrees_sequence = map(degrees -> last(degrees), sequence)
    maxdegrees_sequence = collect(last(degrees) for degrees in sequence)
    maxdegrees_sum = sum(maxdegrees_sequence)
    
    return recu_conditions_check!(maxdegrees_sequence, maxdegrees_sum, sequence)

end

