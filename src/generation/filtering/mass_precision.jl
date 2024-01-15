#=
    The main generation aufbau algorithm works on integer values as masses.
    The results are filtered here and only those that lie within the bounds of the preceise mass are accepted
    as the results. It is faster this way than scaling preceise mass to int on the main function input
    becuase aufbau algorithm runs substancially slower on large mass inputs. And we don't deal with large masses
    anyway.  
=#

function mass_precision_filter(compomer, masses_precise, M_ref, ϵ)
    # if compomer == [36,38,0,19,0,0]
    #     println(sum(compomer .* masses_precise))
    # end
    M_test = sum(compomer .* masses_precise)
    if M_test ≤ M_ref + ϵ && M_test ≥ M_ref - ϵ
        return true
    end
    return false
end