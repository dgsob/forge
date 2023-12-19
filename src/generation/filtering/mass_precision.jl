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