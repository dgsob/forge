#= 
    This filter simply removes sequences that lack C or H or both 
=#
"""
    Returns true on compomers containing C and H. 
"""
function basic_organic_filter(compomer)
    if compomer[1] === 0 || compomer[2] === 0
        return false
    end
    return true
end