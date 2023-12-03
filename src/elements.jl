# element_symbols = ["C", "H", "B", "Br", "Cl", "F", "I", "N", "O", "P", "S", "Si"]
# # order = ["C", "H", "B", "Br", "Cl", "F", "I", "N", "O", "P", "S", "Si"]
# element_masses = [12.011, 1.008, 10.81, 79.904, 35.45, 18.998403, 126.90447, 14.007, 15.999, 30.973761998, 32.06, 28.085]
# element_valences = [[2,4], [1], [3], [1,3,4,5], [1,3,5,7], [1], [1,3,4,5,7], [3,5], [2], [1,3,5], [2,4,6], [2,4]]

element_symbols = ["C", "H", "N", "O", "Cl", "I"]
element_masses = [12.011, 1.008, 14.007, 15.999, 35.45, 126.90447]
element_valences = [[2,4], [1], [3,5], [2], [1,3,5,7], [1,4,7]]

# element_symbols = ["C", "H", "Cl", "F", "N", "O", "P", "S"]
# element_masses = [12, 1, 35, 19, 14, 16, 31, 32]
# element_valences = [[2,4], [1], [1,3,5,7], [1], [3,5], [2], [1,4,5], [2,4,6]]

valences_dict = Dict(zip(element_symbols, element_valences))
symbols_dict = Dict(zip(element_valences, element_symbols))
