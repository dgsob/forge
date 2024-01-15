`forge` (FORmulae GEnerator) is a cli program whose objective is to find molecular formula of an unknown compound based on its mass and accuracy of the measurement. 

Example usage:

With Julia installed type in the terminal:

`julia forge.jl mass-value accuracy-value`

e.g. `julia forge.jl 509.34811 1e-3`

The ouptut is going to be a list of all formulae that are graphically realisable, accounting for all possible valences of each element (or rather the ones defined in src/elements.jl).

Graphically realisable means each listed formula has an existing structural graph. To find the actual structures one can use [`surge`](https://github.com/StructureGenerator/surge) which is a open source chemical graph generator, it takes strings as an input.  

To deal with multiple valences it uses a custom algorithm which can determine graphical realizability of a sequence of sets of possible valences, as opposed to classical SENIOR check which requires each element to take some value as its basic valence.

This program also implements generation algorithm described here:

Li, S., Bohman, B., & Jayatilaka, D. (2022). 
Enumerating Possible Molecular Formulae in Mass Spectrometry Using a Generating Function Based Method. 
MATCH Communications in Mathematical and in Computer Chemistry, 88(2).

Unfortunately the problem of isotope masses is not addressed. By default the masses of elements are defined as the masses of their most abundant isotopes.