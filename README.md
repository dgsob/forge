`forge` (FORmulae GEnerator) is a cli program that outputs organic compound molecular forumla for a given molecular mass.

Example usage:

With Julia installed type in the terminal:

`julia forge.jl mass-value precision-value`

e.g. `julia forge.jl 509.34811 1e-3`

The ouptut is going to be a list of all formulae that are graphically realisable, accounting for all possible valences of each element (or rather the ones defined in src/elements.jl).

Graphically realisable means each listed formula has an existing structural graph. To find the actual structures one can use [`surge`](https://github.com/StructureGenerator/surge) which is a open source chemical graph generator, it takes strings as an input. 


This program also implements generation algorithm described here:

Li, S., Bohman, B., & Jayatilaka, D. (2022). 
Enumerating Possible Molecular Formulae in Mass Spectrometry Using a Generating Function Based Method. 
MATCH Communications in Mathematical and in Computer Chemistry, 88(2).