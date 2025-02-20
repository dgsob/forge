## What is forge?
`forge` (FORmulae GEnerator) is a cli program that finds molecular formula of an unknown compound based on its total mass. The ouptut is a list of all formulae that are graphically realisable, accounting for all possible valences of each chemical element (defined in src/elements.jl). 

It was built on top of the topic of my bachelor thesis as a platform for testing and experimentation.

## The problem
Graph realization problem asks if the given finitie sequence of *integers* can be a degree sequence of a *simple graph*. It is a classical problem in graph theory.

The goal of the thesis was to solve a modified problem, i.e. whether the given sequence of *sets of integers* contains a degree sequence of a *molecular graph*.  
This new problem can be interpreted as whether it is possible for the given sequence of chemical elements to form a compound.

According to my knowledge, at the time of the repository's publication there was no efficient solution to such a "molecular graph realization" problem.

Coming up with and testing different approaches to solve it was presented in my bachelor thesis "Analysis of structures of chemical compounds by means of graph theory" in January 2024.

## The solution
Algorithm with theoretically polynomial time complexity, where the worst case can be approximated with $O(m^2n)$ with $m$ - the maximum size of a set, and $n$ - the number of sets in the sequence.
This is an improvement over $O(m^n)$ for blindly checking every combination.

In practice (when applied to chemical elements), as $m$ is both small and does not change much, the observed complexity resembles a linear one. 

The pseudo-code is available [here](https://github.com/dgsob/forge/blob/main/mgrc.pdf).  
Example implementation in Julia is a part of `forge` and can be found [here](https://github.com/dgsob/forge/blob/main/src/generation/filtering/mgraph.jl). 

## Additional remarks on forge

### Example usage
With Julia installed type in the terminal:

`julia forge.jl mass_value accuracy_value`

e.g. `julia forge.jl 88.15400 1e-1`

The final ouptut in this speciifc case is: 

`["C3H5FSi", "C3H8N2O", "C3H8OSi", "C5H9F", "C5H12O", "CH3B2ClO", "CH4OSi2", "CH4N2OSi", "CHOPSi", "CHFSi2", "CHFN2Si", "C2H5PSi", "C2H8Si2", "C2H8N2Si", "C2H8N4", "C2H4O2Si", "C2H4SSi", "C2HFOSi", "C4H12N2", "C4H12Si"]`


### Limitations
Technically increasing the accuracy should result in reduction of the number of formulas in the output, helping to identify more likely ones. Unfortunately the problem of isotope masses is not addressed. In the program's code - by default - the masses of elements are assumed to be and defined as the masses of their most abundant isotopes. 

In the example above the provided mass is of `C4H12N2`. In the output list of possible formulas, we can see that `C4H12N2` is indeed present at the n-1 position, where n is the size of the list. However, if we increase the accuracy, we will filter out this "correct" formula. This is caused by the most abundant isotope assumption and is dependent on individual case, which is a major drawback for practical use.

### Filtering
Having a mass of an unknow compound, we aim to pinpoint as little number of compounds it can correspond to. It may seem counterintuitive that we essentially increase the number of these compunds using the proposed MGRC algorithm instead of SENIOR check. 

However, as mentioned above, the aim is to include rare compounds overlooked otherwise. Then, under the assumption that the masses of individual elements are known (which is not that simple, as pointed out in "Limitations"), we obtain a list of compound formulas that for sure contains the formula of the compound of interest. 

Theoritically speaking, such a list will already be relatively "small", as solving the molecular graph realizaton problem guarantees that each of the output formulas can be represented as a [`molecular graph`](https://en.wikipedia.org/wiki/Molecular_graph), which alone greatly reduces the number of possible atom combinations. It does not account for a more sophisticated chemical and/or biological context of molecule formation though. Thus, a number of additional filters - or modern machine learning techniques - should be used to further reduce the possible options. 

For examples of heuristic filtering see: 
[`Kind T, Fiehn O. Seven Golden Rules for heuristic filtering of molecular formulas obtained by accurate mass spectrometry. BMC Bioinformatics. 2007 Mar 27;8:105. doi: 10.1186/1471-2105-8-105. PMID: 17389044; PMCID: PMC1851972.`](https://pubmed.ncbi.nlm.nih.gov/17389044/)

### Visualisation
To find the structures one can use [`surge`](https://github.com/StructureGenerator/surge) which is a open source chemical graph generator, it takes strings as an input.

## Usage of other work in this repo
`forge` implements a generation algorithm described here:
[`Li, S., Bohman, B., & Jayatilaka, D. (2022). 
Enumerating Possible Molecular Formulae in Mass Spectrometry Using a Generating Function Based Method. 
MATCH Communications in Mathematical and in Computer Chemistry, 88(2).`](https://match.pmf.kg.ac.rs/electronic_versions/Match88/n2/match88n2_321-350.pdf)
