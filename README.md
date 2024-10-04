## What is forge?
`forge` (FORmulae GEnerator) is a cli program whose objective is to find molecular formula of an unknown compound based on its mass and accuracy of the measurement. The ouptut is going to be a list of all formulae that are graphically realisable, accounting for all possible valences of each element (or rather the ones defined in src/elements.jl). It was built on top of the topic of my bachelor thesis as a platform for testing and experimentation. 

## Graph realization problem in this context
By graphically realisable it is meant that each listed formula has an existing structural graph. To find the actual structures one can use [`surge`](https://github.com/StructureGenerator/surge) which is a open source chemical graph generator, it takes strings as an input.  

To deal with multiple valences it uses a custom [`algorithm`](https://github.com/dgsob/forge/blob/main/src/generation/filtering/mgraph.jl) which can determine graphical realizability of a sequence of sets of possible valences, as opposed to classical "SENIOR check" which requires each element to take some value as its basic valence. "SENIOR check" is commonly utilized by chemical graph generators. In a result they may overlook rare compounds.

According to my knowledge, at the time this repository was published there was no efficient approach to solve molecular graph realization problem, where given a finite sequence of sets of integers it is determined whether there exists a molecular graph such that each degree in its degree sequence is in respective set of integers from the input sequence.

Coming up with and testing different approaches to solve this problem was presented in my bachelor thesis "Analysis of structures of chemical compounds by means of graph theory" in January 2024.

## Usage of other work
This program implements generation algorithm described here:
[`Li, S., Bohman, B., & Jayatilaka, D. (2022). 
Enumerating Possible Molecular Formulae in Mass Spectrometry Using a Generating Function Based Method. 
MATCH Communications in Mathematical and in Computer Chemistry, 88(2).`](https://match.pmf.kg.ac.rs/electronic_versions/Match88/n2/match88n2_321-350.pdf)

## Example usage
With Julia installed type in the terminal:

`julia forge.jl mass-value accuracy-value`

e.g. `julia forge.jl 88.15400 1e-1`

The final ouptut in this speciifc case is: 

`["C3H5FSi", "C3H8N2O", "C3H8OSi", "C5H9F", "C5H12O", "CH3B2ClO", "CH4OSi2", "CH4N2OSi", "CHOPSi", "CHFSi2", "CHFN2Si", "C2H5PSi", "C2H8Si2", "C2H8N2Si", "C2H8N4", "C2H4O2Si", "C2H4SSi", "C2HFOSi", "C4H12N2", "C4H12Si"]`


## Limitation
Technically increasing the accuracy should result in reduction of the number of formulas in the output, helping to identify more likely ones. Unfortunately the problem of isotope masses is not addressed. By default the masses of elements are defined as the masses of their most abundant isotopes. In the example above the provided mass is of `C4H12N2`. 

In the example output list of possible formulas, we can see that `C4H12N2` is indeed present at the n-1 position, where n is the size of the list. However, if we increase the accuracy, we will filter out this "correct" formula. This is caused by the most abundant isotopes assumption and is dependent on individual case, which is a major drawback for practical use.

## Remark
Having a mass of an unknow compound, we aim to pinpoint as little number of compounds it can correspond to. It may seem counterintuitive that we essentially increase the number of these compunds using my algorithm instead of SENIOR check. 

However, as mentioned above, the aim is to include rare compounds overlooked otherwise. Then, under the assumption that the masses of individual elements are known (which is not that simple, as pointed out in "Limitations"), we obtain a list of compound formulas that for sure contains the formula of the compound of interest. 

Theoritically speaking, such a list will already be relatively "small", as solving the molecular graph realizaton problem guarantees that each of the output formulas can be represented as a [`molecular graph`](https://en.wikipedia.org/wiki/Molecular_graph), which alone greatly reduces the number of possible atom combinations. It does not account for a more sophisticated chemical and/or biological context of molecule formation though. Thus, a number of additional filters should be used to further reduce the possible options. 

For example, see: 
[`Kind T, Fiehn O. Seven Golden Rules for heuristic filtering of molecular formulas obtained by accurate mass spectrometry. BMC Bioinformatics. 2007 Mar 27;8:105. doi: 10.1186/1471-2105-8-105. PMID: 17389044; PMCID: PMC1851972.`](https://pubmed.ncbi.nlm.nih.gov/17389044/)
