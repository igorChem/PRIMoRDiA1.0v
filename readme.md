# What is PRIMoRDiA?

PRIMoRDiA ( PRIMoRDiA Macromolecular Reactivity
Descriptors Access ) is a shared memory parallel software
written in C++ for post electronic structure calculations, that
efficiently reads output files from most used quantum mechanics packages, storing molecular information and processing
it to generate several descriptors to calculate the global and
local reactivity of molecular systems. PRIMoRDiA supports
the main reactivity descriptors of the Conceptual Density
Functional Theory, the most famous and used reactivity theory, which works from response variables of the electronic
structure of the molecules, as also other electrostatics properties.

Installation instructions and tutorials are in our [wiki](https://github.com/igorChem/PRIMoRDiA1.0v/wiki)
All theoretical background and working instructions of our program are in our userguide pdf file that can be found in this repository. 


![cover](https://github.com/igorChem/PRIMoRDiA1.0v/blob/master/cover.png)

## PRIMoRDiA 1.0v - List of Features 

1. Calculates Global reactivity descriptors 
  (Hardness, softness, Ionization Potential, Electron Affinity, and more six related with electronic energies)
2. Calculate Local reactivity Descriptors
(Four working methods for Local hardness, Local Softness, hyper softness, local electrophilicity, multiphilcity, Fukui functions, electron density, and other common local electrostatic properties)
3. Caluclate Local Reactivity Descriptors for residues from biological systems ( set in PDB files, such as protein adn DNA fragments )
4. Implements local reactivity descriptor methods adjusted for macromolecules
5. Outputs those reactivity descriptors either in Cube files or in numerical values assigned to each atom.
6. Calculate total electron density.
7. Calculate Molecular Orbitals.
8. Make R scripts to automate Density of States (DOS) plot.


# Dowload and Installation

You can download the executable file or just compule the software usinf CMAKE. 
We strongly recommend the latter to maximize the performance of the software in your machine. 
The detailed instructions of these two intallation options are described in our [wiki](https://github.com/igorChem/PRIMoRDiA1.0v/wiki) home. 

# Tutorials and Userguide

The following links present four basic tutorials on PRIMoRDiA usage. Running PRIMoRDiA is a very simple task, although the visualization of the results in graphical packages can become a complex task as the volume of results produced by PRIMoRDiA increases. Thus, the tutorials are an invitation to users to discuss the meaning of these theoretical quantities and understand how to transform them into useful chemical information.

1. [Tutorial 1](https://github.com/igorChem/PRIMoRDiA1.0v/wiki/Tutorial-1:-Calculating-Frozen-Orbital-Reactivity-Descriptors): Frozen Orbital Method calculations
2. [Tutorial 2](https://github.com/igorChem/PRIMoRDiA1.0v/wiki/Tutorial-2:-Calculating-Finite-Differences-Reactivity-Descripors): Finite Differences Method calculations
3. [Tutorial 3](https://github.com/igorChem/PRIMoRDiA1.0v/wiki/Tutorial-3:-Calculating-Reactivity-Descriptors-for-Macromolecules): Band Reactivity Descriptors
4. [Tutorial 4](https://github.com/igorChem/PRIMoRDiA1.0v/wiki/Tutorial-4:-Electron-Density-and-Molecular-Orbitals-Generation): Electron Density and Molecular Orbital generation

A complete description of the descriptors and their theory can be found in our [userguide](https://github.com/igorChem/PRIMoRDiA1.0v/blob/master/userguide/userguide.pdf) pdf file.


# Cite our Work

Read our published works: 

1. Grillo, I. B., Urquiza‐Carvalho, G. A., Chaves, E. J. F., & Rocha, G. B. (2020). Semiempirical methods do Fukui functions: Unlocking a modeling framework for biosystems. Journal of Computational Chemistry, 41(9), 862-873.

2. Grillo, I. B., Urquiza-Carvalho, G. A., Bachega, J. F. R., & Rocha, G. B. (2020). Elucidating Enzymatic Catalysis Using Fast Quantum Chemical Descriptors. Journal of Chemical Information and Modeling, 60(2), 578-591.

3. Grillo, Igor Barden, et al. "Theoretical characterization of the shikimate 5-dehydrogenase reaction from Mycobacterium tuberculosis by hybrid QC/MM simulations and quantum chemical descriptors." Journal of Molecular Modeling 26.11 (2020): 1-12.

4. Grillo, Igor Barden, Gabriel A. Urquiza-Carvalho, and Gerd Bruno Rocha. "PRIMoRDiA: A Software to Explore Reactivity and Electronic Structure in Large Biomolecules." Journal of Chemical Information and Modeling 60.12 (2020): 5885-5890.



