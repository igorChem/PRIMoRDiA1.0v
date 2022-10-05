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


![cover](https://github.com/igorChem/PRIMoRDiA1.0v/blob/master/Repo_images/cover.png)

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

## PRIMoRDiA 1.1v - Unstable Internal Version

A not released version for the public. 
Features related with  TERACHEM support via molden files.
In the next versions this feature will be included along with the support for other quantum chemical packages

## PRIMoRDiA 1.2v - Important Stable Update  

The newest stable release.
Currently, this version is under tests and will be released soon and will include all features of the 1.0v along with important updates listed below. 

1. All features of the 1.0v.
2. New input format, more intuitive and complete for the new analysis.
3. Reaction path trajectory automatic calculation and generation of scripts for statistical analysis.
4. Molecular Dynamics trajectory automatic calculation and generation of scripts for statistical analysis.
5. Implementation of the global and local composite hardness.
6. Important adjustments in scripts for Pymol grpahical visualization 

## PRIMoRDiA 1.25v - Stable Update (upcoming soon)


# Dowload and Installation

You can download the executable file or just compule the software usinf CMAKE. 
We strongly recommend the latter to maximize the performance of the software in your machine. 
The detailed instructions of these two intallation options are described in our [wiki](https://github.com/igorChem/PRIMoRDiA1.0v/wiki) home. 


# Userguide 

We provide a complete description of the descriptors and how to use the program in a user guide pdf file. 
All the theory basis, references and how to prepare the quantum chemistry calculations are detailed in these files. 
[English](https://github.com/igorChem/PRIMoRDiA1.0v/blob/master/user_guide_EN.pdf)
[Portuguese](https://github.com/igorChem/PRIMoRDiA1.0v/blob/master/user_guide_EN.pd)

# Tutorials

We provide full tutorials in a PDF file of all calculations modes of the software, as pratical applications for acid force and enzymatic catalysis reaction. 
[Portuguese]()

We also provide an specific tutorial to run PRIMoRDiA without needing to configure your local machine. 
This tutoral uses the [Google Colab cloud computing plataform](https://colab.research.google.com) and thus you can run from any operational system machine.
[Tutorial Google Colabs](https://github.com/igorChem/PRIMoRDiA1.0v/blob/master/Tutorial_PRIMoRDiA_Colab.ipynb)

# Cite our Work

Read our published works: 

1. Grillo, I. B., Urquiza‐Carvalho, G. A., Chaves, E. J. F., & Rocha, G. B. (2020). Semiempirical methods do Fukui functions: Unlocking a modeling framework for biosystems. Journal of Computational Chemistry, 41(9), 862-873.

2. Grillo, I. B., Urquiza-Carvalho, G. A., Bachega, J. F. R., & Rocha, G. B. (2020). Elucidating Enzymatic Catalysis Using Fast Quantum Chemical Descriptors. Journal of Chemical Information and Modeling, 60(2), 578-591.

3. Grillo, Igor Barden, et al. "Theoretical characterization of the shikimate 5-dehydrogenase reaction from Mycobacterium tuberculosis by hybrid QC/MM simulations and quantum chemical descriptors." Journal of Molecular Modeling 26.11 (2020): 1-12.

4. Grillo, Igor Barden, Gabriel A. Urquiza-Carvalho, and Gerd Bruno Rocha. "PRIMoRDiA: A Software to Explore Reactivity and Electronic Structure in Large Biomolecules." Journal of Chemical Information and Modeling 60.12 (2020): 5885-5890.



