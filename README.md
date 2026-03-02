# Molecular Docking Project

This repository contains code and input data for studying **molecular docking** using classical simulated annealing, with a focus on **π-stacking interactions**. The project builds upon a quantum-based docking framework originally developed by collaborators as part of the **PNRR research project MOSEGAD**, which explored quantum computing in drug discovery.  

In this repository, the quantum annealing part has been removed, leaving a **classical simulated annealing approach**, and additional modeling of **π-stacking interactions** has been added.

## Repository Structure

The repository is organized into two main folders:

### 1. `simulated_annealing`
Contains Python scripts to perform molecular docking using **simulated annealing**. Key features include:

- Representation of ligand and protein pockets as discrete graphs.  
- Inclusion of **π-stacking interactions** for aromatic residues.  
- Preprocessing scripts to convert protein and ligand structures into suitable input formats.  
- Example runs and utilities to analyze docking results.

> See the `simulated_annealing/README.md` for detailed instructions on how to run the simulations.

### 2. `qchem_inputs`
Contains input files for **Q-Chem simulations** used to calibrate and validate interaction energies, including:

- Ligand and pocket geometries in `.pdb` or `.xyz` format.  
- Example input files for π-stacking energy calculations.  
- Scripts to automate batch calculations of interaction energies along translations or rotations.

## Credits
- Original framework developed by collaborators as part of **MOSEGAD**.  
- Current modifications by **Alessandro Beneventi**, adding π-stacking and focusing on classical simulated annealing.  



