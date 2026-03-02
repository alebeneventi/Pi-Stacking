# Molecular Docking Project

This repository contains code and input data for studying **molecular docking** using classical simulated annealing, with a focus on **π-stacking interactions**. The project builds upon a quantum-based docking framework originally developed by collaborators as part of the **PNRR research project MOSEGAD**, which explored quantum computing in drug discovery.  

In this repository, we adopt a **classical simulated annealing framework**, further enhanced by incorporating an explicit treatment of **π-stacking interactions**.

## Repository Structure

The repository is organized into two main folders:

### 1. `simulated_annealing`
Contains Python scripts to perform molecular docking using **simulated annealing**. Key features include:

- Representation of ligand and protein pockets as discrete graphs.  
- Inclusion of **π-stacking interactions** for aromatic residues.  
- Preprocessing scripts to convert protein and ligand structures into suitable input formats.  
- Example runs and utilities to analyze docking results.

> See the `simulated_annealing/README.md` for detailed instructions on how to run the calculations.

### 2. `qchem_inputs`
Contains Q-Chem (v6.3) input files and automation scripts for calculating $\pi-\pi$ interaction energies between benzene dimers using SAPT/XSAPT protocols.

Based on the computational framework from JACS 2024, the suite includes a Python utility for rigid roto-translations of molecular geometries and a SLURM template optimized for the Leonardo (CINECA) HPC cluster. 
It is designed to streamline energy decomposition analysis and distance/rotation scans for molecular dimers.

> See the `qchem_inputs/README.md` for detailed instructions on how to run the simulations.

## Credits
- Original framework developed by collaborators as part of **MOSEGAD**.  
- Current modifications by **Alessandro Beneventi**, adding π-stacking and focusing on classical simulated annealing.  

## License
This project is licensed under the MIT License – see the LICENSE file for details.



