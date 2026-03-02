# Q-Chem Input Files for π--π Interaction Energy Calculations

## Overview

This repository contains the **Q-Chem input files** used to compute the interaction energy between two benzene molecules within the framework described in the SAPT/XSAPT section of the official Q-Chem Manual (v6.3, Section 12.10.4):

https://manual.q-chem.com/6.3/Ch12.S10.SS4.html

The calculations follow the computational protocol adopted in the associated scientific work and allow for the extraction of the total interaction energy and its physical components.

In addition, a Python script (`electrostatic_vdw.py`) is provided to compute, for each molecular configuration, the **classical electrostatic and van der Waals interaction energy** between the two fragments.  
This classical contribution is subtracted from the total interaction energy obtained from SAPT/XSAPT calculations in order to obtain the **residual interaction energy**, which is subsequently used in the simulated annealing framework.

------------------------------------------------------------------------

## Folder Structure

The folder includes the following Q-Chem input files:

- `dhf.in` --- Density-Fitted Hartree--Fock calculation  
- `sapt.in` --- SAPT(DFT) interaction energy decomposition  
- `xsapt.in` --- XSAPT interaction energy calculation  

Additionally:

- `molecular_geometry.in` --- Contains the Cartesian coordinates of the molecular dimer  
- `rototraslation.py` --- Python script to generate rigid roto-translations of the molecular dimer in Cartesian space  
- `geometry.in` --- Automatically generated geometry file that can be directly used as the molecular block for SAPT/XSAPT calculations  
- `electrostatic_vdw.py` --- Python script to compute classical Coulomb and Lennard-Jones (vdW) interaction energies for a given configuration  
- `job.slurm` --- Example SLURM job script for running the calculations on the Leonardo HPC system (CINECA)

------------------------------------------------------------------------

## Classical Energy Correction (`electrostatic_vdw.py`)

For each intermolecular configuration:

1. The total interaction energy is computed via SAPT or XSAPT.
2. The classical electrostatic and van der Waals contributions are computed using `electrostatic_vdw.py`.
3. The classical contribution is subtracted from the SAPT/XSAPT total energy.
4. The resulting **residual interaction energy** is used as input for the simulated annealing docking model.

------------------------------------------------------------------------

## Important Usage Instructions

Inside each input file (`dhf.in`, `sapt.in`, `xsapt.in`), the placeholder block:

<molecular_geometry>

**must be replaced** with the full content of a geometry file.

You may either:

- Use the original `molecular_geometry.in`, or  
- Use the automatically generated `geometry.in` produced by `rototraslation.py`.

This substitution is required for the calculations to run correctly.

------------------------------------------------------------------------

## Modifying the Molecular Configuration

To obtain energy profiles corresponding to different intermolecular arrangements (e.g., parallel-displaced, T-shaped, rotated configurations, or distance scans), it is possible to:

- Apply rigid roto-translations to the two benzene molecules using `rototraslation.py`  
- Generate a new `geometry.in` file  
- Replace the `<molecular_geometry>` block in the input files accordingly  
- Compute the classical correction using `electrostatic_vdw.py`  

No further modifications to the input files are required, provided the molecular geometry block is correctly replaced.

------------------------------------------------------------------------

## HPC Execution (Leonardo – CINECA)

An example SLURM job script (`job.slurm`) is provided for running the calculations on the **Leonardo HPC system at CINECA**.

The script includes:

- Proper module loading for Q-Chem  
- Scratch directory configuration  
- OpenMP parallelization setup  
- Sequential execution of DHF, SAPT, and XSAPT calculations  
- Automatic extraction of the total interaction energy  

Users running on different HPC systems may need to adapt:

- Partition and QoS settings  
- Account information  
- Module names  
- Number of cores  

------------------------------------------------------------------------

## Scientific Reference

The input files were **kindly provided by the authors** of the following article:

*Substituent and Heteroatom Effects on π--π Interactions: Evidence That Parallel-Displaced π-Stacking is Not Driven by Quadrupolar Electrostatics*,  
Journal of the American Chemical Society.  
https://pubs.acs.org/doi/10.1021/jacs.4c13291

If you use these inputs in your research, please cite the original article accordingly.

------------------------------------------------------------------------

## Software Requirements

- Q-Chem (version 6.3 or compatible)  
- Python 3.x (for running `rototraslation.py` and `electrostatic_vdw.py`)  
- SLURM workload manager (for HPC execution)  
- Sufficient computational resources for SAPT/XSAPT calculations  
