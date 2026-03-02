# Molecular Docking with π-Stacking

This repository builds upon a **quantum-based molecular docking framework** developed by our collaborators, originally designed to solve ligand pose search as a **Quadratic Unconstrained Binary Optimization (QUBO)** problem with both classical and quantum annealing.  

In this version, the code has been **modified to include π-stacking interactions** and to focus exclusively on **classical simulated annealing**, removing the quantum annealing components.

> Original work: part of the **PNRR research project MOSEGAD**, exploring the application of quantum computing in drug discovery.

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Methodology](#methodology)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Performance Notes](#performance-notes)
- [Main References](#main-references)

## Overview

Molecular docking is critical in drug discovery for predicting how small molecules (ligands) bind to a biological target (proteins). This project formulates docking as a **weighted subgraph isomorphism problem**, modeling:

* Ligands as **flexible weighted graphs**
* Protein pockets as **rigid 3D spatial grids**

These representations integrate both **geometric** and **chemical** characteristics, enabling accurate binding pose prediction.

## Key Features

* Ligand and pocket graph construction with detailed chemical and spatial attributes
* QUBO formulation capturing:
  * **Geometric constraints**
  * **Van der Waals forces**
  * **Electrostatic potential energy**
  * **Residual energy (related to pi-stacking)**
* Injective mapping constraints enforced within QUBO
* Solvable with classical **Simulated Annealing (SA)**

## Methodology

1. **Pocket Grid Construction**
   
   Protein binding pockets are discretized using the PASS algorithm into a spatially weighted grid.

2. **Ligand Graph Modeling**
   
   Ligands are converted into graphs where atoms are nodes, and edges encode structural constraints such as bond lengths, fixed bond angles, and chemically restricted (non-rotatable) bonds.

3. **Chemical Embedding**
   
   Chemical information is embedded into graph nodes using precomputed values from the *MMFF94 force field* and molecular interaction patterns derived from *ProLIF*.

4. **QUBO Formulation**
   
   The docking objective and constraints (e.g., injective mapping) are translated into a single QUBO expression.

## Project Structure

```
quantum-molecular-docking/
├── data/                  # Input data files (ligands, proteins, pocket points)
├── src/                   # Source code for docking algorithm and QUBO formulation
├── requirements.txt       # Python dependencies
├── makefile               # Makefile for environment management and commands
└── README.md              # Project documentation
```

## Installation

### Prerequisites

- [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) package manager
- Python 3.13.2 (handled automatically by conda)
- [Git LFS](https://git-lfs.com/) for large file support

### Clone the Repository
```bash
git clone https://github.com/francescomicucci/quantum-molecular-docking.git 
cd quantum-molecular-docking
```

### Git LFS
This project uses Git Large File Storage (LFS) to manage large files. Perform the following steps to ensure LFS is set up correctly:

```bash
git lfs install
git lfs pull
```

### Environment Setup

The project uses a makefile for easy environment management. To get started:

```bash
# Complete setup (creates environment and installs dependencies)
make setup
```

This command will:
1. Create a conda environment named `qmol_env` with Python 3.13.2
2. Install all required dependencies from `requirements.txt`

### Manual Environment Setup

If you prefer to set up the environment step by step:

```bash
# Create conda environment
make create

# Install dependencies
make install
```

### Verify Environment Setup

Check that your environment was created successfully:

```bash
make env-info
```

## Usage

### Available Commands

> **Note**: Configuration parameters for each command can be modified in `src/parameters.py`. See the [src/README.md](src/README.md) for detailed parameter descriptions and customization options for each available command.

#### Data Processing
- `make preprocess` - Run the data preprocessing pipeline on a single protein-ligand pair

#### Molecular Docking - Simulated Annealing
- `make run-sa` - Run simulated annealing on a single protein-ligand pair

#### Environment Management
- `make env-info` - Display current environment information
- `make env-export` - Export environment configuration to `environment.yml`

#### Cleanup
- `make clean` - Remove Python cache files and temporary data
- `make clean-env` - Remove the conda environment
- `make clean-all` - Complete cleanup (cache files + environment)

#### Help
- `make help` - Display all available commands with descriptions

### Example Usage

Here's a typical workflow for running the proposed approach:

```bash
# Initial setup (only needed once)
make setup

# Process your data (optional if using included data)
make preprocess

# Run docking analysis
make run-      sa        # For classical simulated annealing

# Clean up when done
make clean
```


## Main References

1. Emanuele Triuzzi, Riccardo Mengoni, Francesco Micucci, Domenico Bonanni, Daniele Ottaviani, Andrea Beccari, and Gianluca Palermo. Molecular docking via weighted subgraph isomorphism on quantum annealers. arXiv, 2025.
2. Thomas A Halgren. Merck molecular force field. II. MMFF94 van der waals and electrostatic parameters for intermolecular interactions. J. Comput. Chem., 17(5-6):520–552, April 1996.
3. Cédric Bouysset and Sébastien Fiorucci. ProLIF: a library to encode molecular interactions as fingerprints. J. Cheminform., 13(1):72, September 2021.
4. PDBbind 2020 Refined Set
5. J. Chen, T. Stollenwerk, and N. Chancellor. Performance of domain-wall encoding for quantum annealing. IEEE Transactions on Quantum Engineering, 2:1–14, 2021.