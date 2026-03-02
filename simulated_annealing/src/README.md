# Molecular Docking: Parameter Overview

This README describes the parameters defined in `parameters.py` and how they relate to the different Makefile commands used for molecular docking tasks.

## Parameters in `parameters.py`

1. **GENERATE\_OUTPUT\_DIR** (`True` or `False`)

   Whether to create the output directory for storing preprocessing data.

2. **TRAINING\_PROBABILITY** (`float`, 0.0 to 1.0)
   
   Probability that a given protein-ligand pair is included in the training set.

3. **OUTPUT\_BASE\_DIR** (`str`)
   
   Directory where output data will be saved. Must follow format `output/outputN` (e.g., `output/output1`).

4. **SEED** (`int`)
   
   Random seed for reproducibility.

5. **REMOVE\_HYDROGENS** (`True` or `False`)
   
   Whether to remove hydrogen atoms from the ligand during graph construction.

6. **LIGAND\_PATH1**, **LIGAND\_PATH2** (`str`)
   
   Paths to ligand files in `.mol2` format.

7. **POCKET\_PATH** (`str`)
   
   Path to a pocket `.pdb` file.

8. **POCKETGRID\_PDB\_PATH** (`str`)
   
   Path to a pocket grid `.pdb` file.

9. **POCKETGRID\_CSV\_PATH** (`str`)
   
   Path to a preprocessed pocket grid `.csv` file.

10. **LJ\_PARAMS\_FILE\_PATH** (`str`)
   
   Path to the file containing Lennard-Jones parameters for atoms.

11. **MAX\_HA**, **MAX\_RB** (`int`)
   
   Maximum number of heavy atoms (HA) and rotatable bonds (RB) in ligands. Valid combinations:

   * (6, 0)
   * (8, 2)
   * (10, 3)

12. **DATASET\_FILENAME** (`str`)
   
   Filename listing complexes to consider. Options:

   * `train_set.txt`
   * `test_set.txt`
   * `full_set.txt`

13. **GRID\_TYPE** (`str`)
    
    Type of grid spacing or selection. Options:

    * `''` (all points)
    * `_fix_1A`, `_fix_2A`, `_fix_3A`, `_fix_4A`

14. **USE\_AUGMENTED\_GRID** (`True` or `False`)
    
    Whether to include the true ligand positions in the grid.

15. **ENCODING** (`str`)
    
    Encoding strategy for the QUBO problem:

    * `one_hot_encoding`
    * `dw_encoding`
    * `dw_encoding_V2`

16. **LAMBDA\_WEIGHTS** (`list[float]`)
    
    Coefficients used in the Hamiltonian for energy term weighting.

17. **ALPHA\_MULTIPLIER** (`float`)
    
    Multiplier for the alpha coefficient in one-hot encoding.

18. **GAMMA\_1\_MULTIPLIER** (`float`)
    
    Multiplier for the gamma\_1 term in domain wall encoding.

19. **GAMMA\_2\_MULTIPLIER** (`float`)
    
    Multiplier for the gamma\_2 term in domain wall encoding.

20. **DIFFERENCE\_TYPE** (`str`)
    
    How to compute geometric differences:

    * `'absolute'`
    * `'squared'`

21. **VDW\_CUTOFF** (`float` or `None`)
    
    Cutoff distance for Van der Waals interactions.

22. **USE\_LJ\_12\_6** (`True` or `False`)
    
    Whether to use the 12-6 Lennard-Jones potential (otherwise 8-4).

23. **THETA\_GRANULARITY** (`int`)
    
    Number of angular points between 0 and π/2 for hydrogen bond analysis.

24. **USE\_HYDROPHOBIC\_FORM\_V2** (`True` or `False`)
    
    Whether to use the updated hydrophobic interaction form (V2).

25. **NUM\_READS** (`int`)
    
    Number of solutions to sample in simulated annealing.

26. **NUM\_SWEEPS** (`int`)
    
    Number of sweeps in simulated annealing.

27. **NUM\_EMBEDDINGS** (`int`)
    
    Number of embeddings for quantum annealing.

28. **NUM\_READS\_PER\_EMBEDDING** (`int`)
    
    Number of reads per embedding in quantum annealing.

29. **ANNEALING\_TIME** (`int`)
    
    Annealing duration (in μs) for each quantum solution.

30. **DWAVE\_TOKEN** (`str`)
    
    D-Wave API token.

31. **DWAVE\_SOLVER** (`str`)
    
    D-Wave solver name.

## Makefile Command Parameter Dependencies

Below is a mapping of Makefile commands to the parameters you may want to customize for each:

### `make preprocess`
* `LIGAND_PATH1`
* `POCKET_PATH`
* `POCKETGRID_PDB_PATH`
* `LJ_PARAMS_FILE_PATH`
* `VDW_CUTOFF`
* `USE_LJ_12_6`
* `USE_HYDROPHOBIC_FORM_V2`

### `make preprocess-dataset`

* `GENERATE_OUTPUT_DIR`
* `TRAINING_PROBABILITY`
* `OUTPUT_BASE_DIR`
* `SEED`
* `LJ_PARAMS_FILE_PATH`
* `VDW_CUTOFF`
* `USE_LJ_12_6`
* `USE_HYDROPHOBIC_FORM_V2`

### `make generate-embeddings`

* `OUTPUT_BASE_DIR`
* `REMOVE_HYDROGENS`
* `MAX_HA`
* `MAX_RB`
* `DATASET_FILENAME`
* `GRID_TYPE`
* `USE_AUGMENTED_GRID`
* `ENCODING`
* `LAMBDA_WEIGHTS`
* `ALPHA_MULTIPLIER`
* `GAMMA_1_MULTIPLIER`
* `GAMMA_2_MULTIPLIER`
* `NUM_EMBEDDINGS`
* `DWAVE_TOKEN`
* `DWAVE_SOLVER`

### `make run-sa`

* `SEED`
* `REMOVE_HYDROGENS`
* `LIGAND_PATH1`
* `POCKETGRID_CSV_PATH`
* `ENCODING`
* `LAMBDA_WEIGHTS`
* `ALPHA_MULTIPLIER`
* `GAMMA_1_MULTIPLIER`
* `GAMMA_2_MULTIPLIER`
* `DIFFERENCE_TYPE`
* `USE_LJ_12_6`
* `THETA_GRANULARITY`
* `USE_HYDROPHOBIC_FORM_V2`
* `NUM_READS`
* `NUM_SWEEPS`

### `make run-batch-sa`

* `OUTPUT_BASE_DIR`
* `SEED`
* `REMOVE_HYDROGENS`
* `MAX_HA`
* `MAX_RB`
* `DATASET_FILENAME`
* `GRID_TYPE`
* `USE_AUGMENTED_GRID`
* `ENCODING`
* `LAMBDA_WEIGHTS`
* `ALPHA_MULTIPLIER`
* `GAMMA_1_MULTIPLIER`
* `GAMMA_2_MULTIPLIER`
* `DIFFERENCE_TYPE`
* `USE_LJ_12_6`
* `THETA_GRANULARITY`
* `USE_HYDROPHOBIC_FORM_V2`
* `NUM_READS`
* `NUM_SWEEPS`

### `make run-multi-ligand-sa`

* `SEED`
* `REMOVE_HYDROGENS`
* `LIGAND_PATH1`
* `LIGAND_PATH2`
* `POCKETGRID_CSV_PATH`
* `ENCODING`
* `LAMBDA_WEIGHTS`
* `ALPHA_MULTIPLIER`
* `GAMMA_1_MULTIPLIER`
* `GAMMA_2_MULTIPLIER`
* `DIFFERENCE_TYPE`
* `USE_LJ_12_6`
* `THETA_GRANULARITY`
* `USE_HYDROPHOBIC_FORM_V2`
* `NUM_READS`
* `NUM_SWEEPS`

### `make run-batch-qa`

* `OUTPUT_BASE_DIR`
* `REMOVE_HYDROGENS`
* `MAX_HA`
* `MAX_RB`
* `DATASET_FILENAME`
* `GRID_TYPE`
* `USE_AUGMENTED_GRID`
* `ENCODING`
* `LAMBDA_WEIGHTS`
* `ALPHA_MULTIPLIER`
* `GAMMA_1_MULTIPLIER`
* `GAMMA_2_MULTIPLIER`
* `DIFFERENCE_TYPE`
* `USE_LJ_12_6`
* `THETA_GRANULARITY`
* `USE_HYDROPHOBIC_FORM_V2`
* `NUM_EMBEDDINGS`
* `NUM_READS_PER_EMBEDDING`
* `ANNEALING_TIME`
* `DWAVE_TOKEN`
* `DWAVE_SOLVER`

### `make run-qaoa`

* `REMOVE_HYDROGENS`
* `LIGAND_PATH1`
* `POCKETGRID_CSV_PATH`
* `ENCODING`
* `LAMBDA_WEIGHTS`
* `ALPHA_MULTIPLIER`
* `GAMMA_1_MULTIPLIER`
* `GAMMA_2_MULTIPLIER`
* `DIFFERENCE_TYPE`
* `USE_LJ_12_6`
* `THETA_GRANULARITY`
* `USE_HYDROPHOBIC_FORM_V2`

---

> 💡 **Tip:** Always make sure paths exist and format constraints (e.g., `output/outputN`) are respected to avoid runtime errors.

---