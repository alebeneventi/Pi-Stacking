# Configuration Parameters

# General settings
SEED                    = 25    # Random seed for reproducibility
REMOVE_HYDROGENS        = False  # Whether to remove hydrogens from the input structures
ADD_BENZENE_NODES       = True   # Whether to add benzene nodes to the graph representation


# Molecular docking file paths
LIGAND_PATH1            = 'data/ligands/benzene.mol2' 
POCKET_PATH             = 'data/pockets/benzene.pdb'
POCKETGRID_PDB_PATH     = 'data/pocketgrids/pocketgrid_A.pdb'
POCKETGRID_CSV_PATH     = 'data/pocketgrids/pocketgrid_A.csv'
LJ_PARAMS_FILE_PATH     = 'data/pockets/vdw_params.csv' 



# QUBO problem settings
ENCODING                = 'one_hot_encoding'           
LAMBDA_WEIGHTS          =  [1.0, 1.0, 1.0, 1.0]   # Coefficients for the Hamiltonian
ALPHA_MULTIPLIER        = 1.0        # Multiplier for the alpha coefficient in the one-hot encoding
DIFFERENCE_TYPE         = 'squared'  # Type of difference to use in the geometric Hamiltonian. Options: 'absolute' or 'squared'
VDW_CUTOFF              = None        # Van der Waals cutoff distance, None for no cutoff
USE_LJ_12_6             = True        # Whether to use the 12-6 Lennard-Jones potential or the 8-4 potential

# Simulated annealing settings
NUM_READS               = 100        # Number of reads for simulated annealing
NUM_SWEEPS              = 2000        # Number of sweeps for simulated annealing
