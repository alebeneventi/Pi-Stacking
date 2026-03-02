#!/bin/bash

############################
# SLURM Job Configuration
############################

#SBATCH -p boost_usr_prod          # Partition
#SBATCH -q boost_qos_dbg           # QoS
#SBATCH -N 1                       # Number of nodes
#SBATCH -n 32                      # Total MPI tasks (Q-Chem cores)
#SBATCH --cpus-per-task=1          # OpenMP threads per task
#SBATCH -e benzene.err             # Standard error file
#SBATCH -o benzene.out             # Standard output file
#SBATCH -J benzene-total-energy    # Job name
#SBATCH -A <project_account>       # Project account

############################
# Module Loading
############################

module load profile/chem-phys
module load q-chem/6.1.1

source $QC/qcenv.sh

############################
# Environment Configuration
############################

# Define number of processors from SLURM
export QCNPROCS=$SLURM_NTASKS

# Define scratch directories
export QCSCRATCH=$CINECA_SCRATCH/qchem_scratch_${SLURM_JOB_ID}
export QCLOCALSCR=$TMPDIR

# Create scratch directory if not existing
mkdir -p $QCSCRATCH

############################
# Function to Run Q-Chem
############################

run_qchem() {
    input_file=$1
    output_file=${input_file%.in}.out

    echo "========================================="
    echo "Running Q-Chem calculation: $input_file"
    echo "Using $QCNPROCS cores"
    echo "Scratch directory: $QCSCRATCH"
    echo "========================================="

    qchem -np $QCNPROCS $input_file $output_file

    if [ $? -ne 0 ]; then
        echo "Error detected during execution of $input_file"
        exit 1
    fi
}

############################
# Run Calculations Sequentially
############################

run_qchem dhf.in
run_qchem sapt.in
run_qchem xsapt.in

############################
# Extract Final Interaction Energy
############################

echo "Extracting total interaction energy..."

python3 extract_energy.py | tail -n 1 | awk -F',' '{print "Total interaction energy (kcal/mol):", $5}'

############################
# Cleanup
############################

echo "Cleaning scratch directory..."
rm -rf $QCSCRATCH

echo "Job completed successfully."
