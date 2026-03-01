#!/bin/bash

#SBATCH -p boost_usr_prod
# #SBATCH -q boost_qos_dbg
# #SBATCH --time=20:00:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=1
#SBATCH -e qchem_scan.err
#SBATCH -o qchem_scan.out
#SBATCH -J qchem-benzene-scan
#SBATCH -A PHD_benevent

ml purge
ml profile/chem-phys
ml q-chem/6.1.1
source $QC/qcenv.sh

export QCSCRATCH=$CINECA_SCRATCH/qchem_scratch_${SLURM_JOB_ID}
export QCLOCALSCR=$TMPDIR
export nproc=$SLURM_NTASKS

# ─────────────────────────────────────────────────────────────────
# TRANSLATION AXES SELECTION
# Set to true/false to enable or disable translation along each axis.
# ─────────────────────────────────────────────────────────────────
TRANSLATE_X=true
TRANSLATE_Y=true
TRANSLATE_Z=false

# ─────────────────────────────────────────────────────────────────
# SCAN PARAMETERS FOR EACH AXIS (start  step  end)
# These are used only if the corresponding axis is enabled above.
# All values in Angstrom.
# ─────────────────────────────────────────────────────────────────
X_START=0.0;  X_STEP=0.2;  X_END=1.6
Y_START=0.0;  Y_STEP=0.5;  Y_END=2.0
Z_START=0.0;  Z_STEP=0.3;  Z_END=1.2

# ─────────────────────────────────────────────────────────────────
# Build the displacement sequence for each active axis.
# Inactive axes are replaced by a single "0.0" step so the loop
# still works uniformly — those coordinates will simply not change.
# ─────────────────────────────────────────────────────────────────
if $TRANSLATE_X; then
    X_VALUES=$(seq $X_START $X_STEP $X_END)
else
    X_VALUES="0.0"
fi

if $TRANSLATE_Y; then
    Y_VALUES=$(seq $Y_START $Y_STEP $Y_END)
else
    Y_VALUES="0.0"
fi

if $TRANSLATE_Z; then
    Z_VALUES=$(seq $Z_START $Z_STEP $Z_END)
else
    Z_VALUES="0.0"
fi

# Convert each sequence to an indexed array so we can iterate in parallel
X_ARR=($X_VALUES)
Y_ARR=($Y_VALUES)
Z_ARR=($Z_VALUES)

# The number of steps is driven by the longest active axis sequence.
# Shorter sequences are clamped to their last value once exhausted.
len_x=${#X_ARR[@]}
len_y=${#Y_ARR[@]}
len_z=${#Z_ARR[@]}

# Compute the maximum number of steps across all active axes
nsteps=$len_x
[ $len_y -gt $nsteps ] && nsteps=$len_y
[ $len_z -gt $nsteps ] && nsteps=$len_z

# Hartree -> kcal/mol conversion factor
conversion=627.509

# Build a descriptive label for the CSV header
axis_label=""
$TRANSLATE_X && axis_label="${axis_label}x"
$TRANSLATE_Y && axis_label="${axis_label}y"
$TRANSLATE_Z && axis_label="${axis_label}z"

# Output CSV file
echo "Step, dx(A), dy(A), dz(A), Interaction Energy (kcal/mol)" > energy_results.csv

echo "Active axes: [${axis_label}]  —  total steps: $nsteps"

# ─────────────────────────────────────────────────────────────────
# MAIN SCAN LOOP
# At each step i, retrieve dx/dy/dz from the respective arrays.
# If an array is shorter than the current index, reuse its last value.
# ─────────────────────────────────────────────────────────────────
for (( i=0; i<nsteps; i++ )); do

    # Pick the displacement for each axis (clamp to last element if exhausted)
    idx_x=$(( i < len_x ? i : len_x-1 ))
    idx_y=$(( i < len_y ? i : len_y-1 ))
    idx_z=$(( i < len_z ? i : len_z-1 ))

    dx=${X_ARR[$idx_x]}
    dy=${Y_ARR[$idx_y]}
    dz=${Z_ARR[$idx_z]}

    # If an axis is inactive, force its displacement to zero
    $TRANSLATE_X || dx=0.0
    $TRANSLATE_Y || dy=0.0
    $TRANSLATE_Z || dz=0.0

    echo "== STEP $i  |  dx=$dx  dy=$dy  dz=$dz =="

    # Helper function: apply the per-axis displacements to a base coordinate triple.
    # Arguments: base_x  base_y  base_z
    displaced_coords() {
        local bx=$1 by=$2 bz=$3
        local nx ny nz
        nx=$(echo "$bx + $dx" | bc)
        ny=$(echo "$by + $dy" | bc)
        nz=$(echo "$bz + $dz" | bc)
        echo "$nx  $ny  $nz"
    }

    # Compute displaced coordinates for each atom of molecule 1
    # (molecule 2 is kept fixed at its reference position)
    C1=$(displaced_coords  2.198  -0.286  3.500)
    C2=$(displaced_coords  0.852  -0.646  3.500)
    C3=$(displaced_coords  2.559   1.061  3.500)
    C4=$(displaced_coords  0.227   1.686  3.500)
    C5=$(displaced_coords  1.573   2.046  3.500)
    C6=$(displaced_coords -0.134   0.339  3.500)
    H1=$(displaced_coords  2.964  -1.051  3.500)
    H2=$(displaced_coords  0.571  -1.693  3.500)
    H3=$(displaced_coords -1.180   0.059  3.500)
    H4=$(displaced_coords -0.539   2.451  3.500)
    H5=$(displaced_coords  1.854   3.093  3.500)
    H6=$(displaced_coords  3.605   1.341  3.500)

    # Use a unique label for all files generated at this step
    label="${i}_dx${dx}_dy${dy}_dz${dz}"

    # ── Generate the geometry input block ───────────────────────
    cat > geometry_${label}.in << EOF
\$molecule
0 1
--
0 1
C  $C1
C  $C2
C  $C3
C  $C4
C  $C5
C  $C6
H  $H1
H  $H2
H  $H3
H  $H4
H  $H5
H  $H6
--
0 1
C  3.298  -0.286  0.000
C  1.952  -0.646  0.000
C  3.659   1.061  0.000
C  1.327   1.686  0.000
C  2.673   2.046  0.000
C  0.966   0.339  0.000
H  4.064  -1.051  0.000
H  1.671  -1.693  0.000
H -0.080   0.059  0.000
H  0.561   2.451  0.000
H  2.954   3.093  0.000
H  4.705   1.341  0.000
\$end
EOF

    # ── DHF input file ───────────────────────────────────────────
    cat > dhf_${label}.in << EOF
\$rem
jobtype xsapt
method HF
basis 6-311G
scf_convergence 7
thresh 12
\$end

\$xpol
embed none
print 3
\$end

\$sapt
cphf
dscf
print 3
algorithm ao
order 2
basis dimer
Dispersion MBD
\$end

$(cat geometry_${label}.in)
EOF

    # ── SAPT input file ──────────────────────────────────────────
    cat > sapt_${label}.in << EOF
\$rem
jobtype xsapt
method lrc-wpbe
basis def2-tzvpd
lrc_dft true
scf_guess autosad
thresh 12
\$end

\$xpol
embed none
dft-lrc
print 3
\$end

\$sapt
algorithm ao
basis projected
print 3
Dispersion MBD
\$end

\$lrc_omega
340
340
\$end

$(cat geometry_${label}.in)
EOF

    # ── XSAPT input file ─────────────────────────────────────────
    cat > xsapt_${label}.in << EOF
\$rem
jobtype xsapt
method lrc-wpbe
basis def2-tzvpd
lrc_dft true
scf_guess autosad
thresh 12
\$end

\$xpol
embed charges
dft-lrc
print 3
\$end

\$sapt
print 3
algorithm ao
order 2
basis projected
Dispersion MBD
\$end

\$lrc_omega
340
340
\$end

$(cat geometry_${label}.in)
EOF

    # ── Run the three Q-Chem calculations ────────────────────────
    qchem -np $nproc dhf_${label}.in   dhf_${label}.out
    qchem -np $nproc sapt_${label}.in  sapt_${label}.out
    qchem -np $nproc xsapt_${label}.in xsapt_${label}.out

    # Rename outputs to match the filenames expected by the Python parser
    cp dhf_${label}.out   dhf.in.out
    cp sapt_${label}.out  sapt.in.out
    cp xsapt_${label}.out xsapt.in.out

    # Run the Python energy-extraction script and append results to CSV
    energy=$(python3 extract_energy.py | tail -n 1 | awk -F',' '{print $5}')
    echo "$i, $dx, $dy, $dz, $energy" >> energy_results.csv

done

# Clean up the scratch directory
rm -rf $QCSCRATCH