import numpy as np

# ==============================
# Molecular coordinates
# ==============================


# First molecule (z = 3.5 Å)

mol1 = [
    ['C', 2.198000, -0.286000, 3.500000],
    ['C', 0.852000, -0.646000, 3.500000],
    ['C', 2.559000,  1.061000, 3.500000],
    ['C', 0.227000,  1.686000, 3.500000],
    ['C', 1.573000,  2.046000, 3.500000],
    ['C', -0.134000, 0.339000, 3.500000],
    ['H', 2.964000, -1.051000, 3.500000],
    ['H', 0.571000, -1.693000, 3.500000],
    ['H', -1.180000, 0.059000, 3.500000],
    ['H', -0.539000, 2.451000, 3.500000],
    ['H', 1.854000, 3.093000, 3.500000],
    ['H', 3.605000, 1.341000, 3.500000],
]

# Second molecule (z = 0.0 Å)

mol2 = [
    ['C', 3.298000, -0.286000, 0.000000],
    ['C', 1.952000, -0.646000, 0.000000],
    ['C', 3.659000,  1.061000, 0.000000],
    ['C', 1.327000,  1.686000, 0.000000],
    ['C', 2.673000,  2.046000, 0.000000],
    ['C', 0.966000,  0.339000, 0.000000],
    ['H', 4.064000, -1.051000, 0.000000],
    ['H', 1.671000, -1.693000, 0.000000],
    ['H', -0.080000, 0.059000, 0.000000],
    ['H', 0.561000, 2.451000, 0.000000],
    ['H', 2.954000, 3.093000, 0.000000],
    ['H', 4.705000, 1.341000, 0.000000],
]

# ==============================
# Physical constants
# ==============================

k = 8.987552e9                  # Coulomb constant (N·m²/C²)
e_charge = 1.602176634e-19     # Elementary charge (C)

# Partial charges (MMFF94)
q_C = -0.15 * e_charge
q_H =  0.15 * e_charge

# Lennard-Jones parameters
eps_C = 0.06779699304291385
eps_H = 0.0216

angstrom_to_meter = 1e-10

r_C = 4.193078986609192 * angstrom_to_meter
r_H = 2.9698484809835 * angstrom_to_meter

# Conversion factor from Joule to chosen energy units
joule_conversion = 1.4393e20

# ==============================
# Energy calculation
# ==============================

coulomb_energy = 0.0
vdw_energy = 0.0

for atom1 in mol1:
    for atom2 in mol2:

        # Compute interatomic distance in meters
        dx = atom2[1] - atom1[1]
        dy = atom2[2] - atom1[2]
        dz = atom2[3] - atom1[3]

        distance = np.sqrt(dx**2 + dy**2 + dz**2) * angstrom_to_meter

        # Assign atomic parameters
        if atom1[0] == 'C':
            q1 = q_C
            r1 = r_C
            eps1 = eps_C
        else:
            q1 = q_H
            r1 = r_H
            eps1 = eps_H

        if atom2[0] == 'C':
            q2 = q_C
            r2 = r_C
            eps2 = eps_C
        else:
            q2 = q_H
            r2 = r_H
            eps2 = eps_H

        # Mixed Lennard-Jones parameters
        r_mix = 0.5 * (r1 + r2)
        eps_mix = np.sqrt(eps1 * eps2)

        # Coulomb contribution
        coulomb_energy += (k * q1 * q2 / distance) * joule_conversion

        # Lennard-Jones (12-6) contribution
        vdw_energy += eps_mix * ((r_mix / distance)**12 - 2 * (r_mix / distance)**6)

total_energy = coulomb_energy + vdw_energy

# ==============================
# Output
# ==============================

print("Coulomb Energy:", coulomb_energy)
print("Van der Waals Energy:", vdw_energy)
print("Total Interaction Energy:", total_energy)