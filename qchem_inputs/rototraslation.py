import numpy as np

# ─────────────────────────────────────────────────────────────────
# ROTO-TRANSLATION PARAMETERS
# Translation vector in Angstrom
translation = (0.0, 0.0, 0.0)   # (t_x, t_y, t_z)

# Rotation angles in degrees (applied in order: Rz @ Ry @ Rx)
angles = (0.0, 0.0, 0.0)       # (alpha=Rx, beta=Ry, gamma=Rz)
# ─────────────────────────────────────────────────────────────────

# Reference geometry of benzene molecule 1
# Format: (element, x, y, z)
molecule = [
    ("C",  2.198, -0.286,  3.500),
    ("C",  0.852, -0.646,  3.500),
    ("C",  2.559,  1.061,  3.500),
    ("C",  0.227,  1.686,  3.500),
    ("C",  1.573,  2.046,  3.500),
    ("C", -0.134,  0.339,  3.500),
    ("H",  2.964, -1.051,  3.500),
    ("H",  0.571, -1.693,  3.500),
    ("H", -1.180,  0.059,  3.500),
    ("H", -0.539,  2.451,  3.500),
    ("H",  1.854,  3.093,  3.500),
    ("H",  3.605,  1.341,  3.500),
]

# Fixed reference geometry of benzene molecule 2 (not transformed)
molecule2 = [
    ("C",  3.298, -0.286,  0.000),
    ("C",  1.952, -0.646,  0.000),
    ("C",  3.659,  1.061,  0.000),
    ("C",  1.327,  1.686,  0.000),
    ("C",  2.673,  2.046,  0.000),
    ("C",  0.966,  0.339,  0.000),
    ("H",  4.064, -1.051,  0.000),
    ("H",  1.671, -1.693,  0.000),
    ("H", -0.080,  0.059,  0.000),
    ("H",  0.561,  2.451,  0.000),
    ("H",  2.954,  3.093,  0.000),
    ("H",  4.705,  1.341,  0.000),
]


def rotation_matrix(alpha_deg, beta_deg, gamma_deg):
    """
    Build the combined rotation matrix R = Rz(gamma) @ Ry(beta) @ Rx(alpha).
    Angles are in degrees and follow the extrinsic (fixed-axis) convention:
      alpha -> rotation around X
      beta  -> rotation around Y
      gamma -> rotation around Z
    """
    alpha = np.radians(alpha_deg)
    beta  = np.radians(beta_deg)
    gamma = np.radians(gamma_deg)

    # Rotation around X axis
    Rx = np.array([
        [1,           0,            0],
        [0,  np.cos(alpha), -np.sin(alpha)],
        [0,  np.sin(alpha),  np.cos(alpha)],
    ])

    # Rotation around Y axis
    Ry = np.array([
        [ np.cos(beta), 0, np.sin(beta)],
        [0,             1,           0],
        [-np.sin(beta), 0, np.cos(beta)],
    ])

    # Rotation around Z axis
    Rz = np.array([
        [np.cos(gamma), -np.sin(gamma), 0],
        [np.sin(gamma),  np.cos(gamma), 0],
        [0,              0,             1],
    ])

    return Rz @ Ry @ Rx


def rototranslate(mol, R, t):
    """
    Apply rotation R and translation t to every atom in mol.
    Rotation is applied first around the molecular centroid,
    then the translation vector t is added.

    Parameters
    ----------
    mol : list of (element, x, y, z)
    R   : 3x3 numpy rotation matrix
    t   : (t_x, t_y, t_z) translation vector in Angstrom

    Returns
    -------
    list of (element, new_x, new_y, new_z)
    """
    # Compute centroid so the rotation pivots around the molecule's center
    coords = np.array([[x, y, z] for _, x, y, z in mol])
    centroid = coords.mean(axis=0)

    t_vec = np.array(t)
    transformed = []
    for elem, x, y, z in mol:
        # Shift to centroid, rotate, shift back, then translate
        r = np.array([x, y, z]) - centroid
        new_pos = R @ r + centroid + t_vec
        transformed.append((elem, *new_pos))

    return transformed


def write_geometry(mol1_transformed, mol2, filename="geometry.in"):
    """
    Write the Q-Chem $molecule block to filename.
    mol1_transformed : roto-translated coordinates of molecule 1
    mol2             : fixed coordinates of molecule 2
    """
    with open(filename, "w") as f:
        f.write("$molecule\n")
        f.write("0 1\n")

        # Molecule 1 — roto-translated
        f.write("--\n")
        f.write("0 1\n")
        for elem, x, y, z in mol1_transformed:
            f.write(f"{elem:<2s}  {x:10.6f}  {y:10.6f}  {z:10.6f}\n")

        # Molecule 2 — fixed reference
        f.write("--\n")
        f.write("0 1\n")
        for elem, x, y, z in mol2:
            f.write(f"{elem:<2s}  {x:10.6f}  {y:10.6f}  {z:10.6f}\n")

        f.write("$end\n")

    print(f"Geometry written to '{filename}'")


# ── Main ─────────────────────────────────────────────────────────
alpha, beta, gamma = angles
R = rotation_matrix(alpha, beta, gamma)

print(f"Rotation matrix (alpha={alpha}°, beta={beta}°, gamma={gamma}°):")
print(np.round(R, 6))
print(f"Translation vector: t_x={translation[0]}, t_y={translation[1]}, t_z={translation[2]}")

mol1_transformed = rototranslate(molecule, R, translation)
write_geometry(mol1_transformed, molecule2, filename="geometry.in")