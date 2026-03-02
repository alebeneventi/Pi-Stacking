import os
import neal
import re
import warnings
# import graph_utils as gu
import numpy as np
from rdkit import Chem
from constants import ATOMIC_RADII

def read_molecule_list(path, filename):
    """
    Reads a list of molecule names from a file.
    
    Parameters
    ----------
    path : str
        The directory where the file is located.
    filename : str
        The name of the file containing the molecule names.
    
    Returns
    -------
    list
        A list of molecule names.
    """
    with open(os.path.join(path, filename), 'r') as file:
        folder_names = file.readlines()
    return [name.strip() for name in folder_names]

def compute_angle(A, B, C):
    """
    Computes the angle ABC between three points A, B, C (in 3D space) in degrees.
    The angle is formed at point B between vectors BA and BC.
    
    Parameters
    ----------
    A : numpy.ndarray
        Coordinates of point A.
    B : numpy.ndarray
        Coordinates of point B.
    C : numpy.ndarray
        Coordinates of point C.
    
    Returns
    -------
    float
        The angle ABC in degrees.
    """
    BA = A - B
    BC = C - B

    dot_product = np.dot(BA, BC)
    magnitude_BA = np.linalg.norm(BA)
    magnitude_BC = np.linalg.norm(BC)

    cos_theta = np.clip(dot_product / (magnitude_BA * magnitude_BC), -1.0, 1.0)
    
    return np.degrees(np.arccos(cos_theta))

def get_coordinates_from_mol_object(mol):
    """
    Get the coordinates of the atoms in a RDKit molecule object.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.
    
    Returns
    -------
    numpy.ndarray
        Array with the coordinates of the atoms.
    """
    conf = mol.GetConformer()
    if conf is None:
        raise ValueError("No conformer found for the molecule.")

    atom_coords = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atom_coords.append([pos.x, pos.y, pos.z])
    atom_coords = np.array(atom_coords)

    return atom_coords

def pistacking_energy(d, z, a, b):
    """
    Compute the pi-stacking energy considering position f(d,z) and modify the result
    based on the value of f and the angular factor g(a,b).

    If f < 0:
        return f * g
    If f >= 0:
        return f + f * (1 - g)

    Parameters
    ----------
    d : float
        Planar distance component.
    z : float
        Vertical distance component.
    a : float
        First angular parameter.
    b : float
        Second angular parameter.

    Returns
    -------
    float
        Modified pi-stacking energy.
    """
    # f(d, z)
    f = (
        -0.4126    * d**3
        +12.3218   * z**3
        -0.7296    * d**2 * z
        + 4.6375   * d * z**2
        + 5.3407   * d**2
        -150.7733  * z**2
        -33.4231   * d * z
        +55.5721   * d
        +613.0924  * z
        -827.5303
    )

    # g(a, b)
    g = (
        -0.0000089524 * a**3
        + 0.0000156667 * b**3
        + 0.0000133571 * a**2 * b
        - 0.0000025429 * a * b**2
        - 0.0006714286 * a**2
        - 0.0009153333 * b**2
        - 0.0000479286 * a * b
        + 0.0006638095 * a
        + 0.0005151190 * b
        + 1.0018190476
    )


    if f < 0:
        return f * g
    else:
        return f + f * (1 - g)



def get_bonds_from_mol_object(mol):
    """
    Get the bonds of the atoms in a RDKit molecule object.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.
    
    Returns
    -------
    list
        Array with the bonds of the atoms.
    """
    bonds = []
    for bond in mol.GetBonds():
        bonds.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
    return bonds

def get_atom_types_from_mol_object(mol):
    """
    Get the atom types of the atoms in a RDKit molecule object.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.
    
    Returns
    -------
    list
        Array with the atom types of the atoms.
    """
    atom_types = []
    for atom in mol.GetAtoms():
        atom_types.append(atom.GetSymbol())
    return atom_types

def find_substructure_matches(mol, smarts_pattern):
    """
    Finds matches of a substructure (query) within a target molecule.
    
    Parameters
    ----------
    mol : RDKit Mol
        The target molecule to search in.
    smarts_pattern : str
        The SMARTS pattern to search for.
        
    Returns
    -------
    list of tuple
        A list of tuples where each tuple contains the atom indices of the match.
    """    
    query_mol = Chem.MolFromSmarts(smarts_pattern)
    if query_mol is None:
        raise ValueError("Invalid SMARTS pattern for the query molecule.")
    
    matches = mol.GetSubstructMatches(query_mol)
    
    return matches

def solve_qubo(
    model,
    seed=None,
    num_reads=1000,
    num_sweeps=10000
):
    """
    Solve the QUBO problem using the Simulated Annealing sampler.
    
    Parameters
    ----------
    model : pyqubo.Model
        The model to solve.
    seed : int
        Seed for reproducibility.
    num_reads : int
        Number of reads (samples) to generate.
    num_sweeps : int
        Number of sweeps to perform for each read.
    
    Returns
    -------
    tuple
        A tuple containing the energy of the best sample and the best sample itself.
    """
    bqm = model.to_bqm()
    sa  = neal.SimulatedAnnealingSampler()
    sampleset = sa.sample(
        bqm, 
        seed=seed,
        num_reads=num_reads, 
        num_sweeps=num_sweeps
    )
    samples = model.decode_sampleset(sampleset)
    best_sample = min(samples, key=lambda s: s.energy)
    return best_sample.energy, best_sample.sample


def get_mapping_pairs_one_hot_encoding(dictionary):
    """
    This function accepts a dictionary where the keys are strings (e.g., 'X[i][j]') and 
    the values are either 0 or 1. It returns a list of tuples representing the pairs 
    (i, j) where the corresponding value is 1.
    
    Parameters
    ----------
    dictionary : dict
        Dictionary where the keys are strings and the values are either 0 or 1.
    
    Returns
    -------
    list[tuple]
        List of tuples representing the pairs (i, j) where the corresponding value is 1.
    """
    pairs = []
    for key, value in dictionary.items():
        if value == 1:
            i, j = re.findall(r'\d+', key)
            pairs.append((int(i), int(j)))
    return pairs


def compute_rmsd(
    ligand_graph,
    pocketgrid_graph,
    mapping
):
    """
    Compute the Root Mean Square Deviation (RMSD) between the true position of the 
    ligand in the pocket and its position derived from the solution of the QUBO problem.
    
    Parameters
    ----------
    ligand_graph : nx.Graph
        Graph representation of the ligand molecule.
    pocketgrid_graph : nx.Graph
        Graph representation of the pocketgrid.
    list[tuple]
        List of tuples representing the pairs (i, j) corresponding to the mapping
        between the atoms of the ligand and the pocketgrid.
    
    Returns
    -------
    float
        The RMSD between the true position of the ligand and its position in the 
        pocketgrid.
    """
    rmsd = 0.0
    list_ligand_nodes = list(ligand_graph.nodes())
    list_pocketgrid_nodes = list(pocketgrid_graph.nodes())
    for (i, j) in mapping:
        ligand_coords = ligand_graph.nodes[list_ligand_nodes[i]]['coords']
        pocketgrid_coords = pocketgrid_graph.nodes[list_pocketgrid_nodes[j]]['coords']
        rmsd += np.linalg.norm(ligand_coords - pocketgrid_coords)**2
    rmsd = np.sqrt(rmsd / len(mapping)) 
    return rmsd

def compute_rmsd_wrt_exact_mapping(
    ligand_graph,
    pocketgrid_graph,
    mapping
):
    """
    Compute the Root Mean Square Deviation (RMSD) between the position of the ligand
    when mapped to the closest point in the pocketgrid and its position derived from
    the solution of the QUBO problem.
    
    Parameters
    ----------
    ligand_graph : nx.Graph
        Graph representation of the ligand molecule.
    pocketgrid_graph : nx.Graph
        Graph representation of the pocketgrid.
    list[tuple]
        List of tuples representing the pairs (i, j) corresponding to the mapping
        between the atoms of the ligand and the pocketgrid.
    
    Returns
    -------
    float
        The RMSD between the position of the ligand when mapped to the closest point
        in the pocketgrid and its position derived from the solution of the QUBO problem.
    """
    rmsd = 0.0
    list_ligand_nodes = list(ligand_graph.nodes())
    list_pocketgrid_nodes = list(pocketgrid_graph.nodes())
    
    pocketgrid_points = np.array([pocketgrid_graph.nodes[v]['coords'] for v in list_pocketgrid_nodes])
    
    for (i, j) in mapping:
        ligand_coords = ligand_graph.nodes[list_ligand_nodes[i]]['coords']
        pocketgrid_coords = pocketgrid_graph.nodes[list_pocketgrid_nodes[j]]['coords']
        closest_point = pocketgrid_points[np.argmin(np.linalg.norm(pocketgrid_points - ligand_coords, axis=1))]
        # closest_point = pocketgrid_points[np.argmin(np.sum(np.abs(pocketgrid_points - ligand_coords), axis=1))]
        rmsd += np.linalg.norm(closest_point - pocketgrid_coords)**2
    rmsd = np.sqrt(rmsd / len(mapping))
    return rmsd

def compute_adjusted_rmsd(
    ligand_graph,
    pocketgrid_graph,
    mapping
):
    """
    Compute the difference between the RMSD between the true position of the ligand in
    the pocket and its position derived from the solution of the QUBO problem and the
    RMSD between the position of the ligand when mapped to the closest point in the
    pocketgrid and its position derived from the solution of the QUBO problem.
    
    Parameters
    ----------
    ligand_graph : nx.Graph
        Graph representation of the ligand molecule.
    pocketgrid_graph : nx.Graph
        Graph representation of the pocketgrid.
    list[tuple]
        List of tuples representing the pairs (i, j) corresponding to the mapping
        between the atoms of the ligand and the pocketgrid.
    
    Returns
    -------
    float
        The difference between the RMSD between the true position of the ligand in the
        pocket and its position derived from the solution of the QUBO problem and the
        RMSD between the position of the ligand when mapped to the closest point in the
        pocketgrid and its position derived from the solution of the QUBO problem.
    """
    rmsd = 0.0
    rmsd_lower_bound = 0.0
    list_ligand_nodes = list(ligand_graph.nodes())
    list_pocketgrid_nodes = list(pocketgrid_graph.nodes())
    
    pocketgrid_points = np.array([pocketgrid_graph.nodes[v]['coords'] for v in list_pocketgrid_nodes])
    
    for (i, j) in mapping:
        ligand_coords = ligand_graph.nodes[list_ligand_nodes[i]]['coords']
        pocketgrid_coords = pocketgrid_graph.nodes[list_pocketgrid_nodes[j]]['coords']
        closest_point = pocketgrid_points[np.argmin(np.linalg.norm(pocketgrid_points - ligand_coords, axis=1))]
        # closest_point = pocketgrid_points[np.argmin(np.sum(np.abs(pocketgrid_points - ligand_coords), axis=1))]
        rmsd += np.linalg.norm(ligand_coords - pocketgrid_coords)**2
        rmsd_lower_bound += np.linalg.norm(ligand_coords - closest_point)**2
    rmsd = np.sqrt(rmsd / len(mapping))
    rmsd_lower_bound = np.sqrt(rmsd_lower_bound / len(mapping))
    return rmsd - rmsd_lower_bound

def compute_volume_overlap_ratio(
    ligand_graph,
    pocketgrid_graph,
    mapping,
    grid_resolution=50
):
    """
    Compute the volume overlap ratio between the true position of the ligand in the pocket 
    when mapped to the closest point in the pocketgrid and its position derived from
    the solution of the QUBO problem.
    
    Parameters
    ----------
    ligand_graph : nx.Graph
        Graph representation of the ligand molecule.
    pocketgrid_graph : nx.Graph
        Graph representation of the pocketgrid.
    list[tuple]
        List of tuples representing the pairs (i, j) corresponding to the mapping
        between the atoms of the ligand and the pocketgrid.
    grid_resolution : int
        Resolution of the grid to compute the volume overlap.
    
    Returns
    -------
    float
        The volume overlap ratio between the true position of the ligand in the
        pocket when mapped to the closest point in the pocketgrid and its position
        derived from the solution of the QUBO problem.
    """
    list_pocketgrid_nodes = list(pocketgrid_graph.nodes())
    ligand_atom_types = [ligand_graph.nodes[v]['atom_type'] for v in ligand_graph.nodes()]
    ligand_coords = np.array([ligand_graph.nodes[v]['coords'] for v in ligand_graph.nodes()])
    pocketgrid_coords = np.array([pocketgrid_graph.nodes[v]['coords'] for v in pocketgrid_graph.nodes()])
    
    distances = np.linalg.norm(pocketgrid_coords[:, np.newaxis] - ligand_coords, axis=2)
    true_mapping_coords = np.array([
        pocketgrid_coords[np.argmin(distances[:, i]), :] for i in range(len(ligand_coords))
    ])
    mapping_dict = {i: j for i, j in mapping}
    
    found_mapping_coords = np.array([
        pocketgrid_coords[mapping_dict[i]]
        for i in range(len(ligand_coords))
    ])
    
    min_coords = np.array([
        np.min([np.min(true_mapping_coords[:, 0]), np.min(found_mapping_coords[:, 0])]) - max(ATOMIC_RADII.values()),
        np.min([np.min(true_mapping_coords[:, 1]), np.min(found_mapping_coords[:, 1])]) - max(ATOMIC_RADII.values()),
        np.min([np.min(true_mapping_coords[:, 2]), np.min(found_mapping_coords[:, 2])]) - max(ATOMIC_RADII.values())
    ])
    
    max_coords = np.array([
        np.max([np.max(true_mapping_coords[:, 0]), np.max(found_mapping_coords[:, 0])]) + max(ATOMIC_RADII.values()),
        np.max([np.max(true_mapping_coords[:, 1]), np.max(found_mapping_coords[:, 1])]) + max(ATOMIC_RADII.values()),
        np.max([np.max(true_mapping_coords[:, 2]), np.max(found_mapping_coords[:, 2])]) + max(ATOMIC_RADII.values()) 
    ])
    
    # Generate a grid of points
    x = np.linspace(min_coords[0], max_coords[0], grid_resolution)
    y = np.linspace(min_coords[1], max_coords[1], grid_resolution)
    z = np.linspace(min_coords[2], max_coords[2], grid_resolution)
    
    # Compute the midpoints between the grid points along each axis
    x_centers = (x[:-1] + x[1:]) / 2
    y_centers = (y[:-1] + y[1:]) / 2
    z_centers = (z[:-1] + z[1:]) / 2

    # Generate all combinations of the center coordinates
    x_grid, y_grid, z_grid = np.meshgrid(x_centers, y_centers, z_centers, indexing='ij')

    # Combine the center coordinates into a single array of shape (num_cubes, 3)
    cube_centers = np.vstack([x_grid.ravel(), y_grid.ravel(), z_grid.ravel()]).T
    
    # Compute volume overlap
    volume = 0
    volume_overlap = 0
    for center in cube_centers:
        for i in range(len(ligand_coords)):
            if (
                np.linalg.norm(center - true_mapping_coords[i]) <= ATOMIC_RADII.get(ligand_atom_types[i], 0.0)
            ):
                volume += 1
                for j in range(len(ligand_coords)):
                    if (
                        np.linalg.norm(center - found_mapping_coords[j]) <= ATOMIC_RADII.get(ligand_atom_types[j], 0.0)
                    ):
                        volume_overlap += 1
                        break
                break
    
    return volume_overlap / volume

def compute_barycenter(points):
    """Calculate the barycenter (centroid) of a set of points."""
    return np.mean(points, axis=0)

def compute_main_axis(points):
    """Calculate the main axis (direction of maximum variance) of a set of points."""
    shifted_points = points - compute_barycenter(points)    # Shift points by subtracting barycenter
    cov_matrix = np.cov(shifted_points, rowvar=False)       # Covariance matrix
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)  # Eigenvalue decomposition
    main_axis = eigenvectors[:, np.argmax(eigenvalues)]     # Eigenvector with largest eigenvalue
    return main_axis

def compute_rotation_matrix(axis1, axis2):
    """Compute the rotation matrix to align axis1 with axis2 using Rodrigues' formula."""
    axis1 = axis1 / np.linalg.norm(axis1)
    axis2 = axis2 / np.linalg.norm(axis2)
    
    # Compute the cross product (axis of rotation)
    v = np.cross(axis1, axis2)
    v = v / np.linalg.norm(v)
    
    # Compute the dot product (cosine of the angle)
    c = np.dot(axis1, axis2)
    
    # Compute the angle of rotation
    angle = np.arccos(c)
    
    # If the angle is 0, no rotation is needed
    if angle == 0:
        return np.eye(3)
    
    # Skew-symmetric matrix of the rotation axis
    v_skew = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])
    
    # Rodrigues' rotation matrix
    rotation_matrix = np.eye(3) + np.sin(angle) * v_skew + (1 - np.cos(angle)) * np.dot(v_skew, v_skew)
    return rotation_matrix

def align_sets(set1, set2):
    """Align the main axis of set1 to the main axis of set2."""
    # Compute the main axis for both sets
    axis1 = compute_main_axis(set1)
    axis2 = compute_main_axis(set2)
    
    # Compute the rotation matrix to align axis1 with axis2
    rotation_matrix = compute_rotation_matrix(axis1, axis2)
    
    # Apply the rotation to the first set of points
    set1_barycenter = compute_barycenter(set1)
    set1_aligned = np.dot(set1 - set1_barycenter, rotation_matrix) + compute_barycenter(set2)  # Translate to set2's barycenter
    return set1_aligned

def export_sdf(
    mol2_file,
    sdf_file,
    removeHs=False
):
    """
    Read a molecule from a Mol2 file and export it to an SDF file.
    
    Parameters
    ----------
    mol2_file : str
        Path to the Mol2 file.
    sdf_file : str
        Path to the SDF file.
    removeHs : bool (optional)
        Whether to remove hydrogens from the molecule.
    """
    # Read the molecule from the Mol2 file
    mol = Chem.MolFromMol2File(mol2_file, removeHs=removeHs)
    if mol is None:
        raise ValueError("Invalid Mol2 file.")
    
    # Export the molecule to an SDF file
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()

def export_sdf_with_modified_coordinates(
    mol2_file,
    sdf_file,
    new_coords,
    removeHs=False
):
    """
    Read a molecule from a Mol2 file, modify its coordinates, and export it to an SDF file.
    
    Parameters
    ----------
    mol2_file : str
        Path to the Mol2 file.
    sdf_file : str
        Path to the SDF file.
    new_coords : numpy.ndarray
        New coordinates to assign to the atoms.
    removeHs : bool (optional)
        Whether to remove hydrogens from the molecule.
    """
    # Read the molecule from the Mol2 file
    mol = Chem.MolFromMol2File(mol2_file, removeHs=removeHs)
    if mol is None:
        raise ValueError("Invalid Mol2 file.")
    
    # Assign the new coordinates to the atoms
    conf = mol.GetConformer()
    for i, coords in enumerate(new_coords):
        conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(float(coords[0]), float(coords[1]), float(coords[2])))
    
    # Export the molecule to an SDF file
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()