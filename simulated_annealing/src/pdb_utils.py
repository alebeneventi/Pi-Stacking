import os
import numpy as np
import pandas as pd
from Bio import PDB
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from parameters import LJ_PARAMS_FILE_PATH, VDW_CUTOFF, USE_LJ_12_6
from mol2_utils import get_coordinates_from_mol2_path
from utils import get_coordinates_from_mol_object
from utils import find_substructure_matches, compute_angle
from utils import pistacking_energy

def get_coordinates_from_pdb(pdb_file):
    """
    Extract the coordinates of the atoms from a PDB file.

    Parameters
    ----------
    pdb_file : str
        The path to the PDB file.
    
    Returns
    -------
    ndarray
        A NumPy array of coordinates.
    """
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Parse the PDB file
    structure = parser.get_structure(pdb_file, pdb_file)

    # Count total number of atoms
    num_atoms = sum(len(residue) for model in structure for chain in model for residue in chain)

    # Preallocate a NumPy array for the coordinates
    coordinates = np.zeros((num_atoms, 3), dtype=float)

    # Extract the coordinates
    index = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coordinates[index] = atom.get_coord()
                    index += 1
    
    return coordinates

def get_atom_type_from_pdb(pdb_file):
    """
    Extract the atom type of the atoms from a PDB file.

    Parameters
    ----------
    pdb_file : str
        The path to the PDB file.
    
    Returns
    -------
    list
        A list of atom types.
    """
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Parse a PDB file
    structure = parser.get_structure(pdb_file, pdb_file)

    # Extract the atom types
    atom_types = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_types.append(atom.element)
    
    return atom_types

def get_csv_data_from_pbd(pdb_file):
    """
    Initialize the CSV data with the atom IDs and coordinates from a PDB file.
    
    Parameters
    ----------
    pdb_file : str
        The path to the PDB file.
    
    Returns
    -------
    csv_data : dict
        A dictionary to store the CSV data.
    """
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Parse a PDB file
    structure = parser.get_structure(pdb_file, pdb_file)

    # Extract the atom IDS and the coordinates
    atom_IDs = []
    x_coords = []
    y_coords = []
    z_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_IDs.append(atom.get_serial_number())
                    x, y, z = atom.get_coord()
                    x_coords.append(float(x))
                    y_coords.append(float(y))
                    z_coords.append(float(z))

    # Create the CSV data
    csv_data = {
        'atom_ID': atom_IDs,
        'x': x_coords,
        'y': y_coords,
        'z': z_coords
    }

    return csv_data

def get_augmented_csv_data_from_pbd(pdb_file, ligand_path):
    """
    Initialize the augmented CSV data with the coordinates from a PDB file and the
    ligand coordinates from a mol2 file.
    
    Parameters
    ----------
    pdb_file : str
        The path to the PDB file.
    ligand_path : str
        The path to the mol2 file.
    
    Returns
    -------
    csv_data : dict
        A dictionary to store the CSV data.
    """
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Parse a PDB file
    structure = parser.get_structure(pdb_file, pdb_file)

    # Extract the atom IDS and the coordinates
    atom_IDs = []
    x_coords = []
    y_coords = []
    z_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_IDs.append(atom.get_serial_number())
                    x, y, z = atom.get_coord()
                    x_coords.append(float(x))
                    y_coords.append(float(y))
                    z_coords.append(float(z))
    
    next_atom_ID = max(atom_IDs) + 1
    
    # Extract the ligand coordinates
    ligand_coords = get_coordinates_from_mol2_path(ligand_path)
    for coord in ligand_coords:
        atom_IDs.append(next_atom_ID)
        x_coords.append(float(coord[0]))
        y_coords.append(float(coord[1]))
        z_coords.append(float(coord[2]))
        next_atom_ID += 1

    # Create the CSV data
    csv_data = {
        'atom_ID': atom_IDs,
        'x': x_coords,
        'y': y_coords,
        'z': z_coords
    }

    return csv_data

def vdw_12_6_lj(epsilon, r_min, r):
    """
    Calculate the 12-6 Lennard-Jones potential energy between two atoms.

    Parameters
    ----------
    epsilon : float
        The depth of the potential well.
    r_min : float
        The distance at which the potential energy is minimized.
    r : float
        The distance between the two atoms.
    
    Returns
    -------
    float
        The potential energy between the two atoms
    """
    return epsilon * ((r_min / r) ** 12 - 2 * (r_min / r) ** 6)

def vdw_8_4_lj(epsilon, r_min, r):
    """
    Calculate the 8-4 Lennard-Jones potential energy between two atoms.

    Parameters
    ----------
    epsilon : float
        The depth of the potential well.
    r_min : float
        The distance at which the potential energy is minimized.
    r : float
        The distance between the two atoms.
    
    Returns
    -------
    float
        The potential energy between the two atoms
    """
    return epsilon * ((r_min / r) ** 8 - 2 * (r_min / r) ** 4)

def get_element_to_exclude(mol):
    """
    Get the list of elements to exclude from the molecule.
    
    Parameters
    ----------
    mol : Chem.Mol
        The molecule to preprocess.
    
    Returns
    -------
    list
        The list of elements to exclude.
    """
    elements_to_exclude = []
    RDLogger.DisableLog('rdApp.*')
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol != 'H':
            try:
                atom_mol = Chem.MolFromSmiles(symbol)
                atom_mol = Chem.AddHs(atom_mol)
                if (symbol not in elements_to_exclude) and (not AllChem.MMFFHasAllMoleculeParams(atom_mol)):
                    elements_to_exclude.append(symbol) 
            except Exception as e:
                if symbol not in elements_to_exclude:
                    elements_to_exclude.append(symbol)
    RDLogger.EnableLog('rdApp.*')
    return elements_to_exclude

def pocket_preprocess(mol, elements_to_remove):
    """
    Preprocess the molecule by removing the atoms with the specified elements.
    
    Parameters
    ----------
    mol : Chem.Mol
        The molecule to preprocess.
    elements_to_remove : list
        The list of elements to remove.
    
    Returns
    -------
    Chem.Mol
        The preprocessed molecule.
    """
    atoms_to_remove = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in elements_to_remove:
            atoms_to_remove.append(atom.GetIdx())
    for idx in reversed(atoms_to_remove):
        mol = Chem.RWMol(mol)
        mol.RemoveAtom(idx)
    return mol

def vdw_12_6_lj_energy_from_pdb(
    pocket_path,
    pocketgrid_coordinates,
    lj_parameters_csv_path
):
    """
    Calculate the van der Waals energy for each pocketgrid point using the 12-6 Lennard-Jones 
    potential.
    
    Parameters
    ----------
    pocket_path : str
        Path to the PDB file containing the pocket.
    pocketgrid_coordinates : ndarray
        A NumPy array of coordinates for the pocketgrid points.
    lj_parameters_csv_path : str
        Path to the CSV file containing Lennard-Jones parameters.
    
    Returns
    -------
    list[list]
        A list of lists containing van der Waals energy values computed for each pocketgrid point
        and each set of Lennard-Jones parameters.
    """
    # Load the pocket and preprocess it
    pocket = Chem.MolFromPDBFile(pocket_path, removeHs=False)
    elements_to_exclude = get_element_to_exclude(pocket)
    pocket = pocket_preprocess(pocket, elements_to_exclude)

    pocket = Chem.AddHs(pocket, addCoords=True)
    
    # Retrieve the MMFF94 properties for the pocket
    RDLogger.DisableLog('rdApp.*')
    pocket_mmff_props = AllChem.MMFFGetMoleculeProperties(pocket, mmffVariant='MMFF94')
    pocket_coordinates = get_coordinates_from_mol_object(pocket)
    RDLogger.EnableLog('rdApp.*')
    
    lj_parameters = pd.read_csv(lj_parameters_csv_path)
    vdw_energy_values = []
    
    if VDW_CUTOFF is None:
        for pocketgrid_coord in pocketgrid_coordinates:
            temp_list = []
            
            distances = np.linalg.norm(pocket_coordinates - pocketgrid_coord, axis=1)
            epsilon = [
                pocket_mmff_props.GetMMFFVdWParams(i, i)[1]
                for i in range(pocket_coordinates.shape[0])
            ]
            r_min = [
                pocket_mmff_props.GetMMFFVdWParams(i, i)[0]
                for i in range(pocket_coordinates.shape[0])
            ]
            for _, lj_param in lj_parameters.iterrows():
                epsilon_combined = [
                    (epsilon[i] * lj_param['Eps'])**(1/2)
                    for i in range(len(epsilon))
                ]
                r_min_combined = [
                    0.5 * (r_min[i] + lj_param['Rmin'])
                    for i in range(len(r_min))
                ]
                temp_list.append(
                    np.sum(
                        vdw_12_6_lj(epsilon_combined[i], r_min_combined[i], distance)
                        for i, distance in enumerate(distances)
                    )
                )
            vdw_energy_values.append(temp_list)
    else:
        for pocketgrid_coord in pocketgrid_coordinates:
            temp_list = []
            
            distances = np.linalg.norm(pocket_coordinates - pocketgrid_coord, axis=1)
            epsilon = [
                pocket_mmff_props.GetMMFFVdWParams(i, i)[1]
                for i in range(pocket_coordinates.shape[0])
                if distances[i] <= VDW_CUTOFF
            ]
            r_min = [
                pocket_mmff_props.GetMMFFVdWParams(i, i)[0]
                for i in range(pocket_coordinates.shape[0])
                if distances[i] <= VDW_CUTOFF
            ]
            distances = distances[distances <= VDW_CUTOFF]
            for _, lj_param in lj_parameters.iterrows():
                epsilon_combined = [
                    (epsilon[i] * lj_param['Eps'])**(1/2)
                    for i in range(len(epsilon))
                ]
                r_min_combined = [
                    0.5 * (r_min[i] + lj_param['Rmin'])
                    for i in range(len(r_min))
                ]
                temp_list.append(
                    np.sum(
                        vdw_12_6_lj(epsilon_combined[i], r_min_combined[i], distance)
                        for i, distance in enumerate(distances)
                    )
                )
            vdw_energy_values.append(temp_list)
    
    return vdw_energy_values 

def vdw_8_4_lj_energy_from_pdb(
    pocket_path,
    pocketgrid_coordinates,
    lj_parameters_csv_path
):
    """
    Calculate the van der Waals energy for each pocketgrid point using the 8-4 Lennard-Jones 
    potential.
    
    Parameters
    ----------
    pocket_path : str
        Path to the PDB file containing the pocket.
    pocketgrid_coordinates : ndarray
        A NumPy array of coordinates for the pocketgrid points.
    lj_parameters_csv_path : str
        Path to the CSV file containing Lennard-Jones parameters.
    
    Returns
    -------
    list[list]
        A list of lists containing van der Waals energy values computed for each pocketgrid point
        and each set of Lennard-Jones parameters.
    """
    # Load the pocket and preprocess it
    pocket = Chem.MolFromPDBFile(pocket_path, removeHs=False)
    elements_to_exclude = get_element_to_exclude(pocket)
    pocket = pocket_preprocess(pocket, elements_to_exclude)

    
    # Retrieve the MMFF94 properties for the pocket
    RDLogger.DisableLog('rdApp.*')
    pocket_mmff_props = AllChem.MMFFGetMoleculeProperties(pocket, mmffVariant='MMFF94')
    pocket_coordinates = get_coordinates_from_mol_object(pocket)
    RDLogger.EnableLog('rdApp.*')
    
    lj_parameters = pd.read_csv(lj_parameters_csv_path)
    vdw_energy_values = []
    
    if VDW_CUTOFF is None:
        for pocketgrid_coord in pocketgrid_coordinates:
            temp_list = []
            
            distances = np.linalg.norm(pocket_coordinates - pocketgrid_coord, axis=1)
            epsilon = [
                pocket_mmff_props.GetMMFFVdWParams(i, i)[1]
                for i in range(pocket_coordinates.shape[0])
            ]
            r_min = [
                pocket_mmff_props.GetMMFFVdWParams(i, i)[0]
                for i in range(pocket_coordinates.shape[0])
            ]
            for _, lj_param in lj_parameters.iterrows():
                epsilon_combined = [
                    (epsilon[i] * lj_param['Eps'])**(1/2)
                    for i in range(len(epsilon))
                ]
                r_min_combined = [
                    0.5 * (r_min[i] + lj_param['Rmin'])
                    for i in range(len(r_min))
                ]
                temp_list.append(
                    np.sum(
                        vdw_8_4_lj(epsilon_combined[i], r_min_combined[i], distance)
                        for i, distance in enumerate(distances)
                    )
                )
            vdw_energy_values.append(temp_list)
    else:
        for pocketgrid_coord in pocketgrid_coordinates:
            temp_list = []
            
            distances = np.linalg.norm(pocket_coordinates - pocketgrid_coord, axis=1)
            epsilon = [
                pocket_mmff_props.GetMMFFVdWParams(i, i)[1]
                for i in range(pocket_coordinates.shape[0])
                if distances[i] <= VDW_CUTOFF
            ]
            r_min = [
                pocket_mmff_props.GetMMFFVdWParams(i, i)[0]
                for i in range(pocket_coordinates.shape[0])
                if distances[i] <= VDW_CUTOFF
            ]
            distances = distances[distances <= VDW_CUTOFF]
            for _, lj_param in lj_parameters.iterrows():
                epsilon_combined = [
                    (epsilon[i] * lj_param['Eps'])**(1/2)
                    for i in range(len(epsilon))
                ]
                r_min_combined = [
                    0.5 * (r_min[i] + lj_param['Rmin'])
                    for i in range(len(r_min))
                ]
                temp_list.append(
                    np.sum(
                        vdw_8_4_lj(epsilon_combined[i], r_min_combined[i], distance)
                        for i, distance in enumerate(distances)
                    )
                )
            vdw_energy_values.append(temp_list)
    
    return vdw_energy_values 

    
def compute_electric_potential(
    positions, 
    charges, 
    observation_point
):
    """
    Compute the electric potential at a given observation point due to a set of charges.
    
    Parameters
    ----------
    positions : np.ndarray
        A 2D array of shape (n, 3) containing the coordinates of the charges.
    charges : np.ndarray
        A 1D array of shape (n,) containing the charges.
    observation_point : np.ndarray
        A 1D array of shape (3,) containing the coordinates of the observation point.
    
    Returns
    -------
    float
        The electric potential at the observation point.
    """
    # Constants
    epsilon_0 = 8.854187817e-12 # Vacuum permittivity in units of F/m
    conversion_factor = 1e-10   # Convert from Angstrom to meters

    # Compute the distances between the observation point and the charge positions
    distances = np.linalg.norm(positions - observation_point, axis=1)
    distances_m = distances * conversion_factor

    # Compute the electric potential
    potential = np.sum([q / r for q, r in zip(charges, distances_m)], axis=0) / (4 * np.pi * epsilon_0)
    return potential

def get_pocketgrid_electric_potential(
    pdb_pocket,
    pocketgrid_coordinates
):
    """
    Compute the electric potential at each pocketgrid point due to the charges in the pocket.
    
    Parameters
    ----------
    pdb_pocket : str
        The path to the PDB file with the pocket.
    pocketgrid_coordinates : np.ndarray
        A 2D array of shape (n, 3) containing coordinates for the pocketgrid points.
    
    Returns
    -------
    list[float]
        A list of electric potentials at each pocketgrid point.
    """
    elementary_charge = 1.60217662e-19 # Elementary charge in Coulombs
    
    pocket = Chem.MolFromPDBFile(pdb_pocket, removeHs=False)
    if pocket is None:
        raise ValueError("Failed to load molecule from the given .pdb file.")
    elements_to_exclude = get_element_to_exclude(pocket)
    pocket = pocket_preprocess(pocket, elements_to_exclude)
    pocket_coordinates = get_coordinates_from_mol_object(pocket)
    
    pocket_with_Hs = Chem.AddHs(pocket, addCoords=True)
    mmff_props = AllChem.MMFFGetMoleculeProperties(pocket_with_Hs, mmffVariant='MMFF94')
    
    #pocket_charges = np.zeros(len(pocket_coordinates))
    #for atom in pocket.GetAtoms():
    #    charge = mmff_props.GetMMFFPartialCharge(atom.GetIdx()) * elementary_charge
        #charge = -0.062269 * elementary_charge
    #    atom_idx = atom.GetIdx()
    #    pocket_charges[atom_idx] = charge

    pocket_coordinates = get_coordinates_from_mol_object(pocket_with_Hs)
    pocket_charges = np.zeros(len(pocket_coordinates))
    for atom in pocket_with_Hs.GetAtoms():
        charge = mmff_props.GetMMFFPartialCharge(atom.GetIdx()) * elementary_charge
        atom_idx = atom.GetIdx()
        pocket_charges[atom_idx] = charge
    
    electric_potentials = []
    for pocketgrid_coord in pocketgrid_coordinates:
        electric_potentials.append(
            compute_electric_potential(pocket_coordinates, pocket_charges, pocketgrid_coord)
        )
    
    return electric_potentials



def get_pocketgrid_pistacking(
    pdb_pocket,
    pocketgrid_coordinates,
    z_range=(3.2, 4),
    r_range=(1, 3.0),
    cv_range=(1, 2),
    angle_limit=15.0
):

    # from your_module import pistacking_energy

    # Load the pocket and identify the rings
    pocket = Chem.MolFromPDBFile(pdb_pocket, removeHs=False)
    if pocket is None:
        raise ValueError(f"Cannot load PDB from {pdb_pocket!r}")
    rings = pocket.GetRingInfo().AtomRings()

    # Compute centroids and normals for hexagonal aromatic rings
    ring_planes = []
    conf = pocket.GetConformer()
    for ring in rings:
        if len(ring) != 6:
            continue
        if not all(
            pocket.GetAtomWithIdx(i).GetSymbol() == 'C' and pocket.GetAtomWithIdx(i).GetIsAromatic()
            for i in ring
        ):
            continue
        coords = np.array([[conf.GetAtomPosition(i).x,
                            conf.GetAtomPosition(i).y,
                            conf.GetAtomPosition(i).z]
                           for i in ring])
        centroid = coords.mean(axis=0)
        v1 = coords[1] - coords[0]
        v2 = coords[2] - coords[0]
        normal = np.cross(v1, v2)
        normal /= np.linalg.norm(normal)
        ring_planes.append((centroid, normal))

    m = pocketgrid_coordinates.shape[0]
    # labels now come as a list of lists of tuples (V_idx, energy)
    labels = [[] for _ in range(m)]

    # First pass: assign each C to possible planes
    from collections import defaultdict
    C_planes = defaultdict(list)
    for idx, point in enumerate(pocketgrid_coordinates):
        for centroid, normal in ring_planes:
            diff = point - centroid
            dz = abs(np.dot(diff, normal))
            in_plane = diff - np.dot(diff, normal) * normal
            dr = np.linalg.norm(in_plane)
            if z_range[0] <= dz <= z_range[1] and r_range[0] <= dr <= r_range[1]:
                C_planes[idx].append((centroid, normal))

    # Second pass: for each C and for each plane, compute the pairs (V_idx, energy)
    for C_idx, planes in C_planes.items():
        C_point = pocketgrid_coordinates[C_idx]
        energy_list = []
        for centroid, normal in planes:
            diff_C = C_point - centroid
            dz_CC = abs(np.dot(diff_C, normal))
            in_plane_C = diff_C - np.dot(diff_C, normal) * normal
            dr_CC = np.linalg.norm(in_plane_C)

            # define u and v axes in-plane
            u_hat = in_plane_C / dr_CC
            v_hat = np.cross(normal, u_hat)
            v_hat /= np.linalg.norm(v_hat)

            # for each V compute theta_u, theta_v and add if valid
            for V_idx, V_point in enumerate(pocketgrid_coordinates):
                vec = V_point - C_point
                dist = np.linalg.norm(vec)
                if not (cv_range[0] <= dist <= cv_range[1]):
                    continue
                
                if np.dot(vec, normal) < 0:
                    vec = -vec  # Ensure the vector points toward the ring
                c_n = np.dot(vec, normal)
                c_u = np.dot(vec, u_hat)
                c_v = np.dot(vec, v_hat)

                theta_u = abs(np.degrees(np.arctan2(-c_v, c_n)))
                theta_v = np.degrees(np.arctan2(c_u, c_n))

                if theta_u <= angle_limit and -angle_limit <= theta_v <= angle_limit:
                    energy = pistacking_energy(dr_CC, dz_CC, theta_u, theta_v)
                    energy_list.append((V_idx, energy))
                    #energy_list.append((V_idx, dr_CC, dz_CC, theta_u, theta_v))

        labels[C_idx] = energy_list

    return labels



def build_csv_file(pdb_pocketgrid, pdb_pocket, output_dir):
    """
    Build a CSV file with all the major pocketgrid data.

    Parameters
    ----------
    pdb_pocketgrid : str
        The path to the PDB file with the pocketgrid.
    pdb_pocket : str
        The path to the PDB file with the pocket.
    output_dir : str
        The path to the output directory.
    """
    # Generate the output CSV file path.
    os.makedirs(output_dir, exist_ok=True)
    pdb_path_parts = pdb_pocketgrid.split('/')
    csv_output_path = os.path.join(output_dir, f"{pdb_path_parts[-1].replace('.pdb', '.csv')}")

    # Initialize the CSV data
    csv_data = get_csv_data_from_pbd(pdb_pocketgrid)

    # Get the coordinates of the pocketgrid and the pocket
    pocketgrid_coordinates = get_coordinates_from_pdb(pdb_pocketgrid)
    pocket_coordinates     = get_coordinates_from_pdb(pdb_pocket)

    # Add the van der Waals energy to the CSV data
    lj_parameters_csv_path = LJ_PARAMS_FILE_PATH
    if USE_LJ_12_6:
        csv_data['vdw_energy'] = vdw_12_6_lj_energy_from_pdb(
            pdb_pocket,
            pocketgrid_coordinates,
            lj_parameters_csv_path
        )
    else:
        csv_data['vdw_energy'] = vdw_8_4_lj_energy_from_pdb(
            pdb_pocket,
            pocketgrid_coordinates,
            lj_parameters_csv_path
        )


    
    # Add the electric potential to the CSV data
    csv_data['electric_potential'] = get_pocketgrid_electric_potential(
        pdb_pocket,
        pocketgrid_coordinates
    )
    

    # Create DataFrame
    df = pd.DataFrame(csv_data)

    # Save to CSV
    df.to_csv(csv_output_path, index=False)

def build_csv_file_V2(pdb_pocketgrid, pdb_pocket):
    """
    Build a CSV file with all the major pocketgrid data.

    Parameters
    ----------
    pdb_pocketgrid : str
        The path to the PDB file with the pocketgrid.
    pdb_pocket : str
        The path to the PDB file with the pocket.
    """
    # Generate the output CSV file path.
    if pdb_pocketgrid.endswith('.pdb'):
        csv_output_path = pdb_pocketgrid[:-4] + '.csv'
    else:
        raise ValueError("Pocketgrid file must be a PDB file.")

    # Initialize the CSV data
    csv_data = get_csv_data_from_pbd(pdb_pocketgrid)

    # Get the coordinates of the pocketgrid and the pocket
    pocketgrid_coordinates = get_coordinates_from_pdb(pdb_pocketgrid)
    pocket_coordinates     = get_coordinates_from_pdb(pdb_pocket)

    # Add the van der Waals energy to the CSV data
    lj_parameters_csv_path = LJ_PARAMS_FILE_PATH
    if USE_LJ_12_6:
        csv_data['vdw_energy'] = vdw_12_6_lj_energy_from_pdb(
            pdb_pocket,
            pocketgrid_coordinates,
            lj_parameters_csv_path
        )
    else:
        csv_data['vdw_energy'] = vdw_8_4_lj_energy_from_pdb(
            pdb_pocket,
            pocketgrid_coordinates,
            lj_parameters_csv_path
        )
    
    
    # Add the electric potential to the CSV data
    csv_data['electric_potential'] = get_pocketgrid_electric_potential(
        pdb_pocket,
        pocketgrid_coordinates
    )

    csv_data['pistacking'] = get_pocketgrid_pistacking(
        pdb_pocket,
        pocketgrid_coordinates
    )
    

    # Create DataFrame
    df = pd.DataFrame(csv_data)

    # Save to CSV
    df.to_csv(csv_output_path, index=False)
    
def build_csv_augmented_file(
    pdb_pocketgrid,
    pdb_pocket,
    ligand_path,
    output_dir
):
    """
    Build a CSV file with all the major pocketgrid data and augmented with the ligand
    coordinates.

    Parameters
    ----------
    pdb_pocketgrid : str
        The path to the PDB file with the pocketgrid.
    pdb_pocket : str
        The path to the PDB file with the pocket.
    ligand_path : str
        The path to the ligand mol2 file.
    output_dir : str
        The path to the output directory.
    """
    # Generate the output CSV file path.
    os.makedirs(output_dir, exist_ok=True)
    pdb_path_parts = pdb_pocketgrid.split('/')
    csv_output_path = os.path.join(output_dir, f"{pdb_path_parts[-1].replace('.pdb', '_augmented.csv')}")
    
    # Initialize the augmented CSV data
    csv_data = get_augmented_csv_data_from_pbd(pdb_pocketgrid, ligand_path)
    
    # Get the coordinates
    pocket_coordinates = get_coordinates_from_pdb(pdb_pocket)
    ligand_coordinates = get_coordinates_from_mol2_path(ligand_path)
    pocketgrid_coordinates = np.vstack(
        (get_coordinates_from_pdb(pdb_pocketgrid), ligand_coordinates)
    )
    
    # Add the van der Waals energy to the CSV data
    lj_parameters_csv_path = LJ_PARAMS_FILE_PATH
    if USE_LJ_12_6:
        csv_data['vdw_energy'] = vdw_12_6_lj_energy_from_pdb(
            pdb_pocket,
            pocketgrid_coordinates,
            lj_parameters_csv_path
        )
    else:
        csv_data['vdw_energy'] = vdw_8_4_lj_energy_from_pdb(
            pdb_pocket,
            pocketgrid_coordinates,
            lj_parameters_csv_path
        )
    
    
    # Add the electric potential to the CSV data
    csv_data['electric_potential'] = get_pocketgrid_electric_potential(
        pdb_pocket,
        pocketgrid_coordinates
    )

    

    # Create DataFrame
    df = pd.DataFrame(csv_data)

    # Save to CSV
    df.to_csv(csv_output_path, index=False)

#print(get_pocketgrid_pistacking('data/pockets/smallMol_max6HA_max0RB/1l83/1l83_pocket.pdb',get_coordinates_from_pdb('data/pockets/smallMol_max6HA_max0RB/1l83/1l83_pass_pocketgrid.pdb')))