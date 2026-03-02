from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd

def get_coordinates_from_mol2_path(mol2_path, removeHs=False):
    """
    Get the coordinates of the atoms in a mol2 file.

    Parameters
    ----------
    mol2_path : str
        Path to the mol2 file.
    
    Returns
    -------
    numpy.ndarray
        Array with the coordinates of the atoms.
    """
    RDLogger.DisableLog('rdApp.*')
    mol = Chem.MolFromMol2File(mol2_path, removeHs=removeHs)
    if mol is None:
        raise ValueError("Failed to load molecule from the given .mol2 file.")
    
    conf = mol.GetConformer()
    if conf is None:
        raise ValueError("No conformer found for the molecule.")

    atom_coords = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atom_coords.append([pos.x, pos.y, pos.z])
    atom_coords = np.array(atom_coords)

    RDLogger.EnableLog('rdApp.*')
    return atom_coords

def get_atom_type_from_mol2_path(mol2_path, removeHs=False):
    """
    Get the type of the atoms in a mol2 file.

    Parameters
    ----------
    mol2_path : str
        Path to the mol2 file.
    
    Returns
    -------
    list
        List with the type of each atom.
    """
    mol = Chem.MolFromMol2File(mol2_path, removeHs=removeHs)
    if mol is None:
        raise ValueError("Failed to load molecule from the given .mol2 file.")
    
    atom_types = []
    for atom in mol.GetAtoms():
        atom_types.append(atom.GetSymbol())
    
    return atom_types

def get_dataframe_from_mol2_path(
    mol2_path, 
    lj_parameters_csv_path,
    removeHs=False
):
    """
    Get a pandas DataFrame with the atoms data from a mol2 file.
    
    Parameters
    ----------
    mol2_path : str
        Path to the mol2 file.
    lj_parameters_csv_path : str
        Path to the CSV file with the Lennard-Jones parameters.
    removeHs : bool, optional
        Remove hydrogens from the molecule. The default is False.
    
    Returns
    -------
    pandas.DataFrame
        DataFrame with the atoms data.
    """
    mol = Chem.MolFromMol2File(mol2_path, removeHs=removeHs)
    if mol is None:
        raise ValueError("Failed to load molecule from the given .mol2 file.")
    
    conf = mol.GetConformer()
    if conf is None:
        raise ValueError("No conformer found for the molecule.")
    
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
    lj_parameters = pd.read_csv(lj_parameters_csv_path)
    lj_parameters_atom_types = lj_parameters['AtomTypeMMFF'].tolist()
    
    atom_data = []
    for atom in mol.GetAtoms():
        atom_data.append({
            'atom_ID': atom.GetIdx(),   # Indexing starts from 0
            'atom_type': atom.GetSymbol(),
            'mmff_atom_type': mmff_props.GetMMFFAtomType(atom.GetIdx()),
            'mmff_partial_charge': mmff_props.GetMMFFPartialCharge(atom.GetIdx()),
            'x': conf.GetAtomPosition(atom.GetIdx()).x,
            'y': conf.GetAtomPosition(atom.GetIdx()).y,
            'z': conf.GetAtomPosition(atom.GetIdx()).z,
            'lj_param_idx': lj_parameters_atom_types.index(mmff_props.GetMMFFAtomType(atom.GetIdx()))
        })
    
    df = pd.DataFrame(atom_data)
    
    return df

def get_dataframe_from_mol2_path_with_H_selection(
    mol2_path, 
    H_indices_to_preserve=set()
):
    """
    Get a pandas DataFrame with the atoms data from a mol2 file.
    All the hydrogens in the molecule are removed, except for the ones
    with indices in the given set.
    
    Parameters
    ----------
    mol2_path : str
        Path to the mol2 file.
    H_indices_to_preserve : set, optional
        Set with the indices of the hydrogens to preserve. The default is set().
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with the atoms data.
    """
    mol = Chem.MolFromMol2File(mol2_path, removeHs=False)
    if mol is None:
        raise ValueError("Failed to load molecule from the given .mol2 file.")
    
    conf = mol.GetConformer()
    if conf is None:
        raise ValueError("No conformer found for the molecule.")
    
    atom_data = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H' and atom.GetIdx() not in H_indices_to_preserve:
            continue
        atom_data.append({
            'atom_ID': atom.GetIdx(),   # Indexing starts from 0
            'atom_type': atom.GetSymbol(),
            # 'atomic_number': atom.GetAtomicNum(),
            # 'formal_charge': atom.GetFormalCharge(),
            # 'is_aromatic': atom.GetIsAromatic(),
            'x': conf.GetAtomPosition(atom.GetIdx()).x,
            'y': conf.GetAtomPosition(atom.GetIdx()).y,
            'z': conf.GetAtomPosition(atom.GetIdx()).z
        })
    
    df = pd.DataFrame(atom_data)
    
    return df

def get_bonds_dataframe_from_mol2_path(mol2_path, removeHs=False):
    """
    Get a pandas DataFrame with the bonds data from a mol2 file.
    
    Parameters
    ----------
    mol2_path : str
        Path to the mol2 file.
    removeHs : bool, optional
        Remove hydrogens from the molecule. The default is False.
    
    Returns
    -------
    pandas.DataFrame
        DataFrame with the bonds data.
    """
    mol = Chem.MolFromMol2File(mol2_path, removeHs=removeHs)
    if mol is None:
        raise ValueError("Failed to load molecule from the given .mol2 file.")
    
    bond_data = []
    for bond in mol.GetBonds():
        bond_data.append({
            'bond_ID': bond.GetIdx(),   # Indexing starts from 0
            'atom1_ID': bond.GetBeginAtomIdx(),
            'atom2_ID': bond.GetEndAtomIdx(),
            'bond_type': str(bond.GetBondTypeAsDouble())   # (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)
        })
    
    df = pd.DataFrame(bond_data)
    
    return df

def get_bonds_dataframe_from_mol2_path_with_H_selection(
    mol2_path, 
    H_indices_to_preserve=set()
):
    """
    Get a pandas DataFrame with the bonds data from a mol2 file.
    All the bonds involving hydrogens are removed, except for the ones
    having hydrogens with indices in the given set.
    
    Parameters
    ----------
    mol2_path : str
        Path to the mol2 file.
    H_indices_to_preserve : set, optional
        Set with the indices of the hydrogens to preserve. The default is set().
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with the bonds data.
    """
    mol = Chem.MolFromMol2File(mol2_path, removeHs=False)
    if mol is None:
        raise ValueError("Failed to load molecule from the given .mol2 file.")
    
    bond_data = []
    for bond in mol.GetBonds():
        if (
            mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol() == 'H' and
            bond.GetBeginAtomIdx() not in H_indices_to_preserve
        ):
            continue
        if (
            mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol() == 'H' and
            bond.GetEndAtomIdx() not in H_indices_to_preserve
        ):
            continue
        bond_data.append({
            'bond_ID': bond.GetIdx(),   # Indexing starts from 0
            'atom1_ID': bond.GetBeginAtomIdx(),
            'atom2_ID': bond.GetEndAtomIdx(),
            'bond_type': str(bond.GetBondTypeAsDouble())   # (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)
        })
    
    df = pd.DataFrame(bond_data)
    
    return df