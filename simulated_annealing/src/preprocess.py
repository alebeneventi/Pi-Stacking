import os 
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from parameters import LJ_PARAMS_FILE_PATH
from parameters import LIGAND_PATH1, POCKET_PATH, POCKETGRID_PDB_PATH

def generate_vdw_parameters_csv(ligand_path):
    """
    Generates a CSV file with Van der Waals parameters for each atom type 
    in the ligand molecule, relying on MMFF94. Adds a 'RingCenter' entry
    only if the ligand contains at least one aromatic benzene ring.
    """
    # Load molecule and MMFF properties
    ligand = Chem.MolFromMol2File(ligand_path, removeHs=False)
    mmff_props = AllChem.MMFFGetMoleculeProperties(
        ligand, mmffVariant='MMFF94'
    )

    # Detect presence of aromatic benzene rings
    rings = ligand.GetRingInfo().AtomRings()
    has_benzene = any(
        len(ring) == 6 and all(
            ligand.GetAtomWithIdx(i).GetSymbol() == 'C' and ligand.GetAtomWithIdx(i).GetIsAromatic()
            for i in ring
        ) for ring in rings
    )

    # Collect unique atom types and VdW parameters
    seen_types = set()
    records = []
    for i in range(ligand.GetNumAtoms()):
        atom_type_char = ligand.GetAtomWithIdx(i).GetSymbol()
        atom_type_num  = mmff_props.GetMMFFAtomType(i)
        vdw_param = mmff_props.GetMMFFVdWParams(i, i)
        if atom_type_num not in seen_types:
            seen_types.add(atom_type_num)
            records.append({
                'AtomType': atom_type_char,
                'AtomTypeMMFF': atom_type_num,
                'Rmin': vdw_param[0],
                'Eps': vdw_param[1]
            })

    # Create DataFrame and sort
    df = pd.DataFrame(records)
    df = df.sort_values(by=['AtomType', 'AtomTypeMMFF'], ignore_index=True)

    # Add RingCenter row only if benzene present
    if has_benzene and 'GhostAtom' not in df['AtomType'].values:
        ring_center_row = pd.DataFrame([{  
            'AtomType': 'GhostAtom',
            'AtomTypeMMFF': -1,
            'Rmin': 0.0,
            'Eps': 0.0
        }])
        df = pd.concat([df, ring_center_row], ignore_index=True)

    # Write to CSV
    os.makedirs(os.path.dirname(LJ_PARAMS_FILE_PATH), exist_ok=True)
    df.to_csv(LJ_PARAMS_FILE_PATH, index=False)

def main():
    generate_vdw_parameters_csv(LIGAND_PATH1)
    import pdb_utils as pdb
    pdb.build_csv_file_V2(POCKETGRID_PDB_PATH, POCKET_PATH)

if __name__ == "__main__":
    main()