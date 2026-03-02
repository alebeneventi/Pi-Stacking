import numpy as np
import generate_qubo as gq
import graph_utils as gu
import mol2_utils as mu
import utils as ut
from parameters import SEED, REMOVE_HYDROGENS
from parameters import ENCODING, LAMBDA_WEIGHTS, NUM_READS, NUM_SWEEPS
from parameters import LIGAND_PATH1, POCKETGRID_CSV_PATH
from parameters import ADD_BENZENE_NODES

import dimod

ENCODING_FUNCTIONS = {
    'one_hot_encoding': gq.generate_one_hot_model,
}

MAPPING_FUNCTIONS = {
    'one_hot_encoding': ut.get_mapping_pairs_one_hot_encoding,
}

def main():
    print("=" * 80)
    print("Molecular Docking")
    print("=" * 80)
    print(f"Configuration:")
    print(f"  Seed: {SEED}")
    print(f"  Remove Hydrogens: {REMOVE_HYDROGENS}")
    print(f"  Encoding: {ENCODING}")
    print(f"  Lambda Weights: {LAMBDA_WEIGHTS}")
    print(f"  Number of Reads: {NUM_READS}")
    print(f"  Number of Sweeps: {NUM_SWEEPS}")
    print("-" * 80)

    model = ENCODING_FUNCTIONS[ENCODING](
        LIGAND_PATH1, POCKETGRID_CSV_PATH, 
        lambda_=LAMBDA_WEIGHTS, 
        removeHs=REMOVE_HYDROGENS
    )
    qubo, _ = model.to_qubo()
    quadratic_coeffs = dimod.BinaryQuadraticModel.from_qubo(qubo)

    # Count non-zero QUBO elements
    non_zero_elements = {k: v for k, v in qubo.items() if v != 0}

    # Separate linear and quadratic terms
    linear_terms = {k: v for k, v in non_zero_elements.items() if k[0] == k[1]}
    quadratic_terms = {k: v for k, v in non_zero_elements.items() if k[0] != k[1]}

    # Print counts in green
    print(f"  Non-zero linear terms:    {len(linear_terms)}")
    print(f"  Non-zero quadratic terms: {len(quadratic_terms)}")
    print(f"  Total non-zero terms:     {len(non_zero_elements)}")

    best_energy, best_sample = ut.solve_qubo(
        model, 
        seed=SEED,
        num_reads=NUM_READS,
        num_sweeps=NUM_SWEEPS
    )
    decoded_sample = model.decode_sample(best_sample, vartype='BINARY')
    mapping = MAPPING_FUNCTIONS[ENCODING](best_sample)

    print(mapping)

    if (
        decoded_sample.subh['penalty1'] != 0 or
        decoded_sample.subh['penalty2'] != 0
    ):
        print("Invalid Mapping")
        rmsd = float('nan')
        rmsd_wrt_exact_mapping = float('nan')
        rmsd_adjusted = float('nan')
        volume_overlap_ratio = float('nan')
    else:
        ligand_graph, pocketgrid_graph = gu.prepare_graphs(LIGAND_PATH1, POCKETGRID_CSV_PATH, removeHs=REMOVE_HYDROGENS)
        rmsd = ut.compute_rmsd(ligand_graph, pocketgrid_graph, mapping)
        rmsd_wrt_exact_mapping = ut.compute_rmsd_wrt_exact_mapping(ligand_graph, pocketgrid_graph, mapping)
        rmsd_adjusted = ut.compute_adjusted_rmsd(ligand_graph, pocketgrid_graph, mapping)
        volume_overlap_ratio = ut.compute_volume_overlap_ratio(ligand_graph, pocketgrid_graph, mapping)

    print(f"  Best Energy:              {best_energy}")
    print(f"  Geometric Energy:         {decoded_sample.subh['geometric_contribution']}")
    print(f"  vdW Energy:               {decoded_sample.subh['vdw_contribution']}")
    print(f"  Electrostatic Energy:     {decoded_sample.subh['electrostatic_contribution']}")
    if ADD_BENZENE_NODES:
        print(f"  Residual Energy:          {decoded_sample.subh['pistacking_contribution']}")
    print(f"  Penalty1:                 {decoded_sample.subh['penalty1']}")
    print(f"  Penalty2:                 {decoded_sample.subh['penalty2']}")
    print(f"  RMSD:                     {rmsd}")
    print(f"  RMSD wrt Exact Mapping:   {rmsd_wrt_exact_mapping}")
    print(f"  RMSD Adjusted:            {rmsd_adjusted}")
    print(f"  Volume Overlap Ratio:     {volume_overlap_ratio}")

    '''
    mapping = sorted(mapping, key=lambda x: x[0])
    ligand_coords = mu.get_coordinates_from_mol2_path(LIGAND_PATH1, removeHs=REMOVE_HYDROGENS)
    pocketgrid_nodes = list(pocketgrid_graph.nodes)
    new_coords = np.array([pocketgrid_graph.nodes[pocketgrid_nodes[i]]['coords'] for _, i in mapping])
    new_coords = ut.align_sets(ligand_coords, new_coords)
    
    ut.export_sdf(LIGAND_PATH1, 'old_ligand.sdf', REMOVE_HYDROGENS)
    ut.export_sdf_with_modified_coordinates(
        LIGAND_PATH1,
        'new_ligand.sdf',
        new_coords,
        REMOVE_HYDROGENS
    )
    '''

    print("=" * 80)
    return quadratic_coeffs

if __name__ == "__main__":
    main()