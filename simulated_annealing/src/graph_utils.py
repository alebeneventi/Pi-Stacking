import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import mol2_utils as mu
import utils as ut
from parameters import LJ_PARAMS_FILE_PATH, ADD_BENZENE_NODES
from rdkit import Chem

def is_graph_edge_rotatable(G, edge):
    """
    Check if an edge in a graph is rotatable from a geometric point of view.
    
    Parameters
    ----------
    G : networkx.Graph
        Graph to check.
    edge : tuple
        Edge to check.
    
    Returns
    -------
    bool
        True if the edge is rotatable, False otherwise.
    """
    # If removing the edge causes the graph to split into two disconnected components, the edge is rotatable.
    frame_graph = G.copy()
    frame_graph.remove_edge(edge[0], edge[1])
    
    if not nx.is_connected(frame_graph):
        disconnected_subgraph_1 = nx.node_connected_component(frame_graph, edge[0])
        disconnected_subgraph_2 = nx.node_connected_component(frame_graph, edge[1])
        if len(disconnected_subgraph_1) != 1 and len(disconnected_subgraph_2) != 1:
            return True
    return False
    

def add_benzene_nodes(G, mol2_path, removeHs=False, n_connections=2):
    """
    Identify benzene rings and add:
      - A RingCenter node at the centroid of each ring.
      - Edges (n_connections) from the centroid to ring atoms.
      - A second RingCenter node 1Å above the ring plane, connected only to the centroid.

    Parameters
    ----------
    G : networkx.Graph
        Prebuilt ligand graph.
    mol2_path : str
        Path to ligand .mol2 file.
    removeHs : bool
        Whether to ignore hydrogens when detecting rings.
    n_connections : int
        Number of ring atoms to anchor the centroid to (default 2).
    """
    mol = Chem.MolFromMol2File(mol2_path, removeHs=removeHs)
    if mol is None:
        raise ValueError(f"Cannot load MOL2 from {mol2_path}")
    rings = mol.GetRingInfo().AtomRings()

    existing_ids = [n for n in G.nodes() if isinstance(n, int)]
    next_id = max(existing_ids, default=-1) + 1

    for ring in rings:
        if len(ring) != 6 or not all(
            mol.GetAtomWithIdx(i).GetSymbol() == 'C' and mol.GetAtomWithIdx(i).GetIsAromatic()
            for i in ring
        ):
            continue

        # Compute centroid
        coords = np.array([G.nodes[i]['coords'] for i in ring])
        centroid = coords.mean(axis=0)

        # Add primary center
        center_id = next_id
        G.add_node(
            center_id,
            atom_type='GhostAtom',
            mmff_atom_type=-1,
            coords=centroid,
            lj_param_idx=-1,
            hbond_coloring='N',
            charge=0.0,
            hydrophobic_coloring='NH',
            ring_members=tuple(ring),
            benzene_center = 'Y',
            benzene_vertical = 'N'
        )

        # Anchor centroid to ring atoms
        step = len(ring) // n_connections
        for k in range(n_connections):
            atom_idx = ring[k * step]
            length = np.linalg.norm(centroid - G.nodes[atom_idx]['coords'])
            G.add_edge(
                center_id,
                atom_idx,
                bond_type='ring_center',
                length=length,
                rotatable=False
            )

        # Compute ring normal
        v1 = coords[2] - coords[0]
        v2 = coords[1] - coords[0]
        normal = np.cross(v1, v2)
        normal /= np.linalg.norm(normal)
        normal = normal*1.5
        # Add secondary node 1.5Å above plane
        second_pos = centroid + normal
        second_id = center_id + 1
        G.add_node(
            second_id,
            atom_type='GhostAtom',
            mmff_atom_type=-1,
            coords=second_pos,
            lj_param_idx=-1,
            hbond_coloring='N',
            charge=0.0,
            hydrophobic_coloring='NH',
            ring_members=tuple(),
            benzene_center = 'N',
            benzene_vertical = 'Y'
        )

                # Connect secondary node to primary center
        G.add_edge(
            center_id,
            second_id,
            bond_type='ring_center_extension',
            length=1.5,
            rotatable=False
        )

        # Also anchor secondary node to two ring carbon atoms
        # Select two opposite carbons: index 0 and index 3 of the ring
        anchor_indices = [ring[0], ring[1]]
        for atom_idx in anchor_indices:
            length = np.linalg.norm(second_pos - G.nodes[atom_idx]['coords'])
            G.add_edge(
                second_id,
                atom_idx,
                bond_type='ring_center_anchor',
                length=length,
                rotatable=False
            )

        # Update next available ID
        next_id = second_id + 1


def ligand_to_graph(mol2_path, removeHs=False):
    """
    Convert a ligand in a mol2 file to a graph.

    Parameters
    ----------
    mol2_path : str
        Path to the mol2 file.
    removeHs : bool, optional
        Remove hydrogens from the molecule. The default is False.
        
    Returns
    -------
    networkx.Graph
        Graph of the ligand.
    """
    elementary_charge = 1.60217662e-19
    atoms_df = mu.get_dataframe_from_mol2_path(mol2_path, LJ_PARAMS_FILE_PATH)
    bonds_df = mu.get_bonds_dataframe_from_mol2_path(mol2_path)
    
    # Create graph
    G = nx.Graph()
    
    # Add nodes
    for row in atoms_df.itertuples():
        if removeHs and row.atom_type == 'H':
            continue
        G.add_node(
            row.atom_ID, 
            atom_type=row.atom_type, 
            mmff_atom_type=row.mmff_atom_type,
            coords=np.array([row.x, row.y, row.z]),
            lj_param_idx=row.lj_param_idx,
            charge=row.mmff_partial_charge * elementary_charge,
            benzene_center = 'N',
            benzene_vertical = 'N'
        )
    
    # Add edges
    for row in bonds_df.itertuples():
        if removeHs and (
            atoms_df.loc[row.atom1_ID].atom_type == 'H' or
            atoms_df.loc[row.atom2_ID].atom_type == 'H'
        ):
            continue
        G.add_edge(
            row.atom1_ID, 
            row.atom2_ID, 
            bond_type=row.bond_type,
            length=np.linalg.norm(
                G.nodes[row.atom1_ID]['coords'] - G.nodes[row.atom2_ID]['coords']
            )
        )
    
    # Add H distances to the graph nodes
    H_distances = {node: [] for node in list(G.nodes())}
    for row in bonds_df.itertuples():
        if (
            atoms_df.loc[row.atom1_ID].atom_type != 'H' and
            atoms_df.loc[row.atom2_ID].atom_type == 'H'
        ):
            H_distances[row.atom1_ID].append(
                np.linalg.norm(
                    np.array(atoms_df.loc[row.atom1_ID][['x', 'y', 'z']]) - np.array(atoms_df.loc[row.atom2_ID][['x', 'y', 'z']])
                )
            )
        elif (
            atoms_df.loc[row.atom1_ID].atom_type == 'H' and
            atoms_df.loc[row.atom2_ID].atom_type != 'H'
        ):
            H_distances[row.atom2_ID].append(
                np.linalg.norm(
                    np.array(atoms_df.loc[row.atom1_ID][['x', 'y', 'z']]) - np.array(atoms_df.loc[row.atom2_ID][['x', 'y', 'z']])
                )
            )
    nx.set_node_attributes(G, values=H_distances, name='H_distances')
    
    # Find rotatable bonds
    rotatable_bonds = {
        (source, target): False
        for source, target in G.edges()
    }
    for source, target, data in G.edges(data=True):
        if (
            data['bond_type'] not in ['1.5', '2.0'] and
            is_graph_edge_rotatable(G, (source, target))
        ):
            rotatable_bonds[(source, target)] = True
    nx.set_edge_attributes(G, values=rotatable_bonds, name='rotatable')

    if ADD_BENZENE_NODES:
        add_benzene_nodes(G, mol2_path, removeHs)
    
    return G

def get_edges_for_angle_preservation(G):
    """
    Get the edges to be added to a graph to preserve angles.
    
    Parameters
    ----------
    G : networkx.Graph
        Graph to add edges to.
    
    Returns
    -------
    list[tuple]
        List with the edges to add to the graph.
    """
    # Given a node, add an edge between the node and the neighbors of its neighbors.
    to_add = []
    for node in G.nodes():
        for neighbor in G.neighbors(node):
            for neighbor_neighbor in G.neighbors(neighbor):
                if neighbor_neighbor != node:
                    to_add.append((node, neighbor_neighbor))
    return to_add

def get_edges_for_non_rotatable_bonds(G):
    """
    Get the edges to be added to a graph to make bonds non-rotatable when
    they are non-rotatable for chemical reasons.
    
    Parameters
    ----------
    G : networkx.Graph
        Graph to add edges to.
    
    Returns
    -------
    list[tuple]
        List with the edges to add to the graph.
    """
    to_add = []
    for u, v, data in G.edges(data=True):
        if data['bond_type'] in ['1.5', '2.0'] and is_graph_edge_rotatable(G, (u, v)):
            u_neighbors = set(G.neighbors(u))
            v_neighbors = set(G.neighbors(v))
            
            for u_neighbor in u_neighbors:
                for v_neighbor in v_neighbors:
                    if u_neighbor != v_neighbor and not G.has_edge(u_neighbor, v_neighbor):
                        if (u_neighbor, v_neighbor) not in to_add and (v_neighbor, u_neighbor) not in to_add:
                            to_add.append((u_neighbor, v_neighbor))
    return to_add
    
def get_augmented_graph(G):
    """
    Get an augmented version of a graph with added edges to preserve angles.
    
    Parameters
    ----------
    G : networkx.Graph
        Graph to augment.
    
    Returns
    -------
    networkx.Graph
        Augmented graph.
    """
    # Create a copy of the graph
    G_augmented = G.copy()
    
    # Get the edges to add to the graph
    edges_to_add = get_edges_for_angle_preservation(G)
    edges_to_add.extend(get_edges_for_non_rotatable_bonds(G))
    
    # Add the edges to the graph
    for edge in edges_to_add:
        if not G_augmented.has_edge(edge[0], edge[1]):
            G_augmented.add_edge(
                edge[0], 
                edge[1], 
                bond_type='A',
                length=np.linalg.norm(
                    G_augmented.nodes[edge[0]]['coords'] - G_augmented.nodes[edge[1]]['coords']
                ),
                rotatable=False
            )
    
    return G_augmented

def pocketgrid_to_graph(pocketgrid_csv):
    """
    Converts a pocketgrid CSV file into a NetworkX graph.
    
    Parameters
    ----------
    pocketgrid_csv : str
        Path to the pocket grid CSV file.
    
    Returns
    -------
    networkx.Graph
        A graph representation of the pocket grid with nodes and edges.
    """
    pocketgrid_df = pd.read_csv(pocketgrid_csv)
    
    # Create graph
    G = nx.Graph()
    
    # Add nodes
    for row in pocketgrid_df.itertuples():
        G.add_node(
            row.atom_ID, 
            coords=np.array([row.x, row.y, row.z]),
            vdw_energy=eval(row.vdw_energy),
            electric_potential=row.electric_potential,
            pistacking = eval(row.pistacking) 
        )
    
    #Add edges
    nodes = list(G.nodes)
    for i in range(len(nodes)):
       for j in range(i + 1, len(nodes)):
           dist = np.linalg.norm(G.nodes[nodes[i]]['coords'] - G.nodes[nodes[j]]['coords'])
           G.add_edge(nodes[i], nodes[j], length=dist)
    
    return G

def prepare_graphs(
    ligand_path,
    pocketgrid_path, 
    removeHs=True
):
    """
    Prepare ligand and pocket grid graphs from input files.

    Parameters
    ----------
    ligand_path : str
        Path to the ligand file (in Mol2 format).
    pocketgrid_path : str
        Path to the pocket grid file (in CSV format).
    removeHs : bool, optional
        Whether to remove hydrogens from the ligand molecule. Default is True.
    
    Returns
    -------
    tuple
        A tuple containing the ligand graph and the pocket grid graph.
    """
    ligand_graph = ligand_to_graph(ligand_path, removeHs=removeHs)
    ligand_graph = get_augmented_graph(ligand_graph)
    pocketgrid_graph = pocketgrid_to_graph(pocketgrid_path)
    
    return ligand_graph, pocketgrid_graph

def parse_pdb_points(pdb_block):
    """
    Parses ATOM lines from PDB-like block and returns coordinates and labels.
    """
    coords = []
    labels = []
    for ln in pdb_block.strip().splitlines():
        m = re.match(
            r'^ATOM\s+\d+\s+(\S+)\s+\S+\s+\S+\s+\d+\s+'  # atom name
            r'([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)',     # x y z
            ln)
        if not m:
            continue
        atom_name, xs, ys, zs = m.groups()
        coords.append([float(xs), float(ys), float(zs)])
        labels.append(atom_name)
    return np.array(coords), labels







