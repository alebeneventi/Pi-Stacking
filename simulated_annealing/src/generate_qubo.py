import numpy as np
import graph_utils as gu
from parameters import ALPHA_MULTIPLIER
from parameters import DIFFERENCE_TYPE
from parameters import ADD_BENZENE_NODES
from pyqubo import Array, SubH
from utils import compute_angle
from pyqubo import Binary


def generate_geometric_Hamiltonian_one_hot_encoding(
    X,
    ligand_graph,
    pocketgrid_graph,
    lambda_=None
):
    """
    Build Hamiltonian for the geometric part of the QUBO problem.
    
    Parameters
    ----------
    X : pyqubo.Array
        Binary variables for the QUBO model. Here, X[i, j] is 1 if the i-th node in the
        ligand graph is assigned to the j-th node in the pocketgrid graph.
    ligand_graph : networkx.Graph
        Ligand graph.
    pocketgrid_graph : networkx.Graph
        Pocketgrid graph.
    lambda_ : list
        Weights for the different terms in the Hamiltonian.
    
    Returns
    -------
    H_geometric : pyqubo.SubH
        Geometric Hamiltonian.
    max_length_difference : float
        Maximum length difference between ligand and pocketgrid edges.
    """
    list_ligand_nodes = list(ligand_graph.nodes())
    list_pocketgrid_nodes = list(pocketgrid_graph.nodes())
    
    # Compute the optimization term
    H_opt = 0
    max_length_difference = 0
    for u, v, ligand_edge_data in ligand_graph.edges(data=True):
        for u_prime, v_prime, pocketgrid_edge_data in pocketgrid_graph.edges(data=True):
            i = list_ligand_nodes.index(u)
            j = list_ligand_nodes.index(v)
            i_prime = list_pocketgrid_nodes.index(u_prime)
            j_prime = list_pocketgrid_nodes.index(v_prime)
            
            # Compute difference between the length of the ligand edge and the length
            # of the pocketgrid edge.
            if DIFFERENCE_TYPE == 'absolute':
                difference = abs(ligand_edge_data['length'] - pocketgrid_edge_data['length'])
            elif DIFFERENCE_TYPE == 'squared':
                difference = (ligand_edge_data['length'] - pocketgrid_edge_data['length'])**2
            else:
                raise ValueError('Invalid distance type')
            
            H_opt += difference * X[i, i_prime] * X[j, j_prime]
            H_opt += difference * X[i, j_prime] * X[j, i_prime]
            max_length_difference = max(max_length_difference, difference)
    
        for u_prime, pocketgrid_node_data in pocketgrid_graph.nodes(data=True):
            i = list_ligand_nodes.index(u)
            j = list_ligand_nodes.index(v)
            i_prime = list_pocketgrid_nodes.index(u_prime)
            
            if DIFFERENCE_TYPE == 'absolute':
                difference = ligand_edge_data['length']
            elif DIFFERENCE_TYPE == 'squared':
                difference = ligand_edge_data['length']**2
            else:
                raise ValueError('Invalid distance type')
            
            H_opt += difference * X[i, i_prime] * X[j, i_prime]
            max_length_difference = max(max_length_difference, difference)  
    
    
    if lambda_ is None:
        raise ValueError("Set LAMBDA_WEIGHTS in parameters.py to use the one-hot encoding model.")
    
    H_geometric = SubH(lambda_[0] * SubH(H_opt, label='geometric'), label='geometric_contribution')


    print(f"[DEBUG] Numero archi ligando: {ligand_graph.number_of_edges()}")
    print(f"[DEBUG] Numero archi pocketgrid: {pocketgrid_graph.number_of_edges()}")
    print(f"[DEBUG] Numero nodi pocketgrid: {pocketgrid_graph.number_of_nodes()}")

    # Termini attesi dal primo ciclo (con i_prime ≠ j_prime)
    expected_different = ligand_graph.number_of_edges() * pocketgrid_graph.number_of_edges() * 2

    # Termini attesi dal secondo ciclo (con i_prime = i_prime)
    expected_same = ligand_graph.number_of_edges() * pocketgrid_graph.number_of_nodes()

    print(f"[DEBUG] Termini attesi con indici diversi: {expected_different}")
    print(f"[DEBUG] Termini attesi con indici uguali: {expected_same}")


    return H_geometric, max_length_difference

def generate_penalty_terms_one_hot_encoding(
    X,
    max_length_difference,
    num_ligand_nodes,
    num_pocketgrid_nodes
):
    """
    Generate penalty terms for the one-hot encoding model.
    
    Parameters
    ----------
    X : pyqubo.Array
        Binary variables for the QUBO model (DW encoding).
    max_length_difference : float
        Maximum length difference between ligand and pocketgrid edges.
    num_ligand_nodes : int
        Number of nodes in the ligand graph.
    num_pocketgrid_nodes : int
        Number of nodes in the pocketgrid graph.
    
    Returns
    -------
    H_penalty : pyqubo.Add
        Penalty Hamiltonian.
    """
    # Generate penalty weight
    alpha = ALPHA_MULTIPLIER * max_length_difference

    # Generate penalty terms
    H_penalty1 = 0
    for i in range(num_ligand_nodes):
        H_penalty1 += (1 - sum(X[i, j] for j in range(num_pocketgrid_nodes)))**2

    H_penalty2 = 0
    for i in range(num_ligand_nodes):
        for j in range(i + 1, num_ligand_nodes):
            for i_prime in range(num_pocketgrid_nodes):
                H_penalty2 += X[i, i_prime] * X[j, i_prime]
    
    H_penalty = (
        SubH(alpha * SubH(H_penalty1, label='penalty1'), label='penalty1_contribution') +
        SubH(alpha * SubH(H_penalty2, label='penalty2'), label='penalty2_contribution')
    )

    return H_penalty

def generate_chemical_Hamiltonian_one_hot_encoding(
    X,
    ligand_graph,
    pocketgrid_graph,
    lambda_=None
):
    """
    Build Hamiltonian for the chemical part of the QUBO problem.
    
    Parameters
    ----------
    X : pyqubo.Array
        Binary variables for the QUBO model. Here, X[i, j] is 1 if the i-th node in the
        ligand graph is assigned to the j-th node in the pocketgrid graph.
    ligand_graph : networkx.Graph
        Ligand graph.
    pocketgrid_graph : networkx.Graph
        Pocketgrid graph.
    lambda_ : list
        Weights for the different terms in the Hamiltonian.
    
    Returns
    -------
    H_chemical : pyqubo.Add
        Chemical Hamiltonian.
    """
    list_ligand_nodes = list(ligand_graph.nodes())
    list_pocketgrid_nodes = list(pocketgrid_graph.nodes())
    
    # Van der Waals energy term
    H_vdw = 0
    for u, ligand_node_data in ligand_graph.nodes(data=True):
        for u_prime, pocketgrid_node_data in pocketgrid_graph.nodes(data=True):
            i = list_ligand_nodes.index(u)
            i_prime = list_pocketgrid_nodes.index(u_prime)
            lj_param_idx = ligand_node_data['lj_param_idx']
            H_vdw += pocketgrid_node_data['vdw_energy'][lj_param_idx] * X[i, i_prime]
    

    
    # Electrostatic term
    H_electrostatic = 0
    for u, ligand_node_data in ligand_graph.nodes(data=True):
        i = list_ligand_nodes.index(u)
        for u_prime, pocketgrid_node_data in pocketgrid_graph.nodes(data=True):
            i_prime = list_pocketgrid_nodes.index(u_prime)
            H_electrostatic += (1.4393*10**20)*ligand_node_data['charge'] * pocketgrid_node_data['electric_potential'] * X[i, i_prime]

    
    # Pi-stacking term
    if ADD_BENZENE_NODES:
        H_pistacking = 0 
        for u, ligand_node_data_u in ligand_graph.nodes(data=True):
            for v, ligand_node_data_v in ligand_graph.nodes(data=True):
                if u == v:
                    continue
                i = list_ligand_nodes.index(u)
                j = list_ligand_nodes.index(v)
                for u_prime, pocketgrid_node_data in pocketgrid_graph.nodes(data=True):
                    #if len(pocketgrid_node_data['pistacking']) > 0: 
                    for vertex in pocketgrid_node_data['pistacking']:      
                        if ligand_node_data_u['benzene_center'] == 'Y' and ligand_node_data_v['benzene_vertical'] == 'Y':
                            i_prime = list_pocketgrid_nodes.index(u_prime)
                            j_prime = list_pocketgrid_nodes.index(vertex[0])
                            H_pistacking += vertex[1] * X[i, i_prime] * X[j, j_prime]


        # Chemical Hamiltonian
        if lambda_ is None:
            raise ValueError("Set LAMBDA_WEIGHTS in parameters.py to use the one-hot encoding model.")
        
        H_chemical = (
            SubH(lambda_[1] * SubH(H_vdw, label='vdw'), label='vdw_contribution') +
            SubH(lambda_[2] * SubH(H_electrostatic, label='electrostatic'), label='electrostatic_contribution') +
            SubH(lambda_[3]  * SubH(H_pistacking, label='pistacking'), label='pistacking_contribution')
        )

    else:
        if lambda_ is None:
            raise ValueError("Set LAMBDA_WEIGHTS in parameters.py to use the one-hot encoding model.")
        
        H_chemical = (
            SubH(lambda_[1] * SubH(H_vdw, label='vdw'), label='vdw_contribution') +
            SubH(lambda_[2] * SubH(H_electrostatic, label='electrostatic'), label='electrostatic_contribution') 
        )


    
    return H_chemical

def build_one_hot_model(
    ligand_graph,
    pocketgrid_graph,
    lambda_=None
):
    """
    Build the QUBO model.
    
    Parameters
    ----------
    ligand_graph : networkx.Graph
        Ligand graph.
    pocketgrid_graph : networkx.Graph
        Pocketgrid graph.
    lambda_ : list
        Weights for the different terms in the Hamiltonian.
    
    Returns
    -------
    model : pyqubo.Model
        Problem model.
    """
    num_ligand_nodes = ligand_graph.number_of_nodes()
    num_pocketgrid_nodes = pocketgrid_graph.number_of_nodes()
    
    # Define the binary variables for the QUBO model. Here, X[i, j] is 1 if 
    # the i-th node in the ligand graph is assigned to the j-th node in the 
    # pocketgrid graph.
    X = Array.create(
        'X',
        shape=(num_ligand_nodes, num_pocketgrid_nodes),
        vartype='BINARY'   
    )
    
    # Generate the Hamiltonian for the geometric part of the QUBO problem
    H_geometric, max_length_difference = generate_geometric_Hamiltonian_one_hot_encoding(
        X,
        ligand_graph,
        pocketgrid_graph,
        lambda_=lambda_
    )
    
    H_chemical = generate_chemical_Hamiltonian_one_hot_encoding(
        X,
        ligand_graph,
        pocketgrid_graph,
        lambda_=lambda_
    )

    H_penalty = generate_penalty_terms_one_hot_encoding(
        X,
        max_length_difference,
        num_ligand_nodes,
        num_pocketgrid_nodes
    )

    
    # Build the QUBO model
    H = H_geometric + H_penalty + H_chemical 
    model = H.compile()
    return model

def generate_one_hot_model(
    ligand_path, 
    pocketgrid_path,
    lambda_,
    removeHs=True
):
    """
    Generate QUBO using one-hot encoding strategy.

    Parameters
    ----------
    ligand_path : str
        Path to the ligand file.
    pocketgrid_path : str
        Path to the pocketgrid file.
    lambda_ : list
        Weights for the different terms in the Hamiltonian.
    removeHs : bool
        Whether to remove hydrogens from the ligand graph.
    
    Returns
    -------
    model : pyqubo.Model
        Model with one-hot encoding.
    """
    ligand_graph, pocketgrid_graph = gu.prepare_graphs(ligand_path, pocketgrid_path, removeHs)
    model = build_one_hot_model(ligand_graph, pocketgrid_graph, lambda_=lambda_)
    return model
