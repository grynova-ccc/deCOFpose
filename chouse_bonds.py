
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import RWMol
from rdkit.Chem import AllChem
from utils import xyz2AC, AC2mol

def get_all_paths_of_length(mol, start, end, length):
    """Get all paths of a specific length between two atoms using BFS."""
    queue = [[start]]
    paths = []
    
    while queue:
        path = queue.pop(0)
        node = path[-1]

        # Reset visited nodes for the new path
        visited = {atom: False for atom in range(mol.GetNumAtoms())}
        for atom in path:
            visited[atom] = True
        
        if node == end and len(path) == length:
            paths.append(path)
        else:
            for neighbor in [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(node).GetNeighbors()]:
                if not visited[neighbor] and len(path) < length:
                    new_path = list(path)
                    new_path.append(neighbor)
                    queue.append(new_path)
    return paths


def create_graph_from_bonds(mol, bond_ids):
    G = nx.Graph()
    shortest_path_all=float("inf")

    for i, bond_id1 in enumerate(bond_ids):
        for bond_id2 in bond_ids[i+1:]:
            bond1 = mol.GetBondWithIdx(bond_id1)
            bond2 = mol.GetBondWithIdx(bond_id2)

            atom1_idx = bond1.GetOtherAtomIdx(bond1.GetBeginAtomIdx())
            atom2_idx = bond2.GetOtherAtomIdx(bond2.GetBeginAtomIdx())

            shortest_path = Chem.rdmolops.GetShortestPath(mol, atom1_idx, atom2_idx)
            shortest_path_length = len(shortest_path)
            if shortest_path_all>len(shortest_path):
                shortest_path_all=len(shortest_path)

            if shortest_path_length <shortest_path_all+4:
                print(shortest_path_length)

                # Retrieve all paths of that length using the BFS method
                all_paths_of_length = get_all_paths_of_length(mol, atom1_idx, atom2_idx, shortest_path_length)
                
                # If there's only one path, use it. Otherwise, sort to find the deterministic path
                if len(all_paths_of_length) == 1:
                    selected_path = shortest_path
                else:
                    print("yay")
                    all_paths_of_length = sorted(all_paths_of_length, key=lambda path: ''.join([mol.GetAtomWithIdx(atom_id).GetSymbol()[-1] for atom_id in path]))
                    selected_path = all_paths_of_length[0]


                # Remove atoms that are directly part of the bonds
                path = [atom for atom in selected_path if atom != bond1.GetBeginAtomIdx() and atom != bond1.GetEndAtomIdx() and atom != bond2.GetBeginAtomIdx() and atom != bond2.GetEndAtomIdx()]
                
                # Convert this path into an edge label
                edge_label = ''.join([mol.GetAtomWithIdx(atom_id).GetSymbol()[-1] for atom_id in path])
                
                node_label1 = bond1.GetBeginAtom().GetSymbol() + '-' + bond1.GetEndAtom().GetSymbol()
                node_label2 = bond2.GetBeginAtom().GetSymbol() + '-' + bond2.GetEndAtom().GetSymbol()
                
                # Add nodes and edge between them to the graph
                if not G.has_node(bond_id1):
                    G.add_node(bond_id1, label=node_label1)
                if not G.has_node(bond_id2):
                    G.add_node(bond_id2, label=node_label2)
                
                G.add_edge(bond_id1, bond_id2, label=edge_label)
    return G

def create_adjacency_matrix(G):
    nodes = list(G.nodes())
    n = len(nodes)
    matrix = np.empty((n, n), dtype=object)

    for i in range(n):
        for j in range(n):
            if i != j and G.has_edge(nodes[i], nodes[j]):
                matrix[i, j] = G[nodes[i]][nodes[j]]['label']
            else:
                matrix[i, j] = ''

    return matrix, nodes  # Return both matrix and the order of nodes


def map_indices_to_bond_ids(indices, nodes_order):
    # If it's a list of pairs, handle recursively
    if isinstance(indices[0], (list, tuple)):
        return [map_indices_to_bond_ids(pair, nodes_order) for pair in indices]
    
    i, j = indices
    return nodes_order[i], nodes_order[j]


def get_edge_label(matrix, i, j):
    # Get edge label from adjacency matrix directly.
    return matrix[i, j] if i < j else matrix[j, i]



def find_smallest_labels(matrix):
    n = matrix.shape[0]
    
    # Get a list of all labels and their inverses that have at least a length of 1
    labels = [(get_edge_label(matrix, i, j), get_edge_label(matrix, i, j)[::-1])
              for i in range(n) for j in range(n) if i < j and len(get_edge_label(matrix, i, j)) >= 1]
    
    # Flatten the list of tuples to a list of strings
    flat_labels = [label for sublist in labels for label in sublist]
    
    # If there are no such labels, return an empty list
    if not flat_labels:
        return []

    # Identify the smallest label based on its length
    smallest_label = min(flat_labels, key=len)
    print(smallest_label)
    

    # Extract bonds with the smallest label or its inverse
    bonds_with_smallest_label = [(i, j) for i in range(n) for j in range(n) 
                                 if get_edge_label(matrix, i, j) == smallest_label or get_edge_label(matrix, i, j)[::-1] == smallest_label]

    return bonds_with_smallest_label



def find_symmetrical_triangles_or_all(matrix, nodes_order):
    # First, find all the smallest labels
    smallest_label_bond_pairs = find_smallest_labels(matrix)
    
    
    symmetrical_triangles = []

    for i in range(len(smallest_label_bond_pairs)):
        for j in range(i+1, len(smallest_label_bond_pairs)):
            bond1 = smallest_label_bond_pairs[i]
            bond2 = smallest_label_bond_pairs[j]
            common_nodes = set(bond1).intersection(bond2)
            
            if len(common_nodes) == 1:
                common_node = list(common_nodes)[0]
                remaining_nodes = set(bond1 + bond2) - common_nodes
                
                # Now check if the third bond exists
                for node in remaining_nodes:
                    bond3 = tuple(sorted([common_node, node]))
                    if bond3 in smallest_label_bond_pairs:
                        triangle = [bond1, bond2, bond3]
                        triangle_sorted = sorted(triangle, key=lambda x: (x[0], x[1]))
                        if triangle_sorted not in symmetrical_triangles:
                            symmetrical_triangles.append(triangle_sorted)
    
                
    if symmetrical_triangles:
        # Flatten the triangles into a list of pairs
        flattened_triangles = [pair for triangle in symmetrical_triangles for pair in triangle]
        # Convert matrix indices to bond IDs before returning
        return [map_indices_to_bond_ids(pair, nodes_order) for pair in flattened_triangles]

    # If no symmetrical triangles found, return all smallest labels mapped to bond IDs
    return [map_indices_to_bond_ids(pair, nodes_order) for pair in smallest_label_bond_pairs]


def assign_molecule_type(mol):
    """
    Assigns a new property 'type' to the molecule based on the values of the 'used' property of its atoms.
    Also calculates the number of atoms with the 'used' property and sets it as 'number_connections' mol property.
    """
    used_values = [atom.GetProp("used") for atom in mol.GetAtoms() if atom.HasProp("used") and atom.HasProp("broke_id")]
    rule_values = [atom.GetProp("rule") for atom in mol.GetAtoms() if atom.HasProp("rule") and atom.HasProp("broke_id")]
    
    # Count the number of atoms with the 'used' property
    mol.SetProp("number_connections", str(len(used_values)))
    
    # Determine the molecule type
    unique_used_values = set(used_values)
    if unique_used_values == {"1"}:
        mol_type = 1
    elif unique_used_values == {"2"}:
        mol_type = 2
    elif unique_used_values == {"1", "2"}:
        mol_type = 3
    else:
        mol_type = 0

    mol.SetProp("type", str(mol_type))
    mol.SetProp("rule", rule_values[0])

    return mol


def get_fragments_with_unique_ids(frags):
    # Convert fragments to canonical SMILES
    smiles_frags = [(Chem.MolToSmiles(frag,isomericSmiles=False), idx) for idx, frag in enumerate(frags)]

    seen = {}
    smiles_to_idx = {}

    # Count the occurrences of each fragment and keep track of their original indices
    for smiles, idx in smiles_frags:
        if smiles in seen:
            seen[smiles] += 1
        else:
            seen[smiles] = 1
            smiles_to_idx[smiles] = idx

    unique_frags_with_ids = []

    # Create the list of unique fragments with their counts
    for smiles in seen:  # iterate over the unique smiles
        idx = smiles_to_idx[smiles]
        mol = frags[idx]
        mol.SetProp("amount", str(seen[smiles]))
        unique_frags_with_ids.append(assign_molecule_type(mol))

    return unique_frags_with_ids



def guilatine(mol, list_bonds):
    #?
    with Chem.RWMol(mol) as rwmol:
        list_additional_h = []
        broke_id = 0
        
        for b_idx in list_bonds:
            bond = rwmol.GetBondWithIdx(b_idx)
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()

            rwmol.RemoveBond(begin_atom.GetIdx(), end_atom.GetIdx())
            
            first_h = rwmol.AddAtom(Chem.Atom('H'))
            second_h = rwmol.AddAtom(Chem.Atom('H'))
            
            list_additional_h.append(first_h)
            list_additional_h.append(second_h)
            broke_id += 1
            
            begin_atom.SetUnsignedProp("broke_id", broke_id)
            end_atom.SetUnsignedProp("broke_id", broke_id)
        
        count_for_atoms = 0
        for atom in rwmol.GetAtoms():
            if atom.HasProp("broke_id"):
                rwmol.AddBond(atom.GetIdx(), list_additional_h[count_for_atoms], Chem.rdchem.BondType.SINGLE)
                count_for_atoms += 1

        # Set all bond orders to single and remove aromaticity
        for bond in rwmol.GetBonds():
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)
            bond.SetIsAromatic(False)
            
        for atom in rwmol.GetAtoms():
            atom.SetFormalCharge(0)  # Ensure all atoms have a charge of 0
            atom.SetIsAromatic(False)
        
    return rwmol.GetMol()



def frag_func(out_frag):
    
    fragy_list=[]
    for fragy in out_frag:
        print("yps")
        for atom in fragy.GetAtoms():
            if atom.GetAtomicNum() == 1:  # Hydrogen has atomic number 1
                atom.SetIsotope(2)  # Set isotope to 2 for deuterium
        #fragy.RemoveAllConformers()
        Chem.SanitizeMol(fragy)
        
        #status = AllChem.EmbedMolecule(fragy, AllChem.ETKDGv2())
        #status =AllChem.EmbedMolecule(fragy, AllChem.DistanceGeometry())
        

        status = AllChem.EmbedMolecule(fragy, useRandomCoords=True, useBasicKnowledge=False)
        if status != -1:
            AllChem.UFFOptimizeMolecule(fragy)

            #status=AllChem.EmbedMolecule(fragy, useBasicKnowledge=False, useExpTorsionAnglePrefs=False, useMacrocycleTorsion=False, useMacrocycle14config=False)
            #AllChem.UFFOptimizeMolecule(fragy)


            #AllChem.MMFFOptimizeMolecule(fragy)

            #if status == -1:
                #print("Embedding failed!")
        else:
            conf = fragy.GetConformer(0)
        coords_3d = []
        conf = fragy.GetConformer(0)  # Get the generated 3D conformer
        
        for atom_idx in range(fragy.GetNumAtoms()):
            pos = conf.GetAtomPosition(atom_idx)
            coords_3d.append((fragy.GetAtomWithIdx(atom_idx).GetSymbol(), pos.x, pos.y, pos.z))


        xyz_str = "{}\n\n".format(len(coords_3d))  # Number of atoms at the start


        for atom_info in coords_3d:
            xyz_str += "{}  {:.4f}  {:.4f}  {:.4f}\n".format(atom_info[0], atom_info[1], atom_info[2], atom_info[3])

        raw_mol = Chem.MolFromXYZBlock(xyz_str)
        conf = raw_mol.GetConformer()

        # Get atomic positions
        atomic_positions = [conf.GetAtomPosition(i) for i in range(raw_mol.GetNumAtoms())]

        # Convert positions to tuples
        atomic_positions = [(pos.x, pos.y, pos.z) for pos in atomic_positions]

        # Get atomic numbers
        atomic_numbers = [atom.GetAtomicNum() for atom in raw_mol.GetAtoms()]
        adj, xyzs= xyz2AC(atomic_numbers,atomic_positions, 0, )
        rw_mol2 = RWMol(xyzs)
        #atoms_new=[i.GetAtomicNum()for i in rw_mol2.GetAtoms()]

        fragmented = AC2mol(rw_mol2, adj , atomic_numbers, 0,
                        allow_charged_fragments=False,
                        use_graph=True,
                        use_atom_maps=False)
        fragmented=fragmented[0]



        #conn_mol = Chem.Mol(raw_mol)
        #rdDetermineBonds.DetermineBonds(conn_mol, charge=0)
            # Transfer "used" and "broke_id" properties
        for src_atom, tgt_atom in zip(fragy.GetAtoms(), fragmented.GetAtoms()):
            if src_atom.HasProp("used"):
                tgt_atom.SetProp("used", src_atom.GetProp("used"))
            if src_atom.HasProp("broke_id"):
                tgt_atom.SetProp("broke_id", src_atom.GetProp("broke_id"))
            if src_atom.HasProp("rule"):
                tgt_atom.SetProp("rule", src_atom.GetProp("rule"))
        fragmented.RemoveAllConformers()
        Chem.SanitizeMol(fragmented)





        fragy_list.append(fragmented)
    return fragy_list


def combine(list_mol):
    

    for m_atm in list_mol[0].GetAtoms():
        m_atm.SetUnsignedProp("mol_e",1)

    for m_atm in list_mol[1].GetAtoms():
        m_atm.SetUnsignedProp("mol_e",2)



    comb=Chem.rdmolops.CombineMols(list_mol[0],list_mol[1])
    cleaned_1=False
    cleaned_2=False

    for i in comb.GetAtoms():
        if i.HasProp("used")==True:
            if i.GetUnsignedProp("used")==1 and i.HasProp("broke_id"):
                first_atom=i
                mol_t=i.GetUnsignedProp("mol_e")
                print("one")
                break
    for i in comb.GetAtoms():
        if i.HasProp("used")==True:
            if i.GetUnsignedProp("used")==2 and i.GetUnsignedProp("mol_e")!=mol_t and i.HasProp("broke_id"):
                second_atom=i
                
                print("one2")
                break

    with Chem.RWMol(comb) as rwmol:

        rwmol.AddBond(first_atom.GetIdx(),second_atom.GetIdx(),Chem.rdchem.BondType.SINGLE)

        for nbr in first_atom.GetNeighbors():
            if nbr.GetAtomicNum()==1 and cleaned_1==False:
                rwmol.RemoveAtom(nbr.GetIdx())
                cleaned_1=True

        for nbr in second_atom.GetNeighbors():
            if nbr.GetAtomicNum()==1 and cleaned_2==False:
                rwmol.RemoveAtom(nbr.GetIdx())
            
                cleaned_2=True
    return rwmol








