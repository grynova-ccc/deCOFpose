import json
from rdkit import Chem
from rdkit.Chem import AllChem




def mol_from_periodic_smiles(periodic_smi):
    """
    Reconstruct an RDKit molecule from SMILES and apply "periodic" properties to specified bonds.
    
    Args:
    - smiles: SMILES string of the molecule.
    - periodic_bond_idxs: A list of bond indices that are periodic.
    
    Returns:
    - mol: An RDKit molecule object with "periodic" bonds specified.
    """
    if '|' not in periodic_smi:

        mol = Chem.MolFromSmiles(periodic_smi)
        mol.UpdatePropertyCache(strict=False)
        AllChem.Compute2DCoords(periodic_smi)
    else:
    
        # Separate the SMILES and indices
        smiles, indices_str = periodic_smi.split('|')
        periodic_bond_idxs = json.loads(indices_str)
        
        # Construct the molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        mol.UpdatePropertyCache(strict=False)
        AllChem.Compute2DCoords(mol)

        for bond in mol.GetBonds():
            if bond.GetIdx() in periodic_bond_idxs:
                bond.SetProp("isPeriodic", "true")
            else:
                bond.SetProp("isPeriodic", "false")

        
    
    return mol


def mol_to_periodic_smiles(mol,exteneded=True,isomeric=True):
    """
    Convert an RDKit molecule to SMILES and list of periodic bond indices.
    
    Args:
    - mol: An RDKit molecule object.
    - periodic_bond_idxs: A list of bond indices that are periodic.
    
    Returns:
    - smiles: SMILES string of the molecule.
    - periodic_bond_idxs: The same list of periodic bond indices for reconstruction.
    
    """

    periodic_bond_idxs=[]
    

    smiles = Chem.MolToSmiles(mol,isomericSmiles=isomeric,allHsExplicit=True)
    for bond in mol.GetBonds():
        if bond.HasProp("isPeriodic"):
            if bond.GetProp("isPeriodic")== "true":
                periodic_bond_idxs.append(bond.GetIdx())

    if exteneded== True:

        out=f'{smiles}|{periodic_bond_idxs}'
    else:
        out=smiles

    return out


