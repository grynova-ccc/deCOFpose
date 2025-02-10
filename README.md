# deCOFposition

**deCOFpositon** is a Python-based toolkit for decomposing of Covalent Organic Frameworks. It identifies "suspicious" bonds in molecules (using user-defined rules), fragments the molecules accordingly, and offers additional functionalities such as ring detection, adjacency matrix creation, and periodic boundary condition (PBC) handling.



## Folder and File Structure

```plaintext
.
├── 447.cif
├── core.xyz
├── l_n_test/                  # Folder for output images 
├── periodic_smiels.py         # Script for generating/manipulating periodic SMILES
├── utils_pbc.py               # Helper functions for periodic boundary conditions
├── utils.py                   # General utilities (e.g., xyz2AC, AC2mol, etc.)
├── chem_rules.py              # Contains BaseRule and rule child classes (Rule1, Rule2, Rule3, Rule5, etc.)
├── chouse_bonds.py            # Script for bond selection (applies rules from chem_rules.py)
├── mini_tools.py              # Utility scripts for reading/writing data
└──deCOFpose_explanation.ipynb# Jupyter Notebook demonstrating code usage
```


## Implementing Custom Rules

You can easily add custom rules to the analysis pipeline. For example, to create a new rule just subclass the BaseRule Class


```python
class RuleX(BaseRule):
    def __init__(self):
        super().__init__(rule_name="RuleX")

    def _bond_matches_rule(self, rwmol, bond,
                           exclusive_largest_ring_bonds,
                           largest_ring_bonds,
                           smallest_ring_bonds):
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Example condition: double-bonded carbon-oxygen in a 7-membered ring
        if ((atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8) or
            (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6)):
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                if bond.IsInRingSize(7):
                    return True
        return False

    def _update_atom_props(self, rwmol, bond):
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Optionally, set custom flags: always mark O as used=1, C as used=2
        if a1.GetAtomicNum() == 8:
            a1.SetUnsignedProp("used", 1)
            a2.SetUnsignedProp("used", 2)
        else:
            a1.SetUnsignedProp("used", 2)
            a2.SetUnsignedProp("used", 1)
        a1.SetProp("rule", self.rule_name)
        a2.SetProp("rule", self.rule_name)
```
```