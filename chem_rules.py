from abc import ABC, abstractmethod
from rdkit import Chem


def is_n_with_3c_in_6_ring(mol, atom_n):
    """
    Check if the provided nitrogen atom (atom_N) has three carbon neighbors,
    each of these carbons is inside a 6-membered ring, and all three carbon 
    atoms are part of different 6-membered rings.
    """
    if atom_n.GetSymbol() != 'N' or atom_n.GetDegree() != 3:
        return False

    ring_info = mol.GetRingInfo()

    # Helper function from your code
    def get_6_ring_membership(atom, ring_info):
        return [ring for ring in ring_info.AtomRings() if len(ring) == 6 and atom.GetIdx() in ring]

    carbon_rings = []
    for neighbor in atom_n.GetNeighbors():
        if neighbor.GetSymbol() == 'C':
            ring_list = get_6_ring_membership(neighbor, ring_info)
            if not ring_list:
                return False
            carbon_rings.append(ring_list)

    # Flatten the list
    flattened = [tuple(r) for sublist in carbon_rings for r in sublist]
    # Check that the rings are all distinct
    if len(flattened) != len(set(flattened)):
        return False

    return True


class BaseRule(ABC):
    """
    Abstract base class for 'rule' implementations.

    It does:
      1. Iterates over bonds in 'largest_ring_bonds'.
      2. For each bond, calls _bond_matches_rule(...) to see if it meets the child logic.
      3. If True, it calls _update_atom_props(...) and accumulates that bond's idx
         in the returned list.
    """

    def __init__(self, rule_name):
        self.rule_name = rule_name

    def apply(self, mol, exclusive_largest_ring_bonds, largest_ring_bonds, smallest_ring_bonds):
        """
        Orchestrates applying the rule to 'mol'.
        Returns (list_of_bond_ids_that_matched, updated_mol).
        """
        rwmol = Chem.RWMol(mol)
        bond_identified = []

        for b_idx in largest_ring_bonds:
            bond = rwmol.GetBondWithIdx(b_idx)
            matched, b_idx = self._bond_matches_rule(
                rwmol,
                bond,
                exclusive_largest_ring_bonds,
                largest_ring_bonds,
                smallest_ring_bonds
            )
            if matched:
                bond_identified.append(b_idx)
                bond = rwmol.GetBondWithIdx(b_idx)
                self._update_atom_props(rwmol, bond) 

        return bond_identified, rwmol.GetMol()

    @abstractmethod
    def _bond_matches_rule(self, rwmol, bond,
                           exclusive_largest_ring_bonds,
                           largest_ring_bonds,
                           smallest_ring_bonds):
        """
        Child classes must implement the actual check logic.
        Return True if the bond should be flagged, False otherwise and also a bond_id).
        """
        pass

    def _update_atom_props(self, rwmol, bond):
        """
        By default, mark the begin/end atoms:
          - begin_atom used=1, end_atom used=2
          - begin_atom rule=rule_name, end_atom rule=rule_name

        If the child class needs more complex logic (C=1, B=2, etc.),
        it can override this method.
        """
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        a1.SetUnsignedProp("used", 1)
        a2.SetUnsignedProp("used", 2)
        a1.SetUnsignedProp("rule", self.rule_name)
        a2.SetUnsignedProp("rule", self.rule_name)


class Rule1(BaseRule):
    """

    The core logic:
      - Check for a bond that is C-N (6,7).
      - N must have degree == 2, skipping some ring checks or O neighbors
      - Then find a neighbor of the C that is also a carbon, etc.
      - If it passes, we return True.
    """

    def __init__(self):
        super().__init__(rule_name=1)

    def _bond_matches_rule(self, rwmol, bond,
                           exclusive_largest_ring_bonds,
                           largest_ring_bonds,
                           smallest_ring_bonds):
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        # Are we dealing with a C-N pair?
        if ((atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 7) or
            (atom1.GetAtomicNum() == 7 and atom2.GetAtomicNum() == 6)):

            # Identify who is C, who is N
            if atom1.GetAtomicNum() == 7:
                n_atom, c_atom = atom1, atom2
            else:
                n_atom, c_atom = atom2, atom1

            # Condition: N has degree == 2
            if n_atom.GetDegree() != 2:
         
                return False, bond.GetIdx()
            
            # Additional check: skip if neighbor is O with degree=1
            for nbr in c_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                  
                    return False, bond.GetIdx()

            # Next, we attempt to find a neighbor of c_atom that is carbon
            # and the bond with c_atom is in exclusive_largest_ring_bonds,
            # not in ring 5 or 6, etc.  
            for nbr in c_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and c_atom.HasProp("used")==False and nbr.HasProp("used")==False:
                 
                    bond_idx_tmp = rwmol.GetBondBetweenAtoms(nbr.GetIdx(), c_atom.GetIdx()).GetIdx()

                    # Must be in the list of exclusive_largest_ring_bonds
                    if bond_idx_tmp in exclusive_largest_ring_bonds:
                       
                        bond_c_c = rwmol.GetBondBetweenAtoms(nbr.GetIdx(), c_atom.GetIdx())
                        if not bond_c_c.IsInRingSize(5) and not bond_c_c.IsInRingSize(6):
                            c_atom.SetUnsignedProp("used", 2)
                            nbr.SetUnsignedProp("used", 1)

                            return True,bond_c_c.GetIdx()

        return False, bond.GetIdx()

    def _update_atom_props(self, rwmol, bond):

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        a1.SetUnsignedProp("rule", self.rule_name)
        a2.SetUnsignedProp("rule", self.rule_name)


class Rule2(BaseRule):
    """
    The logic:
      - Must be a C-B bond
      - That bond is in exclusive_largest_ring_bonds, not ring size 5 or 6, etc.
    """

    def __init__(self):
        super().__init__(rule_name=2)

    def _bond_matches_rule(self, rwmol, bond,
                           exclusive_largest_ring_bonds,
                           largest_ring_bonds,
                           smallest_ring_bonds):
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        if ((atom1.GetAtomicNum() == 5 and atom2.GetAtomicNum() == 6) or
            (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 5)):

            bond_idx = bond.GetIdx()
            # Must be in exclusive_largest_ring_bonds & not ring 5 or 6
            if bond_idx in exclusive_largest_ring_bonds:
                if (not bond.IsInRingSize(5) and not bond.IsInRingSize(6)):
                    return True, bond_idx

        return False, bond.GetIdx()

    def _update_atom_props(self, rwmol, bond):

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        if a1.GetAtomicNum() == 5:
            b_atom, c_atom = a1, a2
        else:
            b_atom, c_atom = a2, a1

        b_atom.SetUnsignedProp("used", 1)
        b_atom.SetUnsignedProp("rule", self.rule_name)

        c_atom.SetUnsignedProp("used", 2)
        c_atom.SetUnsignedProp("rule", self.rule_name)


class Rule3(BaseRule):
    """
      - Searching for N(deg=3)-C, skipping if "C=O" or ring checks
      - Checking is_n_with_3c_in_6_ring(...) if needed
    """

    def __init__(self):
        super().__init__(rule_name=3)

    def _bond_matches_rule(self, rwmol, bond,
                           exclusive_largest_ring_bonds,
                           largest_ring_bonds,
                           smallest_ring_bonds):
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        mol = rwmol.GetMol()  # For ring checks

        # Condition: N with degree 3 connected to a C
        if (atom1.GetAtomicNum() == 7 and atom1.GetDegree() == 3 and atom2.GetAtomicNum() == 6):
            n_atom, c_atom = atom1, atom2
        elif (atom2.GetAtomicNum() == 7 and atom2.GetDegree() == 3 and atom1.GetAtomicNum() == 6):
            n_atom, c_atom = atom2, atom1
        else:
            return False, bond.GetIdx()

        # Possibly check if the bond is in exclusive_largest_ring_bonds
        bond_idx = bond.GetIdx()

        if bond_idx not in exclusive_largest_ring_bonds:

            return False, bond.GetIdx()
        if bond.IsInRingSize(5) or bond.IsInRingSize(6):

            return False, bond.GetIdx()

        # Check if N actually has 3 carbon neighbors in 6-membered rings
        if is_n_with_3c_in_6_ring(mol, n_atom):

            return False, bond.GetIdx()

        # “C=O check for c_atom if neighbor is O with degree=1
        for nbr in c_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                c_atom.SetUnsignedProp("C=O", 1)
                # If c_atom is "sick_3", we skip

                return False, bond.GetIdx()

        # Also check for “rule_exception” scenario 
        rule_exception = False

        for nbr in n_atom.GetNeighbors():         
            if nbr.GetAtomicNum() == 6 :
                for s_nbr in nbr.GetNeighbors():
                    if s_nbr.GetAtomicNum()==8 and s_nbr.GetDegree()==1:
                        nbr.SetUnsignedProp("C=O",1)

        for nbr in n_atom.GetNeighbors(): #it means CsingleN(h)singleC
            if nbr.GetAtomicNum() == 1 :
                for nbr in n_atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx()!=c_atom.GetIdx():
                        rule_exception=True
        
        if rule_exception == False and c_atom.HasProp("C=O") == False:
            
                return True, bond.GetIdx()
        if rule_exception == True and c_atom.HasProp("C=O") == False :
            if (c_atom.IsInRingSize(6) == True or c_atom.IsInRingSize(5) == True):
                return True, bond.GetIdx()

                    
        return False, bond.GetIdx()

    def _update_atom_props(self, rwmol, bond):

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        if a1.GetAtomicNum() == 7:
            n_atom, c_atom = a1, a2
        else:
            n_atom, c_atom = a2, a1

        n_atom.SetUnsignedProp("used", 1)
        n_atom.SetUnsignedProp("rule", self.rule_name)

        c_atom.SetUnsignedProp("used", 2)
        c_atom.SetUnsignedProp("rule", self.rule_name)


class Rule4(BaseRule):
    """
    The logic:
      - Looking for O(deg=2)-C
      - Checking for neighbors that might be 'exclude_4' if there's a carbonyl, etc.
      - Then ensuring the bond is in exclusive_largest_ring_bonds, not ring 5 or 6.
    """

    def __init__(self):
        super().__init__(rule_name=4)

    def _bond_matches_rule(self, rwmol, bond,
                           exclusive_largest_ring_bonds,
                           largest_ring_bonds,
                           smallest_ring_bonds):
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        # Check for O(deg=2)-C
        if ((atom1.GetAtomicNum() == 8 and atom1.GetDegree() == 2 and atom2.GetAtomicNum() == 6) or
            (atom2.GetAtomicNum() == 8 and atom2.GetDegree() == 2 and atom1.GetAtomicNum() == 6)):

            # Identify who is O, who is C
            if atom1.GetAtomicNum() == 8:
                o_atom, c_atom = atom1, atom2
            else:
                o_atom, c_atom = atom2, atom1

            # Check if there's an O neighbor that is a "exclude_4" scenario 
            #  if the neighbor is a carbon that has a single-bonded O)
            #checking for the "C-C=O" neigbour
            rule_fifield4 = False
            for nbr in o_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    for s_nbr in nbr.GetNeighbors():
                        if s_nbr.GetAtomicNum() == 8 and s_nbr.GetDegree() == 1:
                            nbr.SetUnsignedProp("exclude_4", 1)
                            rule_fifield4 = True

            # Also check if c_atom has "sick_5" neighbors
            for nbr in c_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                    c_atom.SetUnsignedProp("exclude_4", 1)

            # Must be in exclusive_largest_ring_bonds, not ring size 5 or 6
            bond_idx = bond.GetIdx()
            if not rule_fifield4:
                return False , bond.GetIdx()
            if bond_idx not in exclusive_largest_ring_bonds:
                return False, bond.GetIdx()
            if bond.IsInRingSize(5) or bond.IsInRingSize(6):
                return False, bond.GetIdx()

            # Also skip if c_atom has "sick_5" (like your code does)
            if c_atom.HasProp("exclude_4"):
                return False, bond.GetIdx()

            return True, bond.GetIdx()

        return False, bond.GetIdx()

    def _update_atom_props(self, rwmol, bond):

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        if a1.GetAtomicNum() == 8:
            o_atom, c_atom = a1, a2
        else:
            o_atom, c_atom = a2, a1

        o_atom.SetUnsignedProp("used", 1)
        o_atom.SetUnsignedProp("rule", self.rule_name)

        c_atom.SetUnsignedProp("used", 2)
        c_atom.SetUnsignedProp("rule", self.rule_name)