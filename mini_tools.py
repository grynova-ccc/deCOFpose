def get_largest_ring_size(rings):
    return max(len(ring) for ring in rings)

def get_bonds_in_smaller_rings(rings, largest_size):
    smaller_ring_bonds = set()
    for ring in rings:
        if len(ring) < largest_size: #-2
           
            smaller_ring_bonds.update(ring)
    return smaller_ring_bonds

def get_shared_bonds(rings, largest_size):
    largest_ring_bonds = set()
    for ring in rings:
        
        if len(ring) == largest_size:
            largest_ring_bonds.update(ring)
            
    return largest_ring_bonds

def highlight_shared_bonds(mol):
    # Get the ring systems of the molecule for bonds
    ri = mol.GetRingInfo()
    rings = ri.BondRings()
    #[print(len(i)) for i in rings]

    # Identify the largest ring size
    largest_size = get_largest_ring_size(rings)

    # Get bonds in the largest rings
    largest_ring_bonds = get_shared_bonds(rings, largest_size)

    # Get bonds in smaller rings
    smaller_ring_bonds = get_bonds_in_smaller_rings(rings, largest_size)

    # Exclude bonds that are in smaller rings from the set of largest ring bonds
    exclusive_largest_ring_bonds = largest_ring_bonds - smaller_ring_bonds
    

    return exclusive_largest_ring_bonds,largest_ring_bonds,smaller_ring_bonds
