from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import re
vocab = ["pad", "start", "end", "sep",
                          "/", "\\", "@", "@@", "H",
                          "1", "2", "3", "4", "5", "6",
                          "#", "=", "-", "(", ")", "[", "]", "[*]", "([*])", "+"  
                          ]

def findCuttable(mol) -> list:
    """
    Determine cuttable single bonds that meets certain criteria. 
    
    Parameters:
    - mol (rdkit.Chem.Mol): The molecule for which to find cuttable bonds.

    Returns:
    - list of rdkit.Chem.Bond: List of cuttable single bonds.

    """
    cuttable_bonds = []

    for bond in mol.GetBonds():

        # if the bond is NOT part of a ring
        is_not_in_ring = not bond.IsInRing()

        # if neither end of the bond is a hydrogen atom
        no_hydrogen_at_ends = all(atom.GetAtomicNum() != 1 for atom in (bond.GetBeginAtom(), bond.GetEndAtom()))

        # if at least one end of the bond is attached to a ring
        attached_to_ring = any(atom.IsInRing() for atom in (bond.GetBeginAtom(), bond.GetEndAtom()))

        # if all conditions met
        if is_not_in_ring and no_hydrogen_at_ends and attached_to_ring:
            cuttable_bonds.append(bond)
    
    return cuttable_bonds

def visualizeFrag(mol,outfile):
    """
    Help visualise where the cuttable bonds are

    Parameters:
    - mol (rdkit.Chem.Mol)

    """

    cuttable_bonds = findCuttable(mol)
    Chem.Kekulize(mol)
    bond_indices = [bond.GetIdx() for bond in cuttable_bonds]

    img = Draw.MolToImage(mol, size=(600, 300), kekulize=False, wedgeBonds=True,
                          #highlightAtoms=highlight_atoms, 
                          highlightBonds=bond_indices)
    
    img.save(outfile)

def getRings(mol)->list:
    """
    Identify indices of atoms in a ring structure
    """
    ri = mol.GetRingInfo()
    return [list(ri.AtomRings()[i]) for i in range(ri.NumRings())]

def isRing_fragment(fragment, rings):
    """
    Check if an atom is in a ring structure
    """
    if fragment is None:
        return False
        if any(atom.GetIdx() in ring for ring in rings):
            return True
    return False

def getMatches(sequence, vocab) -> str:
    """
    For Fsmiles, the size information of each atom is concatenated 
    so as the other symbols such as wildcard substitution ([*]).
    This def helps identifying those symbols and add '_0' indicating no size information is found. 
    
    Returns:
    - Intermediate Fsmile presentation with size information
    """
    pattern = re.compile(r'\([^)]*\)|\[[^\]]*\]|[a-zA-Z]|[0-9]|[-+=#@]')
    result = pattern.findall(sequence)
    out = ""
    for s in result:
        if s in vocab:
            s += '_0'
        out+=s

    return out

def annotate(fragment_sequence, is_ring):
    # Add 1 to the start and end of the fragment if it's a ring
    if is_ring:
        results = fragment_sequence[0] + '1' + fragment_sequence[1:-1] + fragment_sequence[-1] + '1'
    else:
        results = fragment_sequence

    # Use the getMatches function to find matches in the vocabulary
    out = getMatches(results, vocab)

    # Substitute letters with _6 for non-ring fragments, and _0 for ring fragments
    output_sequence = re.sub(r'([a-zA-Z])', r'\1_0' if not is_ring else r'\1_6', out)
    return output_sequence
