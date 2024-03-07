from rdkit import Chem
from utils.helper import visualizeFrag,findCuttable,getRings,annotate
import os,sys

def fsmlsConverter(mol): 
    """
    Convert a molecule into FSMILES representation. 
    This function finds cuttable bonds, extracts fragments, and generates FSMILES representation

    Parameters:
    - mol (rdkit.Chem.Mol): The molecule to convert into FSMILES.
    
    """
    
    cuttable_bonds = findCuttable(mol)

    neighbors = set([bond.GetBeginAtom().GetIdx() for bond in cuttable_bonds] +
                    [bond.GetEndAtom().GetIdx() for bond in cuttable_bonds])

    double_bonds = set([bond.GetBeginAtom().GetIdx() for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2.0] +
                       [bond.GetEndAtom().GetIdx() for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2.0])

    # Create a list to track visited atoms
    visited_atoms = set()

    def dfsearch(atom, fragment_atoms):
        """
        Traverse the molecule and annotate fragments

        Parameters:
        - atom (rdkit.Chem.Atom): The current atom in the traversal.
        - fragment_atoms (list): List to store the atoms and their annotations in the fragment.
        
        Returns:
        - None: Modifies the fragment_atoms list in place.
        """
        visited_atoms.add(atom.GetIdx())
        for neighbor_atom in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor_atom.GetIdx())
            if bond and bond.GetIdx() not in [b.GetIdx() for b in cuttable_bonds]:
                if neighbor_atom.GetIdx() not in visited_atoms:
                    fragment_atoms.append((neighbor_atom.GetSymbol(), '([*])' if neighbor_atom.GetIdx() in neighbors else '', '=' if neighbor_atom.GetIdx() in double_bonds else ''))
                    dfsearch(neighbor_atom, fragment_atoms)

    # Iterate through atoms and print sequences for fragments
    fsmiles = ""
    rings = getRings(mol)
    for atom in mol.GetAtoms():
        isRing = any(atom.GetIdx() in ring for ring in rings)
        if atom.GetIdx() not in neighbors and atom.GetIdx() not in visited_atoms:
            fragment_atoms = [(atom.GetSymbol(), '', '')]
            dfsearch(atom, fragment_atoms)
            fragment_sequence = ''.join(f"{symbol}{sub}{double}" for symbol, sub, double in fragment_atoms)
            
            output_sequence = annotate(fragment_sequence,isRing)
            if not fsmiles:
                fsmiles += "'start_0'"
                
            else: 
                fsmiles += "'sep_0'"

            fsmiles += output_sequence
        
    fsmiles += "'sep_0''end_0'"
    return fsmiles

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python main.py <your_input_filepath>.sdf")
        sys.exit(1)

    insdf= sys.argv[1]
    try:
        suppl = Chem.SDMolSupplier(insdf)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
    

    # Write the FSMILES to txt file
    output_folder = 'output'
    os.makedirs(output_folder, exist_ok=True)

    m = suppl[0]
    visualise = True
    if visualise:
        outpath= os.path.join(output_folder, 'example_bondvis.png')
        visualizeFrag(m,outpath)

    fsmiles= fsmlsConverter(m)
    with open(os.path.join(output_folder, 'example_out.txt'), 'w') as f:
        f.write(fsmiles)       


