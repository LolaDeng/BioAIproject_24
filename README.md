# FSMILES converter

This repository provides a tool for converting from an SDF file to FSMILES representation. FSMILES is a modified version of SMILES that reorganises molecules into fragments using a specific syntax.

Original Paper: Feng, W., Wang, L., Lin, Z. et al. Generation of 3D molecules in pockets via a language model. Nat Mach Intell 6, 62–73 (2024). https://doi.org/10.1038/s42256-023-00775-6 


## About FSMILES

FSMILES reorganises molecules into fragments using the standard SMILES syntax for each fragment. The entire molecule is then constructed by combining these fragments using a specific syntax in a fragment first then depth-first manner. 

The FSMILES construction process occurs in two steps:

1. Fragment Division: The ligand is divided into fragments according to the cutting rule.

2. Ring Information Embedding: Ring information is embedded in each FSMILES token, with the number following the element type's underscore indicating the ring size. The symbol "*" denotes the connection points of a fragment, and the preceding atom indicates the connection position. 

The key advantages of FSMILES are outlined as follows:

1. Enhanced 2D Pattern Learning: FSMILES utilises symbols to represent fragments and local structures, leading to improved 2D pattern learning.

2. Prioritization of Ring Closure: FSMILES prioritizes ring closure during molecule generation, ensuring the accurate generation of molecules with precise ring structures and bond angles.



## Usage

1. **Prep**

conda create -n fsmls python=3.8
conda activate fsmls
conda install conda-forge::rdkit

2. **Run the Converter**
Execute the converter script by providing the SDF file as an argument:
```bash
python main.py <your_input_file>.sdf
```

3.**Output**
The converted FSMILES will be saved into a txt file found in the output folder

## Notes

- I coded this little thing in one night it might contain hardcoded parts and will need further improvement.
- Feel free to report issues, or suggest enhancements.
