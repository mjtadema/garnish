# pycg_bonds
Allow a cg structure to be visualized in pymol like an atomistic structure.

# Screenshot
![Screenshot](/screenshots/screenshots.png?raw=true "lysozyme 2VB1")
Lysozyme 2VB1 from left to right: crystal structure in cartoon representation, licorice representation with backbone bonds drawn, elastic network visualized using line representation.

# Usage: 
in PyMOL do: `run pycg_bonds.py`

# Provides:
cg_bonds *selection*, *tpr_file*

- selection   : any selection to act upon (default: all)
- tpr_file : a .tpr file to extract bond information from (default: None)

Without a tpr file, this function only adds bonds between the backbone beads so they can be
nicely visualized using line or stick representation.
A tpr file provides topology information that can be used to draw side chain and elastic bonds. 

# TIPS

Pymol is extremely slow at creating many bonds at the same when the loaded structures contain huge amounts of atoms.
If pycg_bonds is very slow for you, try to **remove all the waters**.