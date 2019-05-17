# pycg_bonds
Allow a cg structure to be visualized in pymol like an atomistic structure.

# Usage: 
in PyMOL do: `run pycg_bonds.py`

# Provides:
cg_bonds [selection], [tpr_file]

selection       : any selection to act upon (default: all)
tpr_file        : a tpr file containing the same number of atoms as the selection

Without a selection, the script will work on everything.
Without a tpr file, only backbone bonds will be drawn.
With a tpr file, bonds between all the atoms, including elastic bonds, will be drawn.
