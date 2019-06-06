# pycg_bonds
Allow a cg structure to be visualized in pymol like an atomistic structure.

# Usage: 
in PyMOL do: `run pycg_bonds.py`

# Provides:
cg_bonds *selection*, *tpr_file*

- selection   : any selection to act upon (default: all)
- tpr_file : a .tpr file to extract bond information from (default: None)

Without a tpr file, this function only adds bonds between the backbone beads so they can be
nicely visualized using line or stick representation.
A tpr file provides topology information that can be used to draw side chain and elastic bonds. 
