# pycg_bonds
Allow a cg structure to be visualized in pymol like an atomistic structure.

# Usage: 
in PyMOL do: `run pycg_bonds.py`

# Provides:
cg_bonds [selection], [aa_template], [bfile]

selection   : any selection to act upon (default: all)
aa_template : an aa pdb file to take ss and bfactors from (default: None)
bfile       : alternative file containing a list of backbone b factors (default: None)

Without an aa_template, this function only adds bonds between the backbone beads 
so they can be nicely visualized using line or stick representation.
An aa_template provides a secondary structure assignment that can be used to draw a cartoon representation
The cartoon representation requires "cartoon_trace_atoms" to be set 
because the backbone beads are not recognized as amino acids by pymol.
Sadly this causes the cartoon representations of all structures to also include non backbone atoms.
Therefore this script provides the 'cg_cartoon' function to represent only the backbone atoms as cartoon.

NOTE: Dealing with separate chains is not (yet) implemented.
