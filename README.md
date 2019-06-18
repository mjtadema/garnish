# pycg_bonds
Allow a cg structure to be visualized in pymol like an atomistic structure by drawing bonds and elastic network.

# Screenshot
![Screenshot](/screenshots/screenshots.png?raw=true "lysozyme 2VB1")
Lysozyme 2VB1 from left to right: crystal structure in cartoon representation, licorice representation with backbone bonds drawn, elastic network visualized using line representation.

# Set up
To access the function, in PyMOL do `run pycg_bonds.py` or add `run pycg_bonds.py`
to your pymolrc.

# Usage 
```
cg_bonds [file [, selection [, gmx]]]
```
- file = a tpr or topology file to extract bond information from (default: None)
- selection = any selection to act upon (default: all)
- gmx = gmx executable path (default: inferred by `which gmx`)

Without a top/tpr file, this function only adds bonds between the backbone beads
so they can be nicely visualized using line or stick representation.
Adding a top/tpr file provides topology information that can be used
to draw side chain and elastic bonds.

## Examples

To nicely show your martini `system` in pymol:
```
cg_bonds topol.top
```
... and wait. No, it did not crash.

# TIPS
Pymol is extremely slow at creating many bonds at the same when the loaded structures contain huge amounts of atoms (we're waiting for upstream to fix this).
If pycg_bonds is very slow for you, try to **remove all the waters** and any other atom you don't care about (ions...).
