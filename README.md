# Garnish Accurately Represents a Network-Interconnected System in a Hearbeat
GARNISH allows a coarse grained structure to be visualized in pymol like an atomistic structure by drawing bonds and elastic network.

# Screenshot
![Screenshot](/screenshots/screenshots.png?raw=true "lysozyme 2VB1")
Lysozyme 2VB1 from left to right: crystal structure in cartoon representation, licorice representation with backbone bonds drawn, elastic network visualized using line representation.

# Installation
```
pip install git+git://github.com/mjtadema/garnish.git#egg=garnish
```

To access the function, in PyMOL do `import garnish`. 

# Usage 
```
garnish [file [, selection [, gmx]]]
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
garnish topol.top
```
