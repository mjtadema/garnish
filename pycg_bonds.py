# MIT License
# 
# Copyright (c) 2019 Matthijs Tadema
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from pymol import cmd, stored
def cg_cartoon(selection):
    cmd.cartoon("automatic", selection)
    cmd.show_as("cartoon", selection+" and (name BB or name CA)")
def norm_blist(blist):
    minb = min(blist)
    maxb = max(blist)
    try:
        return [
            (x - minb) / (maxb - minb)
            for x in blist
        ]
    except ZeroDivisionError:
        # Apparently the bfactors are zero
        # Just return the original list
        return blist
def cg_bonds(selection='(all)', aa_template=None, bfile=None):
    """
    Allow a cg structure to be visualized in pymol like an atomistic structure.

    Usage: cg_bonds [selection], [aa_template], [bfile]

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

    """
    cmd.hide("everything", selection)
    cmd.show_as("lines", selection+" and name BB")
    cmd.color("green", selection)
    stored.bb_atoms = []
    cmd.iterate(str(selection)+" and name BB", "stored.bb_atoms.append(ID)")
    stored.bfactors = []
    for i in range(len(stored.bb_atoms)-1):
        bb = stored.bb_atoms[i]
        bb_next = stored.bb_atoms[i+1]
        cmd.bond("ID {}".format(bb), "ID {}".format(bb_next))
    if aa_template:
        cmd.load(aa_template, "aa_template")
        stored.ss = []
        cmd.iterate("aa_template and name CA", "stored.ss.append(ss)")
        cmd.iterate("aa_template and name CA", "stored.bfactors.append(b)")
        for bb, ss in zip(stored.bb_atoms, stored.ss):
            cmd.alter("ID {}".format(bb), 'ss="{}"'.format(ss))
        cmd.delete("aa_template")
        cmd.center(selection)
        cmd.set("cartoon_trace_atoms")
        cg_cartoon(selection)
        cmd.extend('cg_cartoon', cg_cartoon)
    if bfile:
        with open(bfile, 'r') as f:
            stored.bfactors = [
                float(b)
                for b in f
            ]
    if aa_template or bfile:
        stored.bfactors = norm_blist(stored.bfactors)
        for bb, b in zip(stored.bb_atoms, stored.bfactors):
            cmd.alter("ID {}".format(bb), "b={}".format(b))
cmd.extend('cg_bonds', cg_bonds)
print(cg_bonds.__doc__)
