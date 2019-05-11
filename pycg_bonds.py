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

# Order might be important
cmd.set("retain_order", 1)


def cg_bonds(selection='(all)', aa_template=None):
    """
    Allow a cg structure to be visualized in pymol like an atomistic structure.

    Usage: cg_bonds [selection], [aa_template]

    selection   : any selection to act upon (default: all)
    aa_template : an aa pdb file to take ss and bfactors from (default: None)

    Without an aa_template, this function only adds bonds between the backbone beads
    so they can be nicely visualized using line or stick representation.
    An aa_template provides a secondary structure assignment that can be used to draw a cartoon representation
    The cartoon representation requires "cartoon_trace_atoms" to be set
    because the backbone beads are not recognized as amino acids by pymol.
    Sadly this causes the cartoon representations of all structures to also include non backbone atoms.
    Therefore this script provides the 'cg_cartoon' function to represent only the backbone atoms as cartoon.

    """

    # Fix the view nicely
    cmd.hide("everything", selection)
    cmd.show_as("lines", selection+" and name BB")
    #cmd.color("green", selection)
    cmd.util.cbc(selection)

    # Get all the chain identifiers
    stored.chains = []
    cmd.iterate(str(selection)+" and resi 1 and name BB", "stored.chains.append(chain)")
    chain_bb = {}

    # Store the bb atom IDs for each chain
    for c in stored.chains:
        stored.c_bbs = []
        cmd.iterate(str(selection)+" and name BB and chain {}".format(c), "stored.c_bbs.append(ID)")
        chain_bb[c] = stored.c_bbs

    # For each chain, draw bonds between BB beads
    for c, bbs in chain_bb.items():
        for i in range(len(bbs)-1):
            bb = bbs[i]
            bb_next = bbs[i+1]
            cmd.bond("ID {}".format(bb), "ID {}".format(bb_next))

    # If an atomistic template was also given, extract ss information
    if aa_template:
        cmd.load(aa_template, "aa_template")
        stored.ss = []
        stored.bfactors = []
        cmd.iterate("aa_template and name CA", "stored.ss.append(ss)")
        cmd.iterate("aa_template and name CA", "stored.bfactors.append(b)")
        for bb, ss in zip(stored.bfactors, stored.ss):
            cmd.alter("ID {}".format(bb), 'ss="{}"'.format(ss))
        cmd.delete("aa_template")
        cmd.center(selection)
        cmd.set("cartoon_trace_atoms")
        cg_cartoon(selection)
        cmd.extend('cg_cartoon', cg_cartoon)


def cg_cartoon(selection):
    cmd.cartoon("automatic", selection)
    cmd.show_as("cartoon", selection+" and (name BB or name CA)")


cmd.extend('cg_bonds', cg_bonds)
print(cg_bonds.__doc__)
