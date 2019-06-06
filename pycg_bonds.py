# MIT License
#
# Copyright (c) 2019 Matthijs Tadema, Lorenzo Gaifas
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

import os.path, sys
sys.path.append(os.path.dirname(__file__))

from pymol import cmd, stored
import networkx as nx
from pathlib import Path
import re
import subprocess
import shutil

def get_chain_bb(selection):
    """
    returns dictionary with format {chain: list-of-bb-atoms}
    """
    chains = cmd.get_chains(selection)
    chain_bb = {}
    for c in chains:
        # if chain is empty string, put it in the "all" bin
        if not c:
            c = "*"
        # "chain all" didn't actually select anything
        bb_id = cmd.identify(selection + f" and chain {c} and name BB")
        chain_bb[c] = bb_id
    return chain_bb


def get_gmx(gmx_bin):
    """
    if gmx binary is not given, find it. If it can't be found, raise an exception
    """
    if not gmx_bin:
        gmx_bin = shutil.which('gmx')
    if not gmx_bin:
        raise FileNotFoundError('no gromacs executable found.'
                                'Add it manually with gmx="PATH_TO_GMX"')
    return gmx_bin


def parse_tpr(tpr_file, gmx=False):
    """
    Parses the gmx dump output of a tpr file into a networkx graph representation
    of the connectivity within the system

    input: a filename pointing to a tpr file

    returns: a dictionary of graphs representing different connection types

    bond_graphs {
        bonds:      nx.Graph
        constr:     nx.Graph
        harmonic:   nx.Graph
    }

    """
    tpr = Path(tpr_file)
    assert tpr.is_file()

    gmx = get_gmx(gmx)

    gmxdump = subprocess.run([gmx, "dump", "-s", str(tpr.absolute())],
                             capture_output=True, text=True)

    regexp_header = re.compile("^\s+moltype\s+\((\d+)\):")

    regexp_data = {
        'atomnames': re.compile("^\s+atom\[(\d+)\]=\{name=\"BB\"\}"),
        'bonds': re.compile("^\s+\d+\s\w+=\d+\s\(BONDS\)\s+(\d+)\s+(\d+)"),
        'constr': re.compile("^\s+\d+\s\w+=\d+\s\(CONSTR\)\s+(\d+)\s+(\d+)"),
        'harmonic': re.compile("^\s+\d+\s\w+=\d+\s\(HARMONIC\)\s+(\d+)\s+(\d+)")
    }

    bond_graphs = {
        'bonds': [],
        'constr': [],
        'harmonic': []
    }

    backbone = []

    looking_for_header = True
    for line in gmxdump.stdout.split('\n'):
        # Skip as much lines as possible to be faster
        if looking_for_header:
            matched = regexp_header.match(line)
            if matched:
                # if molecule header was found, start looking for actual data
                looking_for_header = False
        elif not looking_for_header:
            for k, p in regexp_data.items():
                matched = p.match(line)
                if matched:
                    if k == 'atomnames':
                        # save backbone beads for later fix of short elastic bonds
                        backbone.append(int(matched.group(1)))
                    else:
                        bond = tuple( int(b) for b in matched.group(1, 2))
                        bond_graphs[k].append(bond)

    # move short range elastic bonds from `bonds` to `elastic` (protein fix)
    short_elastic = []
    for b in bond_graphs['bonds']:
        # check if both beads are in backbone list
        if b[0] in backbone and b[1] in backbone:
            # if they're not adjacent, move bond to elastic and remove from here
            if abs(backbone.index(b[0]) - backbone.index(b[1])) > 1:
                short_elastic.append(b)

    for b in short_elastic:
        bond_graphs['bonds'].remove(b)
        bond_graphs['harmonic'].append(b)

    # Convert the lists of bonds to networkx graphs
    for bondtype, bonds in bond_graphs.items():
        g = nx.Graph()
        g.add_edges_from(list(bonds))
        bond_graphs[bondtype] = g

    return bond_graphs


def rel_atom(selection):
    # Make a dict of all the atoms (to get effective relative atom numbering)
    rel_atom_dict = {}
    atoms = cmd.get_model(selection)
    for i, at in enumerate(atoms.atom):
        rel_atom_dict[i] = at.index
    return rel_atom_dict


def cg_bonds(*args, **kwargs): #selection='(all)', tpr_file=None): #aa_template=None):
    """
    Allow a cg structure to be visualized in pymol like an atomistic structure.

    Usage: cg_bonds [selection], [aa_template]

    selection   : any selection to act upon (default: all)
    aa_template : an aa pdb file to take ss and bfactors from (default: None)

    Without an aa_template, this function only adds bonds between the backbone beads
    so they can be nicely visualized using line or stick representation.
    An aa_template provides a secondary structure assignment that can be used to draw a
    cartoon representation.
    The cartoon representation requires "cartoon_trace_atoms" to be set
    because the backbone beads are not recognized as amino acids by pymol.
    Sadly this causes the cartoon representations of all structures to also include
    non backbone atoms.
    Therefore this script provides the 'cg_cartoon' function to represent only
    the backbone atoms as cartoon.
    """

    selection = 'all'
    tpr_file = None

    for arg in args:
        maybe_tpr = Path(str(arg))
        if maybe_tpr.is_file() and maybe_tpr.suffix == ".tpr":
            tpr_file = maybe_tpr
        else:
            # arg is probably a selection
            selection = arg

    # Order might be important
    cmd.set("retain_order", 1)

    # Fix the view nicely
    cmd.hide("everything", selection)
    cmd.show_as("lines", selection + " and name BB")
    cmd.util.cbc(selection)

    if tpr_file: # Draw all the bonds based on tpr file
        # Get bond graphs
        bond_graphs = parse_tpr(tpr_file)
        # Create dummy object to draw elastic bonds in
        elastics_selector = selection+"_elastics"
        cmd.create(elastics_selector, selection)
        # Make a dict of all the atoms (to get effective relative atom numbering)
        rel_atom_selection = rel_atom(selection)
        # Draw all the bonds
        for btype in ['bonds', 'constr']:
            for a, b in bond_graphs[btype].edges:
                a = rel_atom_selection[a]
                b = rel_atom_selection[b]
                cmd.bond(f"{selection} and ID {a}", f"{selection} and ID {b}")
            # Get relative atoms for elastics object
        rel_atom_elastics = rel_atom(elastics_selector)
        atoms = cmd.get_model(elastics_selector)
        for i, at in enumerate(atoms.atom):
            rel_atom_elastics[i] = at.index
        # Draw elastic network
        for a, b in bond_graphs['harmonic'].edges:
            a = rel_atom_elastics[a]
            b = rel_atom_elastics[b]
            cmd.bond(f"{elastics_selector} and ID {a}", f"{elastics_selector} and ID {b}")
        cmd.color("orange", elastics_selector)

    else: # Draw simple bonds between backbone beads
        chain_bb = get_chain_bb(selection)
        # For each chain, draw bonds between BB beads
        for _, bbs in chain_bb.items():
            bonds = [ (bbs[i], bbs[i+1]) for i in range(len(bbs) - 1) ]
            for a, b in bonds:
                cmd.bond(f"ID {a}", f"ID {b}")

    ## If an atomistic template was also given, extract ss information
    #if aa_template:
    #    cmd.load(aa_template, "aa_template")
    #    stored.ss = []
    #    stored.bfactors = []
    #    cmd.iterate("aa_template and name CA", "stored.ss.append(ss)")
    #    cmd.iterate("aa_template and name CA", "stored.bfactors.append(b)")
    #    for bb, ss in zip(stored.bfactors, stored.ss):
    #        cmd.alter(f"ID {bb}", f'ss="{ss}"')
    #    cmd.delete("aa_template")
    #    cmd.center(selection)
    #    cmd.set("cartoon_trace_atoms")
    #    cg_cartoon(selection)
    #    cmd.extend('cg_cartoon', cg_cartoon)

#
#def cg_cartoon(selection):
#    cmd.cartoon("automatic", selection)
#    cmd.show_as("cartoon", selection + " and (name BB or name CA)")


cmd.extend('cg_bonds', cg_bonds)
