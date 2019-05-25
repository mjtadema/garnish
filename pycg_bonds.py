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
import io
import subprocess
import shlex
import shutil
from pdb import set_trace

# Order might be important
cmd.set("retain_order", 1)      # TODO: move to a better place


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


def parse_tpr(tpr_file, gmx=False):
    """
    Parses the gmx dump output of a tpr file into a networkx graph representation of the connectivity within the system

    input: a filename pointing to a tpr file

    returns: a dictionary of molecules, each with a dictonary of graphs representing different connection types

    molecules {
        molid: bondtypes {
            bonds:      nx.Graph
            constr:     nx.Graph
            harmonic:   nx.Graph
            }
        }

    """
    tpr = Path(tpr_file)
    assert tpr.is_file()

    # get gmx executable
    if not gmx:
        gmx = shutil.which('gmx')
    if not gmx:
        raise FileNotFoundError('no gromacs executable found. Add it manually with gmx="PATH_TO_GMX"')

    gmxdump = gmx + " dump -s " + str(tpr.absolute())
    gmxdump = subprocess.Popen(shlex.split(gmxdump), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Regex like cg_bonds to get relevant info
    p_grep = re.compile(".*\#atoms|.*\#beads.*|.*moltype.*|.*\#molecules.*|.*\(BONDS\).*|.*\(CONSTR\).*|.*\(HARMONIC\).*")

    regexp_all = {
        'molid': re.compile("^\s+moltype\s+=\s+(\d+)"),
        'occs': re.compile("^\s+\#molecules\s+=\s+(\d+)"),
        'bead_nr': re.compile("^\s+\#atoms_mol\s+=\s+(\d+)"),
        'total_beads': re.compile("^\s+\#atoms+\s+=+\s+(\d+)"),
    }
    regexp_bonds = {
        'bonds': re.compile("^\s+\d+\s\w+=\d+\s\(BONDS\)\s+(\d+)\s+(\d+)"),
        'constr': re.compile("^\s+\d+\s\w+=\d+\s\\(CONSTR\)\s+(\d+)\s+(\d+)"),
        'harmonic': re.compile("^\s+\d+\s\w+=\d+\s\\(HARMONIC\)\s+(\d+)\s+(\d+)")
    }
    regexp_is_mol = re.compile("^\s+moltype\s+\((\d+)\):")
    
    # Should probably return this as well
    regex_data = {
        k: []
        for k in regexp_all.keys()
    }

    molecules = {}
    
    reading_header = True
    for line in io.TextIOWrapper(gmxdump.stdout, encoding="utf-8"):
        # Filter for the relevant info
        if p_grep.match(line):
            # when looking for a header, only care about these regexes
            if reading_header:
                # Parse the meta info
                for k, p in regexp_all.items():
                    matched = p.match(line)
                    if matched:
                        regex_data[k] = matched.group(1)
                matched_mol = regexp_is_mol.match(line)
                if matched_mol:
                    # If it started to describe a molecule, flag reading_header False
                    # and initialise the first molecule
                    reading_header = False
                    molid = matched_mol.group(1)
                    bonds = {
                        k: []
                        for k in regexp_bonds
                    }
            else:
                # parse bond data
                for k, p in regexp_bonds.items():
                    matched = p.match(line)
                    if matched:
                        bond = matched.group(1, 2)
                        # Cast to int
                        bond = tuple( int(b) for b in bond )
                        bonds[k].append(bond)
                        # no need to parse for everything
                        break
                # if none of the above was found, look for a header
                else:
                    molecules[molid] = bonds
                    molid = regexp_is_mol.search(line).group(1)
                    bonds = {
                        k: []
                        for k in regexp_bonds
                    }

    # Convert the lists of bonds to graphs
    for molid, molecule in molecules.items():
        for bondtype, bonds in molecule.items():
            g = nx.Graph()
            g.add_edges_from(list(bonds))
            molecules[molid][bondtype] = g
    
    return molecules


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
    An aa_template provides a secondary structure assignment that can be used to draw a cartoon representation
    The cartoon representation requires "cartoon_trace_atoms" to be set
    because the backbone beads are not recognized as amino acids by pymol.
    Sadly this causes the cartoon representations of all structures to also include
    non backbone atoms.
    Therefore this script provides the 'cg_cartoon' function to represent only
    the backbone atoms as cartoon.
    """

    selection = '(all)'
    tpr_file = None

    for arg in args:
        maybe_tpr = Path(str(arg))
        if maybe_tpr.is_file() and maybe_tpr.suffix == ".tpr":
            tpr_file = maybe_tpr
        else:
            # arg is probably a selection
            selection = arg

    #if not tpr_file and Path("./topol.tpr").exists():
    #    tpr_file = Path("./topol.tpr")


    # Order might be important
    cmd.set("retain_order", 1)

    # Fix the view nicely
    cmd.hide("everything", selection)
    cmd.show_as("lines", selection + " and name BB")
    cmd.util.cbc(selection)

    ## Get all the chain identifiers and all the atoms
    #chains = cmd.get_chains(selection)
    #atoms_per_chain = {}
    #for chain in chains:
    #    model = cmd.get_model(selection+" and chain "+chain)
    #    atoms_per_chain[chain] = [
    #                at.index
    #                for at in model.atom
    #            ]

    if tpr_file:
        # Get bond graphs
        bond_graphs = parse_tpr(tpr_file)
        # Create dummy object to draw elastic bonds in
        elastics_selector = selection+"_elastics"
        cmd.create(elastics_selector, selection)
        # Make a dict of all the atoms (to get effective relative atom numbering)
        rel_atom_selection = rel_atom(selection)
        # Draw all the bonds
        for mol in bond_graphs.values():
            for btype in ['bonds','constr']:
                for a, b in mol[btype].edges:
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

    else:
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


def cg_cartoon(selection):
    cmd.cartoon("automatic", selection)
    cmd.show_as("cartoon", selection + " and (name BB or name CA)")


cmd.extend('cg_bonds', cg_bonds)
