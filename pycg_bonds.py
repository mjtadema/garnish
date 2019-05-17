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
import string
import networkx as nx
from pathlib import Path
import re,io
import subprocess, shlex

# Order might be important
cmd.set("retain_order", 1)

def get_chain_bb(selection, chains):
    """
    returns dictionary with format {chain: list-of-bb-atoms}
    """
    chain_bb = {}
    for c in chains:
        # if chain is empty string, put it in the "all" bin
        if not c:
            c = "all"
        chain_bb[c] = cmd.identify(selection + f" and chain {c} and name BB")
    return chain_bb

#def get_chain_bb(selection, chains):
#    chain_bb = {}
#    for c in chains:
#        if c in string.ascii_letters:
#            stored.c_bbs = []
#            cmd.iterate(str(selection)+" and name BB and chain {}".format(c), "stored.c_bbs.append(ID)")
#            chain_bb[c] = stored.c_bbs
#        # If there are no ids, put them together
#        else:
#            stored.c_bbs = []
#            cmd.iterate(str(selection)+" and name BB", "stored.c_bbs.append(ID)")
#            chain_bb["all"] = stored.c_bbs
#            break
#    return chain_bb

def parse_tpr(tpr_file):
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

    gmxdump = "/usr/bin/gmx dump -s "+str(tpr.absolute())
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
            if reading_header:
                # Parse the meta info
                for k, p in regexp_all.items():
                    if p.match(line):
                        regex_data[k] = p.findall(line)[0]
                # If it started to describe a molecule, flag reading_header False
                if regexp_is_mol.match(line):
                    reading_header = False
                    molid = regexp_is_mol.findall(line)[0]
                    bonds = {
                        k:[]
                        for k in regexp_bonds
                    }
            else:
                # Check if a line is a new molecule
                if not regexp_is_mol.match(line):
                    for k, p in regexp_bonds.items():
                        if p.match(line):
                            # Cast to int and increment with one to match the numbering in pymol
                            bond = p.findall(line)[0]
                            bond = tuple( int(b)+1 for b in bond )
                            bonds[k].append(bond)
                # If not, parse the bonds
                else:
                    molecules[molid] = bonds
                    molid = regexp_is_mol.findall(line)[0]
                    bonds = {
                        k:[]
                        for k in regexp_bonds
                    }
                    
    # Convert the lists of bonds to graphs
    for molid, molecule in molecules.items():
        for bondtype, bonds in molecule.items():
            g = nx.Graph()
            g.add_edges_from(list(bonds))
            molecules[molid][bondtype] = g
    
    return molecules


def cg_bonds(selection='(all)', tpr_file=None): #aa_template=None):
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

    # Get molecules
    molecules = parse_tpr(tpr_file)
    
    if tpr_file:
        # Draw all the bonds
        for mol in molecules.values():
            for btype in mol.values():
                for a, b in btype.edges:
                    cmd.bond(f"ID {a}", f"ID {b}")

    else:
        chain_bb = get_chain_bb(selection, chains)
        # For each chain, draw bonds between BB beads
        for _, bbs in chain_bb.items():
            for i in range(len(bbs)-1):
                bb = bbs[i]
                bb_next = bbs[i+1]
                cmd.bond(f"ID {bb}", f"ID {bb_next}")

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
