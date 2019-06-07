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
from pymol import cmd, stored
import networkx as nx
from pathlib import Path
import re
import subprocess
import shutil


sys.path.append(os.path.dirname(__file__))


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

    # dump tpr info in a string
    gmxdump = subprocess.run([gmx, "dump", "-s", str(tpr.absolute())],
                             capture_output=True, text=True)

    # define regex patterns for later use
    regexp_info = {
        'molblock': re.compile('^\s+molblock\s+\((\d+)'),
        'moltype': re.compile('^\s+moltype\s+=\s+(\d+)'),
        'molcount': re.compile('^\s+#molecules\s+=\s+(\d+)'),
        'endinfo': re.compile('^\s+ffparams')
    }

    regexp_header = re.compile("^\s+moltype\s+\((\d+)\):")

    regexp_data = {
        'atomnames': re.compile("^\s+atom\[(\d+)\]=\{name=\"(\w+)"),
        'bonds': re.compile("^\s+\d+\s\w+=\d+\s\(BONDS\)\s+(\d+)\s+(\d+)"),
        'constr': re.compile("^\s+\d+\s\w+=\d+\s\(CONSTR\)\s+(\d+)\s+(\d+)"),
        'harmonic': re.compile("^\s+\d+\s\w+=\d+\s\(HARMONIC\)\s+(\d+)\s+(\d+)")
    }

    # initialize some stuff
    mol_blocks = {}
    mols_atom_n = {}
    mols_bonds = {}
    backbone = []
    info_section = True

    # split the dump in lines
    for line in gmxdump.stdout.split('\n'):
        # parse for molecule information in the first section
        if info_section:
            for k, p in regexp_info.items():
                matched = p.match(line)
                if matched:
                    if k == 'molblock':
                        # save current molecule block id
                        curr_block_id = matched.group(1)
                    if k == 'moltype':
                        # save current molecule type id
                        curr_mol_type = matched.group(1)
                    if k == 'molcount':
                        # create new entry in mol_blocks as {block_id: (moltype, molcount)}
                        mol_blocks[curr_block_id] = (curr_mol_type, int(matched.group(1)))
                    if k == 'endinfo':
                        # stop parsing these patterns
                        info_section = False
        else:
            # look for a new molecule definition
            matched = regexp_header.match(line)
            if matched:
                # if molecule header was found, initialize some stuff
                # molecule type id
                curr_mol_type = matched.group(1)
                # corresponding bond dictionary
                mols_bonds[curr_mol_type] = {
                    'bonds': [],
                    'constr': [],
                    'harmonic': []
                }
                # corresponding number of atoms
                mols_atom_n[curr_mol_type] = 0
            for k, p in regexp_data.items():
                matched = p.match(line)
                if matched:
                    if k == 'atomnames':
                        mols_atom_n[curr_mol_type] += 1
                        # save backbone beads for later fix of short elastic bonds
                        if matched.group(2) == "BB":
                            backbone.append(int(matched.group(1)))
                    else:
                        bond = tuple(int(b) for b in matched.group(1, 2))
                        mols_bonds[curr_mol_type][k].append(bond)

    # now use the gathered info to generate a graph with all the bonds
    bonds_dict = {}
    offset = 0
    # iterate over each block of molecules
    for _, data in mol_blocks.items():
        mol_type = data[0]
        mol_count = data[1]
        # initialize the bonds for mol_type, unless they already exists
        bonds_dict[mol_type] = bonds_dict.get(mol_type, {})
        # iterate based on molecule count
        for i in range(mol_count):
            # save bonds, then add to offset based on number of atoms in mol_type
            for bond_type, bond_list in mols_bonds[mol_type].items():
                # initialize list of bonds for bond_type, unless it already exists
                bonds_dict[mol_type][bond_type] = bonds_dict[mol_type].get(bond_type, [])
                for bond in bond_list:
                    shifted_bond = tuple(atom_id + offset for atom_id in bond)
                    bonds_dict[mol_type][bond_type].append(shifted_bond)
            offset += mols_atom_n[mol_type]

    # move short range elastic bonds from `bonds` to `elastic` (protein fix)
    short_elastic = {
        mol: [] for mol in mol_blocks
    }

    for mol, bonds in bonds_dict.items():
        for b in bonds['bonds']:
            # check if both beads are in backbone list
            if b[0] in backbone and b[1] in backbone:
                # if they're not adjacent, move bond to elastic and remove from here
                if abs(backbone.index(b[0]) - backbone.index(b[1])) > 1:
                    short_elastic[mol].append(b)

    for mol, bonds in short_elastic.items():
        for b in bonds:
            bonds_dict[mol]['bonds'].remove(b)
            bonds_dict[mol]['harmonic'].append(b)

    # Convert the lists of bonds to networkx graphs
    bond_graphs = {}
    for mol, bonds in bonds_dict.items():
        bond_graphs[mol] = {}
        for bondtype, bond_list in bonds.items():
            g = nx.Graph()
            g.add_edges_from(list(bond_list))
            bond_graphs[mol][bondtype] = g

    return bond_graphs


def rel_atom(selection):
    # Make a dict of all the atoms (to get effective relative atom numbering)
    rel_atom_dict = {}
    atoms = cmd.get_model(selection)
    for i, at in enumerate(atoms.atom):
        rel_atom_dict[i] = at.index
    return rel_atom_dict


def cg_bonds(*args, **kwargs):  # selection='(all)', tpr_file=None): #aa_template=None):
    """
    Allow a cg structure to be visualized in pymol like an atomistic structure.

    Usage: cg_bonds [selection], [tpr_file]

    selection   : any selection to act upon (default: all)
    tpr_file : a .tpr file to extract bond information from (default: None)

    Without a tpr file, this function only adds bonds between the backbone beads so they can be
    nicely visualized using line or stick representation.
    A tpr file provides topology information that can be used to draw side chain and elastic bonds.
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
    cmd.show_as("lines", selection)
    cmd.util.cbc(selection)

    # Draw all the bonds based on tpr file
    if tpr_file:
        # warn at the end if some atoms in the tpr are missing from the structure
        warn = False
        # Get bond graphs
        bond_graphs = parse_tpr(tpr_file)
        # Create dummy object to draw elastic bonds in
        elastics_selector = selection+"_elastics"
        cmd.create(elastics_selector, selection)
        # Make a dict of all the atoms (to get effective relative atom numbering)
        rel_atom_selection = rel_atom(selection)
        # Draw all the bonds
        for btype in ['bonds', 'constr']:
            for _, bonds in bond_graphs.items():
                for a, b in bonds[btype].edges:
                    try:
                        a = rel_atom_selection[a]
                        b = rel_atom_selection[b]
                        cmd.bond(f"{selection} and ID {a}", f"{selection} and ID {b}")
                    except KeyError:
                        warn = True
            # Get relative atoms for elastics object
        rel_atom_elastics = rel_atom(elastics_selector)
        atoms = cmd.get_model(elastics_selector)
        for i, at in enumerate(atoms.atom):
            rel_atom_elastics[i] = at.index
        # Draw elastic network
        for _, bonds in bond_graphs.items():
            for a, b in bonds['harmonic'].edges:
                try:
                    a = rel_atom_elastics[a]
                    b = rel_atom_elastics[b]
                    cmd.bond(f"{elastics_selector} and ID {a}", f"{elastics_selector} and ID {b}")
                except KeyError:
                    warn = True
        cmd.color("orange", elastics_selector)

        # warn about missing atoms if needed.
        if warn:
            print('WARNING: some atoms present in the tpr file were not found in the loaded '
                  'structure.\n Bonds containing those atoms were not drawn.')

    # Draw simple bonds between backbone beads
    else:
        chain_bb = get_chain_bb(selection)
        # For each chain, draw bonds between BB beads
        for _, bbs in chain_bb.items():
            bonds = [(bbs[i], bbs[i+1]) for i in range(len(bbs) - 1)]
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


#def cg_cartoon(selection):
#    cmd.cartoon("automatic", selection)
#    cmd.show_as("cartoon", selection + " and (name BB or name CA)")


cmd.extend('cg_bonds', cg_bonds)
