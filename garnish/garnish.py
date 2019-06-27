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

import os.path
from pymol import cmd  # stored
import networkx as nx
import numpy as np
from glob import glob

# local imports
from .parse_tpr import parse_tpr
from .parse_top import parse_top
from .utils import clean_path, get_chain_bb


def make_graphs(system):
    """
    uses data gathered from file parsing to correctly identify bonds for each molecule

    returns graph representations of all the bonds in the system
    """

    # Iterate over molecules, fix elastic bonds where necessary
    for key in system['topology'].keys():
        molecule = system['topology'][key]
        bond_list = molecule['connectivity']['bonds']
        tmp_harmonics = molecule['connectivity']['harmonic']
        tmp_bonds = []
        backbone = molecule['backbone']

        while bond_list:
            try:
                bond = bond_list.pop()
            except IndexError as e:
                raise e
            bond_type = tmp_bonds
            # if both atoms in a bond are labeled as backbone, go deeper
            if all(atom in backbone for atom in bond):
                at_a, at_b = bond
                at_a, at_b = int(at_a), int(at_b)
                atom1_idx = backbone.index(at_a)
                atom2_idx = backbone.index(at_b)
                # if the bond is between two non-adjacent backbone beads, move it to `elastic`
                if abs(atom1_idx - atom2_idx) > 1:
                    bond_type = tmp_harmonics
            bond_type.append(bond)
        # update old dictionaries with new info
        molecule['connectivity']['bonds'] = tmp_bonds
        molecule['connectivity']['harmonic'] = tmp_harmonics

    # Convert the lists of bonds to networkx graphs
    bond_graphs = {}
    offset = 0
    for block_id, block in system['blocks'].items():
        moltype = block['moltype']
        n_mol = block['n_molecules']
        n_at = system['topology'][moltype]['n_atoms']
        connectivity = system['topology'][moltype]['connectivity']
        # transform bonds in numpy arrays to easily apply offset
        connectivity = {btype: np.array(bonds) for btype, bonds in connectivity.items()}

        # repeat for each occurence of molecule in this block
        for i in range(n_mol):
            # unique key
            key = f"{block_id}_{i}"
            bond_graphs[key] = {}
            for btype, bonds in connectivity.items():
                # create a graph based on connectivity and offset it to match atom numbers
                g = nx.Graph()
                g.add_edges_from(bonds + offset)
                bond_graphs[key][btype] = g
            # shift offset by how many atoms this molecule has
            offset += n_at

    return bond_graphs


def garnish(file=None, selection='all', gmx=None):
    """
DESCRIPTION

    Allow a coarse grained structure to be visualized in pymol like an atomistic structure
    by drawing bonds and elastic network.

    Without a top/tpr file, this function only adds bonds between the backbone beads
    so they can be nicely visualized using line or stick representation.
    Adding a top/tpr file provides topology information that can be used
    to draw side chain and elastic bonds.

USAGE

    garnish [file [, selection [, gmx]]]

ARGUMENTS

    file = a tpr or topology file to extract bond information from (default: None)
    selection = any selection to act upon (default: all)
    gmx = gmx executable path (default: inferred by `which gmx`)
    """
    tpr_file = None
    top_file = None
    aa_template = None

    if file == None:
        pass
    else:
        maybe_file = clean_path(file)
        if os.path.isfile(maybe_file):
            ext = os.path.splitext(maybe_file)[1]
            if ext == ".tpr":
                tpr_file = maybe_file
            elif ext == ".top" or ext == ".itp":
                top_file = maybe_file
            elif ext == ".pdb":
                aa_template = maybe_file
            else:
                raise TypeError(f'"{ext}" is not a supported format.')
        else:
            raise FileNotFoundError(f'{maybe_file} does not exist.')

    # Retain order so pymol does not sort the atoms, giving a different result when saving the file
    cmd.set("retain_order", 1)

    # Draw all the bonds based on tpr or top file
    if tpr_file or top_file:
        # warn at the end if some atoms in the tpr are missing from the structure
        warn = False
        # Get bond graphs
        if tpr_file:
            parsed_data = parse_tpr(tpr_file, gmx)
        elif top_file:
            parsed_data = parse_top(top_file)
        bond_graphs = make_graphs(parsed_data)
        selection_objects = cmd.get_object_list(selection)
        for obj in selection_objects:
            elastics_obj = obj+"_elastics"
            # Create dummy object to draw elastic bonds in
            cmd.copy(elastics_obj, obj)
            # Make a dict of all the atoms (to get effective relative atom numbering)
            # Draw all the bonds
            for btype in ['bonds', 'constr']:
                for _, bonds in bond_graphs.items():
                    for a, b in bonds[btype].edges:
                        try:
                            a += 1
                            b += 1
                            try:
                                cmd.add_bond(obj, a, b)
                            except AttributeError:
                                cmd.bond(f"({obj} and ID {a})", f"({obj} and ID {b})")
                        except KeyError:
                            warn = True
            # Get relative atoms for elastics object
            atoms = cmd.get_model(elastics_obj)
            # Draw elastic network
            for _, bonds in bond_graphs.items():
                for i, (a, b) in enumerate(bonds['harmonic'].edges):
                    try:
                        a += 1
                        b += 1
                        try:
                            cmd.add_bond(elastics_obj, a, b)
                        except AttributeError:
                            cmd.bond(f"({elastics_obj} and ID {a})", f"({elastics_obj} and ID {b})")
                    except KeyError:
                        warn = True
            cmd.color("orange", elastics_obj)

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

    # Fix the view nicely
    cmd.hide("everything", selection)
    cmd.show_as("sticks", selection)
    cmd.show_as("lines", '*_elastics')


#    # If an atomistic template was given, extract ss information
#    if aa_template:
#        cmd.load(aa_template, "aa_template")
#        stored.ss = []
#        stored.bfactors = []
#        cmd.iterate("aa_template and name CA", "stored.ss.append(ss)")
#        cmd.iterate("aa_template and name CA", "stored.bfactors.append(b)")
#        for bb, ss in zip(stored.bfactors, stored.ss):
#            cmd.alter(f"ID {bb}", f'ss="{ss}"')
#        cmd.delete("aa_template")
#        cmd.center(selection)
#        cmd.set("cartoon_trace_atoms")
#        cmd.cartoon("automatic", selection)
#        cmd.show_as("cartoon", selection + " and (name BB or name CA)")


def extend_garnish():
    cmd.extend('garnish', garnish)

    # tab completion for the garnish command
    useful_file_sc = lambda: cmd.Shortcut(
         glob('*/') + glob('*.tpr') + glob('*.top') + glob('*.itp')  # + glob('*.pdb')
    )
    cmd.auto_arg[0]['garnish'] = [useful_file_sc, 'input file', ', ']
    # here object_sc is more informative than selection_sc
    cmd.auto_arg[1]['garnish'] = [cmd.object_sc, 'selection', '']


## make sure it can be run as a script for simplicity
#if __name__ == 'pymol':
#    load()
