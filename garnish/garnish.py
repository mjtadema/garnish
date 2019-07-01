# Copyright 2019-2019 the garnish authors. See copying.md for legal info.

import os.path
from pymol import cmd  # stored
from glob import glob

# local imports
from .parse_tpr import parse_tpr
from .parse_top import parse_top
from .system import System
from .utils import clean_path, get_chain_bb


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
        system = System(parsed_data)
        system.draw(selection)

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
