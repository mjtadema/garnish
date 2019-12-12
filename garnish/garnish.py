# Copyright 2019-2019 the garnish authors. See copying.md for legal info.

from pymol import cmd  # stored
from glob import glob

# local imports
from .parse import parse
from .system import System
from .utils import get_chain_bb


def garnish(file=None, selection='all', gmx=None, fix_elastics=1, guess_bb=1, show=1):
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
    fix_elastics = fix elastic bonds based on atom id. Disable if your beads are
                   numbered non-sequentially (default: 1)
    guess_bb = if file is not present, simply draw bonds between backbone atoms (default: 1)
    show = adjust representation after drawing bonds (default: 1)
    """
    fix_elastics = bool(int(fix_elastics))
    guess_bb = bool(int(guess_bb))
    show = bool(int(show))

    # Retain order so pymol does not sort the atoms, giving a different result when saving the file
    cmd.set("retain_order", 1)

    if file:
        # parse the file
        sys_dict = parse(file, gmx)
        # create System object and draw all the bonds
        system = System(sys_dict, fix_elastics=fix_elastics)
        system.draw_bonds(selection)
        system.transfer_attributes(selection)

    elif guess_bb:
        bb_beads = get_chain_bb(selection)
        # For each object and chain, draw bonds between BB beads
        for obj, chains in bb_beads.items():
            for _, bbs in chains.items():
                # create bond tuples for "adjacent" backbone beads
                bonds = [(bbs[i], bbs[i+1]) for i in range(len(bbs) - 1)]
                for a, b in bonds:
                    try:
                        cmd.add_bond(obj, a, b)
                    except AttributeError:
                        cmd.bond(f"{obj} and ID {a}", f"{obj} and ID {b}")

    else:
        # show as spheres if no info on bonds is present
        if show:
            cmd.show_as('spheres', selection)
        return

    if show:
        cmd.hide("everything", selection)
        cmd.show_as("sticks", selection)
        # Fix the view for elastics
        cmd.color('orange', '*_elastics')
        cmd.show_as("lines", '*_elastics')


def extend_garnish():
    cmd.extend('garnish', garnish)

    # tab completion for the garnish command
    def useful_file_sc():
        return cmd.Shortcut(glob('**/*.tpr', recursive=True) +
                            glob('**/*.top', recursive=True) +
                            glob('**/*.itp', recursive=True))

    cmd.auto_arg[0]['garnish'] = [useful_file_sc, 'input file', ', ']
    # here object_sc is more informative than selection_sc
    cmd.auto_arg[1]['garnish'] = [cmd.object_sc, 'selection', '']
