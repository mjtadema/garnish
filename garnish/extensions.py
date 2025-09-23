"""
Contains PyMOL extensions
Is automatically imported when run from PyMOL
Also automatically extends garnish to be called in PyMOL
"""
# Copyright 2019-2019 the garnish authors. See copying.md for legal info.

from pymol import cmd
from glob import glob
import logging

# local imports
from . import Parser


def bonds_from_system(selection, system, _self=cmd, quiet=1):
    """
    Take a System object and draw all the bonds in selection
    :param system:
    :param _self:
    :param quiet:
    :return:
    """

    for obj in _self.get_object_list(selection):
        elastics_obj = obj + "_elastics"
        _self.copy(elastics_obj, obj)
        idx_shift = 0
        for mol, n in system.topology:
            for i in range(n):
                draw_bonds(obj, mol.bonds, idx_shift)
                draw_bonds(elastics_obj, mol.elastics, idx_shift)
                idx_shift += len(mol)


def draw_bonds(target, connections, idx_shift,
               _self=cmd, quiet=1):
    for a, b in connections:
        a += idx_shift
        b += idx_shift
        _self.add_bond(target, a, b)


@cmd.extend
def garnish(file="topol.top", selection='all', show=1, _self=cmd, quiet=1):
    """
DESCRIPTION

    Allow a coarse grained structure to be visualized in pymol like an atomistic structure
    by drawing bonds and elastic network.

    Using a top/tpr file provides topology information that can be used
    to draw side chain and elastic bonds.

USAGE

    garnish [file [, selection [, gmx]]]

ARGUMENTS

    file = a tpr or topology file to extract bond information from (default: None)
    selection = any selection to act upon (default: all)
    gmx = gmx executable path (default: inferred by `which gmx`)
    fix_elastics = fix elastic bonds based on atom id. Disable if your beads are
                   numbered non-sequentially (default: 1)
    guess_prot = if file is not present, simply draw bonds between backbone atoms of a protein (default: 1)
    show = adjust representation after drawing bonds (default: 1)
    """
    show = bool(int(show))

    # Retain order so pymol does not sort the atoms, giving a different result when saving the file
    _self.set("retain_order", 1)

    if int(quiet) < 0:
        logging.basicConfig(level="DEBUG")
        logging.debug("Starting debug logging")

    if file:
        # parse the file and create system object
        system = Parser(top_file=file).run()
        # Draw all the bonds in target
        bonds_from_system(selection, system)
    else:
        # show as spheres if no info on bonds is present
        if show:
            _self.show_as('spheres', selection)
        return

    if show:
        _self.hide("everything", selection)
        _self.show_as("sticks", selection)
        # Fix the view for elastics
        _self.color('orange', '*_elastics')
        _self.show_as("lines", '*_elastics')


def useful_file_sc():
    """tab completion for the garnish command"""
    return cmd.Shortcut(glob('**/*.top', recursive=True) +
                        glob('**/*.itp', recursive=True))


cmd.auto_arg[0]['garnish'] = [useful_file_sc, 'input file', ', ']
# here object_sc is more informative than selection_sc
cmd.auto_arg[1]['garnish'] = [cmd.object_sc, 'selection', '']
