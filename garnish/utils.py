# Copyright 2019-2019 the garnish authors. See copying.md for legal info.

from pymol import cmd
import shutil
import collections
from pathlib import Path
import sys
import warnings


def get_chain_bb(selection):
    """
    returns nested dictionary with format {object: {chain: list-of-bb-atoms}}
    """
    bb_name = "BB"
    bb_beads = {}

    # get list of objects in selection
    objects = cmd.get_names(selection=selection)

    for obj in objects:
        chains = cmd.get_chains(obj)
        bb_beads[obj] = {}
        for c in chains:
            # if chain is empty string, put it in the "*" bin
            if not c:
                c = "*"
            id_list = cmd.identify(f"{obj} and chain {c} and name {bb_name}")
            bb_beads[obj][c] = id_list
    return bb_beads


def get_gmx(gmx_path):
    """
    if gmx binary is not given, find it. If it can't be found, raise an exception
    """
    if not gmx_path:
        gmx_path = shutil.which('gmx')
    if not gmx_path:
        raise FileNotFoundError('no gromacs executable found.'
                                'Add it manually with gmx="PATH_TO_GMX"')
    return clean_path(gmx_path)


def update_recursive(base_dict, input_dict):
    """
    similar to builtin dict.update, but recursively updates all sub-dictionaries
    """
    for k, v in input_dict.items():
        if isinstance(v, collections.Mapping):
            base_dict[k] = update_recursive(base_dict.get(k, {}), v)
        else:
            base_dict[k] = v
    return base_dict


def clean_path(path):
    """
    cleans up paths and resolves ~ and symlinks
    """
    realpath = Path(path).expanduser().resolve()
    if not realpath.exists():
        #FIXME this is a shitty solution, we should deal with conditional include statements
        warnings.warn(realpath.name+" doesn't exist")
        return None
    return realpath
