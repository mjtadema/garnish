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

from pymol import cmd
import shutil
import collections
import os
import sys


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
    resolves path variables, `~`, symlinks and returns a clean absolute path
    """
    return os.path.realpath(os.path.expanduser(os.path.expandvars(path)))

def extension(loading_func):
    """
    Decorator for pymol extension functions.
    Will only return the loading function if called by pymol, else returns an empty function so no errors are raised.
    These functions can then be called in __init__.py to extend pymol functions to pymol.
    """
    try:
        # check if module was called by pymol, if yes run skewer.main
        main_initfile = sys.modules['__main__'].__file__
    
        # get the name of the parent module
        main_modulename = os.path.basename(clean_path(os.path.join(main_initfile, os.pardir)))
    
    except AttributeError:
        # importing from an interpreter like ipython raises this error
        main_modulename = None
    
    
    if main_modulename == 'pymol':
        return loading_func
    else:
        # Just return a passing lambda doing nothing, to avoid errors further downstream
        return lambda: None
