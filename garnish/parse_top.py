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

import re
import sys
import os
from pymol import cmd

# local imports
from garnish.utils import update_recursive, clean_path


def parse_top(top_file):
    """
    parses a topology file and returns useful information on the system

    returns a dictionary describing the system
    """
    # define regex patterns for later use
    regexp_include = re.compile('^#include\s+\"(.*?)\"')

    regexp_header = re.compile('^\s*\[\s*(\w+)\s*\]')

    regexp_info = re.compile('^\s*(\w+)\s+(\d+)')
    regexp_moltype = re.compile('^\s*(\S+)\s+(\d+)')
    regexp_atom = re.compile('^\s*(\d+)\s+\w+\s+\d+\s+\w+\s+(\w+)')
    regexp_bond = re.compile('^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d\.]+)')
    #regexp_bond = re.compile('^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d\.]+)\s+([\d\.]+)')
    regexp_constr = re.compile('^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d\.]+)')

    # initialize some stuff
    ID = 0
    system = {
        'blocks': {},
        'topology': {},
    }
    section = False
    # we need to shift all the atom ids, because topologies start from 1 and not 0
    id_fix = -1

    # read the file as lines and parse data
    with open(top_file, 'r') as f:
        for line in f.readlines():
            # look for included topologies
            match = regexp_include.match(line)
            if match:
                include_file = match.group(1)
                include_path = clean_path(os.path.join(os.path.dirname(top_file), include_file))
                # recursive call for included topologies and recursive update of system dictionary
                update_recursive(system, parse_top(include_path))
            # look for a header in the form of `[something]`
            match = regexp_header.match(line)
            if match:
                section = match.group(1)
            if section == 'molecules':
                match = regexp_info.match(line)
                if match:
                    # save molecule info
                    curr_mol_type = match.group(1)
                    n_molecules = int(match.group(2))
                    curr_block_id = str(ID)
                    # add molecule info to system
                    system['blocks'][curr_block_id] = {}
                    system['blocks'][curr_block_id]['moltype'] = curr_mol_type
                    system['blocks'][curr_block_id]['n_molecules'] = n_molecules
                    # update current block id. Needed because, differently from tpr,
                    # top does not define an id for blocks
                    ID += 1
            if section == 'moleculetype':
                match = regexp_moltype.match(line)
                if match:
                    # if molecule header was found, initialize its dictionary
                    # molecule type id
                    curr_mol_type = match.group(1)
                    # main molecule dictionary
                    system['topology'][curr_mol_type] = {
                        'connectivity': {
                            'bonds': [],
                            'constr': [],
                            'harmonic': []
                        },
                        'n_atoms': 0,
                        'backbone': []
                    }
            if section == 'atoms':
                match = regexp_atom.match(line)
                if match:
                    atom_id = int(match.group(1))
                    atom_name = match.group(2)
                    # increment atom count for current molecule type
                    system['topology'][curr_mol_type]['n_atoms'] += 1
                    # save backbone beads for later fix of short elastic bonds
                    if atom_name == "BB":
                        system['topology'][curr_mol_type]['backbone'].append(atom_id + id_fix)
            if section == 'bonds':
                match = regexp_bond.match(line)
                if match:
                    # save bond, fixing the atom ids to match those of pymol
                    bond = tuple(int(b) + id_fix for b in match.group(1, 2))
                    system['topology'][curr_mol_type]['connectivity']['bonds'].append(bond)
            if section == 'constraints':
                match = regexp_constr.match(line)
                if match:
                    # save constraint, fixing the atom ids to match those of pymol
                    constr = tuple(int(b) + id_fix for b in match.group(1, 2))
                    system['topology'][curr_mol_type]['connectivity']['constr'].append(constr)

    return system


if __name__ == '__main__':
    args = sys.argv[1:]
    parse_top(*args)

if __name__ == 'pymol':
    cmd.extend('parse_top', parse_top)
