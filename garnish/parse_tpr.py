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

import subprocess
import re
import sys
from pymol import cmd

# local imports
from garnish.utils import get_gmx, clean_path


def parse_tpr(tpr_file, gmx=None):
    """
    parses the gmx dump of a tpr file and returns useful information on the system

    returns a dictionary describing the system
    """
    gmx = get_gmx(gmx)

    # dump tpr info in a string
    gmxdump = subprocess.run([gmx, "dump", "-s", clean_path(tpr_file)],
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
    info_section = True
    system = {
        'blocks': {},
        'topology': {},
    }

    # split the dump in lines
    for line in gmxdump.stdout.split('\n'):
        # parse for molecule information in the first section
        if info_section:
            for k, p in regexp_info.items():
                match = p.match(line)
                if match:
                    if k == 'molblock':
                        # save current molecule block id
                        curr_block_id = match.group(1)
                        system['blocks'][curr_block_id] = {}
                    if k == 'moltype':
                        # save molecule type id of current block and init topology for it
                        curr_mol_type = match.group(1)
                        system['blocks'][curr_block_id]['moltype'] = curr_mol_type
                        system['topology'][curr_mol_type] = {}
                    if k == 'molcount':
                        # save number of molecules in current block
                        n_molecules = int(match.group(1))
                        system['blocks'][curr_block_id]['n_molecules'] = n_molecules
                    if k == 'endinfo':
                        # stop parsing these patterns
                        info_section = False
        else:
            # look for a new molecule definition
            match = regexp_header.match(line)
            if match:
                # if molecule header was found, initialize some stuff
                # molecule type id
                curr_mol_type = match.group(1)
                # corresponding bond dictionary
                system['topology'][curr_mol_type]['connectivity'] = {
                    'bonds': [],
                    'constr': [],
                    'harmonic': []
                }
                # corresponding number of atoms
                system['topology'][curr_mol_type]['n_atoms'] = 0
                # corresponding backbone list
                system['topology'][curr_mol_type]['backbone'] = []
            # start looking for useful data in the line
            for key, p in regexp_data.items():
                match = p.match(line)
                if match:
                    if key == 'atomnames':
                        # update atom number for current molecule
                        system['topology'][curr_mol_type]['n_atoms'] += 1
                        # save backbone beads for later fix of short elastic bonds
                        at_nr = int(match.group(1))
                        at_name = match.group(2)
                        if at_name == "BB":
                            system['topology'][curr_mol_type]['backbone'].append(at_nr)
                    else:
                        # other matches are all connectivity information
                        bond_type = key
                        bond = tuple(int(b) for b in match.group(1, 2))
                        system['topology'][curr_mol_type]['connectivity'][bond_type].append(bond)
    return system


if __name__ == '__main__':
    args = sys.argv[1:]
    parse_tpr(*args)

if __name__ == 'pymol':
    cmd.extend('parse_tpr', parse_tpr)
