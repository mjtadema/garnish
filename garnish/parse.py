# Copyright 2019-2019 the garnish authors. See copying.md for legal info.

import re
import sys
import subprocess
from pymol import cmd

# local imports
from .utils import get_gmx, update_recursive, clean_path


def parse_top(top_file):
    """
    parses a topology file and returns useful information on the system

    returns a dictionary describing the system
    """
    top_file = clean_path(top_file)

    # define regex patterns for later use
    regexp_include = re.compile('^#include\s+\"(.*?)\"')

    regexp_header = re.compile('^\s*\[\s*(\w+)\s*\]')

    regexp_info = re.compile('^\s*(\S+)\s+(\d+)')
    regexp_moltype = re.compile('^\s*(\S+)\s+(\d+)')
    regexp_atom = re.compile('^\s*(\d+)\s+(\S+)\s+\d+\s+\S+\s+(\S+)')
    regexp_bond = re.compile('^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d\.]+)')
    #regexp_bond = re.compile('^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d\.]+)\s+([\d\.]+)')
    regexp_constr = re.compile('^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d\.]+)')

    # initialize some stuff
    block_id = 0
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
            if match := regexp_include.match(line):
                include_file = match.group(1)
                include_path = top_file.parent/include_file
                # recursive call for included topologies and recursive update of system dictionary
                update_recursive(system, parse_top(include_path))
            # look for a header in the form of `[something]`
            if match := regexp_header.match(line):
                section = match.group(1)
            if section == 'molecules':
                if match := regexp_info.match(line):
                    # save molecule info
                    curr_mol_type = match.group(1)
                    n_molecules = int(match.group(2))
                    curr_block_id = str(block_id)
                    # add molecule info to system
                    system['blocks'][curr_block_id] = {}
                    system['blocks'][curr_block_id]['moltype'] = curr_mol_type
                    system['blocks'][curr_block_id]['n_molecules'] = n_molecules
                    # update current block id. Needed because, differently from tpr,
                    # top does not define an id for blocks
                    block_id += 1
            if section == 'moleculetype':
                if match := regexp_moltype.match(line):
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
                        'atomtypes': {},
                        'backbone': []
                    }
            if section == 'atoms':
                if match := regexp_atom.match(line):
                    atom_id = int(match.group(1))
                    atom_type = match.group(2)
                    atom_name = match.group(3)
                    # increment atom count for current molecule type
                    system['topology'][curr_mol_type]['n_atoms'] += 1
                    # also save all the atom types in system
                    system['topology'][curr_mol_type]['atomtypes'][atom_id+id_fix] = atom_type
                    # save backbone beads for later fix of short elastic bonds
                    if atom_name == "BB":
                        system['topology'][curr_mol_type]['backbone'].append(atom_id + id_fix)
            if section == 'bonds':
                if match := regexp_bond.match(line):
                    # save bond, fixing the atom ids to match those of pymol
                    bond = tuple(int(b) + id_fix for b in match.group(1, 2))
                    system['topology'][curr_mol_type]['connectivity']['bonds'].append(bond)
            if section == 'constraints':
                if match := regexp_constr.match(line):
                    # save constraint, fixing the atom ids to match those of pymol
                    constr = tuple(int(b) + id_fix for b in match.group(1, 2))
                    system['topology'][curr_mol_type]['connectivity']['constr'].append(constr)

    return system


def parse_tpr(tpr_file, gmx=None):
    """
    parses the gmx dump of a tpr file and returns useful information on the system

    returns a dictionary describing the system
    """
    tpr_file = clean_path(tpr_file)
    gmx = get_gmx(gmx)

    # dump tpr info in a string
    gmxdump = subprocess.run([gmx, 'dump', '-s', clean_path(tpr_file)],
                             capture_output=True, text=True)

    # define regex patterns for later use
    regexp_info = {
        'molblock': re.compile('^\s+molblock\s+\((\d+)\)'),
        'moltype': re.compile('^\s+moltype\s+=\s+(\d+)'),
        'molcount': re.compile('^\s+#molecules\s+=\s+(\d+)'),
        'endinfo': re.compile('^\s+ffparams')
    }

    regexp_header = re.compile('^\s+moltype\s+\((\d+)\):')

    regexp_data = {
        'atomnames': re.compile('^\s+atom\[(\d+)\]=\{name=\"(\S+?)'),
        'atomtypes': re.compile('^\s+type\[(\d+)\]=\{name=\"(\S+)\",'),
        'bonds': re.compile('^\s+\d+\s\w+=\d+\s\(BONDS\)\s+(\d+)\s+(\d+)'),
        'constr': re.compile('^\s+\d+\s\w+=\d+\s\(CONSTR\)\s+(\d+)\s+(\d+)'),
        'harmonic': re.compile('^\s+\d+\s\w+=\d+\s\(HARMONIC\)\s+(\d+)\s+(\d+)')
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
                if match := p.match(line):
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
            if match := regexp_header.match(line):
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
                # corresponding atom types dictionary
                system['topology'][curr_mol_type]['atomtypes'] = {}
                # corresponding backbone list
                system['topology'][curr_mol_type]['backbone'] = []
            # start looking for useful data in the line
            for key, p in regexp_data.items():
                if match := p.match(line):
                    if key == 'atomnames':
                        # update atom number for current molecule
                        system['topology'][curr_mol_type]['n_atoms'] += 1
                        # save backbone beads for later fix of short elastic bonds
                        at_nr = int(match.group(1))
                        at_name = match.group(2)
                        if at_name == "BB":
                            system['topology'][curr_mol_type]['backbone'].append(at_nr)
                    elif key == 'atomtypes':
                        # also save all the atom types in system
                        at_nr = int(match.group(1))
                        at_type = match.group(2)
                        system['topology'][curr_mol_type]['atomtypes'][at_nr] = at_type
                    else:
                        # other matches are all connectivity information
                        bond_type = key
                        bond = tuple(int(b) for b in match.group(1, 2))
                        system['topology'][curr_mol_type]['connectivity'][bond_type].append(bond)
    return system


# TODO
#def parse_pdb(pdb_file):
#    return None
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


def parse(file=None, gmx=None):
    """
    dispatches parsing duty to the appropriate function depending on the file extension

    returns a dictionary describing the system
    """
    maybe_file = clean_path(file)
    ext = maybe_file.suffix
    # return a differently parsed system depending on the provided file
    if ext == ".tpr":
        return parse_tpr(maybe_file, gmx)
    elif ext in (".top", ".itp"):
        return parse_top(maybe_file)
    else:
        raise TypeError(f'"{maybe_file.suffix}" is not a supported format.')


if __name__ == '__main__':
    args = sys.argv[1:]
    parse(*args)

if __name__ == 'pymol':
    cmd.extend('parse', parse)
