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
import sys
from pymol import cmd  # stored
import networkx as nx
import numpy as np
from pathlib import Path
import re
import subprocess
import shutil
from glob import glob


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

def parse_tpr(tpr_file, gmx=None):
    """
    parses the gmx dump of a tpr file and returns useful information on the system

    returns: - mol_blocks: number of molecules per type
             - mols_atom_n: number of atoms in each molecule type
             - mols_bonds: bonds dictionary
             - backbones: backbone beads
    """
    gmx = get_gmx(gmx)

    # dump tpr info in a string
    gmxdump = subprocess.run([gmx, "dump", "-s", str(tpr_file.absolute())],
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
    backbone = {}
    info_section = True

    system = {}

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
                    if k == 'moltype':
                        # save current molecule type id
                        curr_mol_type = match.group(1)
                        system[curr_mol_type] = {
                                'id': int(curr_block_id),
                                }

                    if k == 'molcount':
                        # create new entry in mol_blocks as {block_id: (moltype, molcount)}
                        n_molecules = int(match.group(1))
                        system[curr_mol_type]['n_molecules'] = n_molecules
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
                system[curr_mol_type]['connectivity'] = {
                    'bonds': [],
                    'constr': [],
                    'harmonic': []
                }
                # corresponding number of atoms
                system[curr_mol_type]['n_atoms'] = 0
                # corresponding backbone list
                system[curr_mol_type]['backbone'] = []

            for key, p in regexp_data.items():
                match = p.match(line)
                if match:
                    if key == 'atomnames':
                        system[curr_mol_type]['n_atoms'] += 1
                        # save backbone beads for later fix of short elastic bonds
                        at_nr = int(match.group(1))
                        at_name = match.group(2)
                        if at_name == "BB":
                            system[curr_mol_type]['backbone'].append(at_nr)
                    else:
                        bond_type = key
                        bond = tuple(int(b) for b in match.group(1, 2))
                        system[curr_mol_type]['connectivity'][bond_type].append(bond)
    return system


def parse_top(top_file):
    """
    parses a topology file and returns useful information on the system

    returns: - mol_blocks: number of molecules per type
             - mols_atom_n: number of atoms in each molecule type
             - mols_bonds: bonds dictionary
             - backbones: backbone beads
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
    mol_blocks = []
    mols_atom_n = {}
    mols_bonds = {}
    backbone = {}
    section = False

    # we need to shift all the atom ids, because topologies start from 1 and not 0
    id_fix = -1

    # Matthijs: Will try to make parse_top work internally with a dictionary
    # After all the includes are done, we take the molecules and everything as defined by the moleculetypes in the top file

    # Initialize an empty dictionary
    # keys between in "[]" are variable
    # The dictionary structure will be:
    #
    # included = {
    #    [moleculetype]: {
    #       id:             int(id nr to check with tpr),
    #       n_molecules:    int(#molecules),
    #       n_atoms:        int(#atoms),
    #       connectivity: {
    #           bonds:      list[ ( int(at#), int(at#) ) ],
    #           constr:     list[ ( int(at#), int(at#) ) ],
    #           harmonics:  list[ ( int(at#), int(at#) ) ],
    #       },
    #       backbone: list[ ( int(at#), int(at#) ) ] 
    #    }
    # }
    #           
    #           

    # read the file as lines and parse data
    with open(top_file, 'r') as f:
        included = {} 
        ID = 0 
        system = {}
        for line in f.readlines():
            # look for included topologies
            match = regexp_include.match(line)
            if match:
                include_file = match.group(1)
                include_path = Path(str(top_file.parent)+"/"+include_file)
                # recursive call for included topologies
                # add all the data we found to the main containers
                included.update(parse_top(include_path))

            # look for a header in the form of `[something]`
            match = regexp_header.match(line)
            if match:
                section = match.group(1)
            if section == 'molecules':
                # save amount of molecules for each block
                match = regexp_info.match(line)
                if match:
                    # Add this molecule to the system
                    moleculetype = match.group(1)
                    n_molecules = int(match.group(2))
                    system[moleculetype] = included[moleculetype]
                    system[moleculetype]['n_molecules'] = n_molecules
                    system[moleculetype]['id'] = ID
                    ID += 1
            if section == 'moleculetype':
                match = regexp_moltype.match(line)
                if match:
                    # if molecule header was found, initialize some stuff
                    # molecule type id
                    curr_mol_type = match.group(1)
                    # main molecule dictionary
                    included[curr_mol_type] = {
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
                    included[curr_mol_type]['n_atoms'] += 1
                    # save backbone beads for later fix of short elastic bonds
                    if atom_name == "BB":
                        included[curr_mol_type]['backbone'].append(atom_id + id_fix)
            if section == 'bonds':
                match = regexp_bond.match(line)
                if match:
                    # save bond, fixing the atom ids to match those of pymol
                    bond = tuple(int(b) + id_fix for b in match.group(1, 2))
                    included[curr_mol_type]['connectivity']['bonds'].append(bond)
            if section == 'constraints':
                match = regexp_constr.match(line)
                if match:
                    # save constraint, fixing the atom ids to match those of pymol
                    constr = tuple(int(b) + id_fix for b in match.group(1, 2))
                    included[curr_mol_type]['connectivity']['constr'].append(constr)

    if not system:
        return included
    else:
        return system

def make_graphs(system):
    """
    uses data gathered from file parsing to correctly identify bonds for each molecule

    returns graph representations of all the bonds in the system
    """

    # Iterate over molecules, fix elastic bonds where necessary
    for key in system.keys():
        molecule = system[key]
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
            if all(atom in backbone for atom in bond):
                at_a, at_b = bond
                at_a, at_b = int(at_a), int(at_b)
                atom1_idx = backbone.index(at_a)
                atom2_idx = backbone.index(at_b)
                # if it's also between two non-adjacent beads, move it to `elastic`
                if abs(atom1_idx - atom2_idx) > 1:
                    bond_type = tmp_harmonics
            bond_type.append(bond)
        molecule['connectivity']['bonds'] = tmp_bonds
        molecule['connectivity']['harmonic'] = tmp_harmonics

    # Convert the lists of bonds to networkx graphs
    bond_graphs = {}

    offset = 0
    for moleculetype, molecule in system.items():
        n_mol = molecule['n_molecules']
        n_at = molecule['n_atoms']
        connectivity = molecule['connectivity']
        # Need to have numpy arrays to apply offset
        connectivity = { btype: np.array(bonds) for btype, bonds in connectivity.items() }

        for i in range(n_mol):
            key = str(molecule['id'])+f"_{i}"
            bond_graphs[key] = {}
            for btype, bonds in connectivity.items():
                g = nx.Graph()
                g.add_edges_from(bonds + offset)
                bond_graphs[key][btype] = g
            offset += n_at

    return bond_graphs


def rel_atom(selection):
    # Make a dict of all the atoms (to get effective relative atom numbering)
    rel_atom_dict = {}
    atoms = cmd.get_model(selection)
    for i, at in enumerate(atoms.atom):
        rel_atom_dict[i] = at.index
    return rel_atom_dict


def cg_bonds(file=None, selection='all'):
    """
    Allow a cg structure to be visualized in pymol like an atomistic structure.

    Usage: cg_bonds [selection], [tpr_file]

    selection   : any selection to act upon (default: all)
    tpr_file : a .tpr file to extract bond information from (default: None)

    Without a tpr file, this function only adds bonds between the backbone beads so they can be
    nicely visualized using line or stick representation.
    A tpr file provides topology information that can be used to draw side chain and elastic bonds.
    """
    tpr_file = None
    top_file = None
    aa_template = None
    maybe_file = Path(str(file))

    if file == None:
        pass
    elif maybe_file.is_file():
        if maybe_file.suffix == ".tpr":
            tpr_file = maybe_file
        elif maybe_file.suffix == ".top" or maybe_file.suffix == ".itp":
            top_file = maybe_file
        elif maybe_file.suffix == ".pdb":
            aa_template = maybe_file
        else:
            raise TypeError(f'"{maybe_file.suffix}" is not a supported format.')
    else:
        raise FileNotFoundError(f'{maybe_file} does not exist.')

    # Order might be important
    cmd.set("retain_order", 1)  # TODO: is it really though?
    # Yes, pymol otherwise sorts the atoms, if you then save the file you get a different result.

    # Fix the view nicely
    cmd.hide("everything", selection)
    cmd.show_as("lines", selection)
    cmd.util.cbc(selection)

    # Draw all the bonds based on tpr or top file
    if tpr_file or top_file:
        # warn at the end if some atoms in the tpr are missing from the structure
        warn = False
        # Get bond graphs
        if tpr_file:
            parsed_data = parse_tpr(tpr_file)
        elif top_file:
            parsed_data = parse_top(top_file)
        bond_graphs = make_graphs(parsed_data)
        selection_objects = cmd.get_object_list(selection)
        for obj in selection_objects:
            elastics_obj = selection+"_elastics"
            # Create dummy object to draw elastic bonds in
            cmd.copy(elastics_obj, obj)
            # Make a dict of all the atoms (to get effective relative atom numbering)
            rel_atom_selection = rel_atom(obj)
            # Draw all the bonds
            for btype in ['bonds', 'constr']:
                for _, bonds in bond_graphs.items():
                    for a, b in bonds[btype].edges:
                        try:
                            a = rel_atom_selection[a]
                            b = rel_atom_selection[b]
                            cmd.bond(f"({obj} and ID {a})", f"({obj} and ID {b})")
                        except KeyError:
                            warn = True
                # Get relative atoms for elastics object
            rel_atom_elastics = rel_atom(elastics_obj)
            atoms = cmd.get_model(elastics_obj)
            for i, at in enumerate(atoms.atom):
                rel_atom_elastics[i] = at.index
            # Draw elastic network
            for _, bonds in bond_graphs.items():
                for a, b in bonds['harmonic'].edges:
                    try:
                        a = rel_atom_elastics[a]
                        b = rel_atom_elastics[b]
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

cmd.extend('cg_bonds', cg_bonds)

# tab completion for the cg_bonds command
useful_file_sc = lambda: cmd.Shortcut(
     glob('*/') + glob('*.tpr') + glob('*.top') + glob('*.itp')  # + glob('*.pdb')
)
cmd.auto_arg[0]['cg_bonds'] = [useful_file_sc, 'input file', ', ']
# here object_sc is more informative than selection_sc
cmd.auto_arg[1]['cg_bonds'] = [cmd.object_sc, 'selection', '']
