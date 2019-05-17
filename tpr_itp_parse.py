from pymol import cmd
import os
import subprocess
import shutil
import re


def get_bonds_top(filepath):
    """
    Parses a topology file and returns a dictionary with as keys the moleculetypes
    and as values all the bonds and constraints for that moleculetype

    recursive call for #includes
    """
    # TODO: find a decent way to recognize elastic bonds. I probably involves using
    #       comments (ugly and unreliable) or distances (hard)
    bonds = {}

    # get top file directory for recursion
    dir = os.path.normpath(os.path.dirname(filepath))

    with open(filepath) as f:
        use_data = False
        new_molecule = False
        current_molecule = ''
        # line number added for easier debugging
        for ln, line in enumerate(f.readlines()):
            # remove comments
            line = line.partition(';')[0]
            # tokenize and join line
            tokens = line.split()
            line = ''.join(tokens)

            # recursive call for included topologies
            if line.startswith('#include'):
                include = tokens[1].strip("'").strip('"')
                include_bonds = get_bonds_top(f'{dir}/{include}')
                # update() replaces existing molecules with the same name, which is what we want
                bonds.update(include_bonds)

            # check if we got to a new header and if it's useful
            if line.startswith('['):
                if line == '[moleculetype]':
                    new_molecule = True
                    continue
                if line in ['[bonds]', '[constraints]']:
                    use_data = True
                    continue
                else:
                    use_data = False
                    continue

            # new molecule? Add it to the dictionary!
            if new_molecule and line:
                try:
                    current_molecule = tokens[0]
                    bonds[current_molecule] = []
                    new_molecule = False
                except:
                    raise ValueError(f'could not parse line {ln+1} of the file {filepath}.')

            # when in a useful section and line contains good data, store it
            if use_data and line:
                try:
                    bonds[current_molecule].append((tokens[0], tokens[1]))
                except:
                    raise ValueError(f'could not parse line {ln+1} of the file {filepath}.')
    return bonds


def get_bonds_tpr(filepath, gmx=None):
    """
    Uses gmx dump to parse a tpr file and returns a dictionary with as keys the
    moleculetypes and as values all the bonds and constraints for that moleculetype
    """
    # get gmx executable
    if not gmx:
        gmx = shutil.which('gmx')
    if not gmx:
        raise FileNotFoundError(f'no gromacs executable found. Add it manually with gmx="PATH_TO_GMX"')

    dump = subprocess.run(['gmx', 'dump', '-s', filepath], capture_output=True, text=True)
    lines = dump.stdout.split('\n')

    # utility function to search several regexes at once
    def search_all(target, *args):
        ret = []
        for re in args:
            match = re.search(target)
            ret.append(match)
        return ret

    # regex objects to parse the dumped data
    re_mol_def = re.compile(r'(moltype\s*=\s*)(\d*)(\s*")(\w*)')     # molecule definition
    re_atom = re.compile(r'(atom\[)(\d*)(\]\=\{name\=")(BB)')       # atom name definition
    re_mol = re.compile(r'(moltype\s\()(\d*)')                      # new molecule section
    re_bond = re.compile(r'(\(BONDS\)\s*)(\d*)(\s*)(\d*)')          # bond definition
    re_con = re.compile(r'(\(CONSTR\)\s*)(\d*)(\s*)(\d*)')          # constraint definition
    re_el = re.compile(r'(\(HARMONIC\)\s*)(\d*)(\s*)(\d*)')         # elastic network definition

    molecules = {}
    current_molecule = 'intermolecular'
    bb_beads = {current_molecule: []}
    bonds = {current_molecule: []}
    elastic = {current_molecule: []}

    for line in lines:
        # search all regexes in current line
        md, m, a, b, c, e = search_all(line, re_mol_def, re_mol, re_atom, re_bond, re_con, re_el)
        if md:
            # new moltype definitions
            mol_id, mol_name = md.group(2, 4)
            molecules[mol_id] = mol_name
        if m:
            # new section delimited by a known molecule type
            # initialize everything we need
            current_molecule = molecules[m.group(2)]
            bb_beads[current_molecule] = []
            bonds[current_molecule] = []
            elastic[current_molecule] = []
        if a:
            # store backbone atoms to later fix bond type in proteins
            bb_beads[current_molecule].append(a.group(2))
        if b:
            # bonds
            bonds[current_molecule].append(b.group(2, 4))
        if c:
            # constraints
            bonds[current_molecule].append(c.group(2, 4))
        if e:
            # elastic network
            elastic[current_molecule].append(e.group(2, 4))

    bonds_old, elastic_old = bonds, elastic
    # move short range elastic bonds from `bonds` to `elastic` (protein fix)
    for mol, b_list in bonds.items():
        for b in b_list:
            # check if both beads are in backbone list
            if b[0] in bb_beads[mol] and b[1] in bb_beads[mol]:
                # if they're not adjacent, move bond to elastic and remove from here
                if abs(bb_beads[mol].index(b[0]) - bb_beads[mol].index(b[1])) > 1:
                    elastic[mol].append(b)
                    bonds[mol].remove(b)

    return bonds, elastic


def find_molecules(molecules, selection='all'):
    """
    Given a list of molecule names, returns all the occurrences of those molecules
    in the loaded structure as a dict: {mol1: offset1, mol2...}
    """
    # look in the loaded stuff
    def get_mol_data(atom):
        print(atom.id)
    cmd.iterate(selection, get_mol_data())


def cg_bonds(bonds, elastic, selection='all'):
    """
    draws bonds and elastic network on a coarse grained structure
    """
    # TODO: store for later if someone wants to delete them
    for mol, b_list in bonds.items():
        for b in b_list:
            first_bead = int(b[0]) + 1
            second_bead = int(b[1]) + 1
            cmd.bond(f'id {first_bead}', f'id {second_bead}')
#    for mol, e_list in elastic.items():
#        for e in e_list:
#            first_bead = int(e[0]) + 1
#            second_bead = int(e[1]) + 1
#            cmd.bond(f'id {first_bead}', f'id {second_bead}')


cmd.extend('cg_bonds', cg_bonds)


bonds, elastic = get_bonds_tpr('em.tpr')
print(bonds, elastic)
cg_bonds(bonds, elastic)
#print(get_bonds_top('cg.top'))

