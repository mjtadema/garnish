# Copyright 2019-2019 the garnish authors. See copying.md for legal info.

import networkx as nx
import numpy as np
from pymol import cmd


class System:
    def __init__(self, sys_dict):
        self.sys_dict = sys_dict
        self.graph = nx.Graph()

        # fix wrong elastic bonds
        self.fix_elastics()

        self.make_graph()

    def fix_elastics(self):
        # Iterate over molecules, fix elastic bonds where necessary
        for key in self.sys_dict['topology'].keys():
            molecule = self.sys_dict['topology'][key]
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
                # if both atoms in a bond are labeled as backbone, go deeper
                if all(atom in backbone for atom in bond):
                    at_a, at_b = bond
                    at_a, at_b = int(at_a), int(at_b)
                    atom1_idx = backbone.index(at_a)
                    atom2_idx = backbone.index(at_b)
                    # if the bond is between two non-adjacent backbone beads, move it to `elastic`
                    if abs(atom1_idx - atom2_idx) > 1:
                        bond_type = tmp_harmonics
                bond_type.append(bond)
            # update old dictionaries with new info
            molecule['connectivity']['bonds'] = tmp_bonds
            molecule['connectivity']['harmonic'] = tmp_harmonics

    def make_graph(self):
        # unfold the topology to create a graph with atoms as nodes
        offset = 0
        for block_id, block in self.sys_dict['blocks'].items():
            moltype = block['moltype']
            n_mol = block['n_molecules']

            n_at = self.sys_dict['topology'][moltype]['n_atoms']
            a_types = self.sys_dict['topology'][moltype]['atomtypes']

            connectivity = self.sys_dict['topology'][moltype]['connectivity']
            # transform bonds in numpy arrays to easily apply offset
            connectivity = {btype: np.array(bonds) for btype, bonds in connectivity.items()}

            # repeat for each occurrence of molecule in this block
            for i in range(n_mol):
                # create nodes and add atom data to each
                # TODO: there should be a way to vectorize this, I feel, but I can't figure it out
                # loop through all the atoms in the molecule
                for local_atom_id in range(n_at):
                    atom_id = local_atom_id + offset
                    at_type = a_types[local_atom_id]
                    self.graph.add_node(atom_id, moltype=moltype, block=block_id, atomtype=at_type)
                # add edges to graph
                for btype, bonds in connectivity.items():
                    # create a graph based on connectivity and offset it to match atom numbers
                    self.graph.add_edges_from(bonds + offset, type=btype)

                # shift offset by how many atoms this molecule has
                offset += n_at

    def draw_bonds(self, selection):
        warn = False
        selection_objects = cmd.get_object_list(selection)
        for obj in selection_objects:
            elastics_obj = obj+"_elastics"
            # Create dummy object to draw elastic bonds in
            cmd.copy(elastics_obj, obj)
            for a, b, data in self.graph.edges(data=True):
                # draw non-elastic bonds
                if data['type'] in ['bonds', 'constr']:
                    try:
                        a += 1
                        b += 1
                        try:
                            cmd.add_bond(obj, a, b)
                        except AttributeError:
                            cmd.bond(f"({obj} and ID {a})", f"({obj} and ID {b})")
                    except KeyError:
                        warn = True
                # Draw elastic network
                if data['type'] == 'harmonic':
                    try:
                        a += 1
                        b += 1
                        try:
                            cmd.add_bond(elastics_obj, a, b)
                        except AttributeError:
                            cmd.bond(f"({elastics_obj} and ID {a})", f"({elastics_obj} and ID {b})")
                    except KeyError:
                        warn = True
            cmd.color("orange", elastics_obj)
        # warn about missing atoms if needed.
        if warn:
            print('WARNING: some atoms present in the tpr file were not found in the loaded '
                  'structure.\n Bonds containing those atoms were not drawn.')

    def transfer_attributes(self, selection):
        data = self.graph.nodes(data=True)
        # create a namespace to feed to pymol
        tmp_namespace = {'data': data}
        cmd.alter(selection=f'{selection}', space=tmp_namespace,
                  expression=f'elem=data[ID]["atomtype"]')
