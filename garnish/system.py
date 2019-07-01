import networkx as nx
import numpy as np


class System:
    def __init__(self, sys_dict):
        self.sys_dict = sys_dict
        self.graph = nx.Graph()
        self.data = {}

        # fix wrong elastic bonds
        self.fix_elastics()

        self.make_graph()
        #self.make_data(sys_dict)

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
        for block_id, block in sys_dict['blocks'].items():
            moltype = block['moltype']
            n_mol = block['n_molecules']
            n_at = sys_dict['topology'][moltype]['n_atoms']
            connectivity = sys_dict['topology'][moltype]['connectivity']
            # transform bonds in numpy arrays to easily apply offset
            connectivity = {btype: np.array(bonds) for btype, bonds in connectivity.items()}

            # repeat for each occurrence of molecule in this block
            for i in range(n_mol):
                for btype, bonds in connectivity.items():
                    # create a graph based on connectivity and offset it to match atom numbers
                    self.graph.add_edges_from(bonds + offset)
                # shift offset by how many atoms this molecule has
                offset += n_at

    def make_data(self):
        pass


if __name__ == '__main__':
    from garnish.parse_top import parse_top
    import sys
    from pprint import pprint as print
    sys_dict = parse_top(sys.argv[1])
    s = System(sys_dict)

    print(s.graph.size())
