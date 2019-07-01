# Copyright 2019-2019 the garnish authors. See copying.md for legal info.

import networkx as nx

class System():
    def __init__(self, args*, kwargs**):
        """
DESCRIPTION:
    The system class holds a collection of molecules, represented by graphs.
    Additionally it provides several helper methods to manipulate and get information about the system.
        """
        pass

class Bead():
    def __init__(self, args*, kwargs**):
        """
DESCRIPTION:
    Holds information about beads/atoms.
    Can hold keys for looking up upstream information in the parent "System" instance.
        """
        self.atom_id = args.pop()

    def __repr__(self):
        return self.atom_id

def test_system(test_data):
    for 

test_data = [
        (0,1),
        (1,2),
        (2,3),
        (3,4)
        ]
test_atom = 1

if __name__ == "__main__":
    #test_system(test_data)
