#!/usr/bin/env python3

from pycg_bonds import parse_tpr
from pycg_bonds import parse_top
from pycg_bonds import make_graphs
from pathlib import Path
import networkx as nx

def test_parse():
    """
    Crude integration test for parsing logic.
    Since we have two different parsing functions, we can compare results.
    The resulting graphs from these should be identical.
    If they are not identical the test will raise an exception
    """
    mismatches = []
    tpr_file = Path("lysotest/lyso.tpr")
    top_file = Path("lysotest/lyso.top")
    #mol_blocks, mols_atom_n, mols_bonds, backbone
    tpr_output = parse_tpr(tpr_file)
    top_output = parse_top(top_file)
    tpr_graphs = make_graphs(tpr_output)
    top_graphs = make_graphs(top_output)

    # Dictionaries should be sorted in python3 but probably
    # we will have to come up with something more robust...
    for tpr, top in zip(tpr_graphs.values(), top_graphs.values()):
        for bt in tpr.keys():
            if not all(bond in tpr[bt].edges() for bond in top[bt].edges()):
                 raise Exception("Molecules are not equal!")

test_parse()
print("Tests passed successfully")
