import pycg_bonds as target
from pymol import cmd, stored
import os

def test_get_chains(selection):
    chains = target.get_chains(selection)
    if len(chains) < 1:
        raise ValueError("got no chains")
    elif selection == "cg_oligomer" and len(chains) <= 1:
        raise ValueError("not enough chains for oligomer")

def run_tests():
    # Fail count
    global fail_count
    fails = []

    # Move to testdir
    os.chdir("test")

    # Load test_structs
    cmd.load("cg_monomer.pdb")
    cmd.load("cg_oligomer.pdb")
    to_test = [ "cg_monomer", "cg_oligomer"]

    # Run all the tests, add exception to fails
    try:
        for selection in to_test:
            test_get_chains(selection)
    except Exception as e:
        fails.append(e)
    
    # Report results
    print(str(len(fails))+" tests failed")
    if fails:
        for e in fails:
            print(e)
        exit(1)

run_tests()
