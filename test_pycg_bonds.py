import pycg_bonds as target
from pymol import cmd, stored, CmdException
import os
from pdb import set_trace as bp

def test_get_chains():
    cmd.reinitialize()
    for selection in ["cg_monomer","cg_oligomer"]:
        cmd.load(selection+".pdb")
        chains = target.get_chains(selection)
        if len(chains) < 1:
            raise ValueError("got no chains")
        elif selection == "cg_oligomer" and len(chains) <= 1:
            raise ValueError("not enough chains for oligomer")

def test_get_chain_bb():
    cmd.reinitialize()
    selection = "cg_oligomer_nochains"
    cmd.load(selection+".pdb")
    chains = target.get_chains(selection)
    bp()
    target.get_chain_bb(selection, chains)

test_list = [
        test_get_chains,
        test_get_chain_bb
        ]

def run_tests(test_functions):
    # Fail count
    global fail_count
    fails = []

    # Move to testdir
    os.chdir("test")

    # Run all the tests, add exception to fails
    try:
        for test_func in test_functions:
            test_func()
    except Exception as e:
        fails.append(e)
    except CmdException as e:
        fails.append(e)
    
    # Report results
    print(str(len(fails))+" tests failed")
    if fails:
        for e in fails:
            print(e)
        exit(1)

run_tests(test_list)
