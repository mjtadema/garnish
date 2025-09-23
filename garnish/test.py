import pytest
from pymol import cmd
import pkg_resources
from pathlib import Path

from .extensions import garnish
from . import Parser


@pytest.fixture
def structure():
    return Path(pkg_resources.resource_filename(__name__, "data/4TSY_cg_membrane.gro"))


@pytest.fixture
def topology():
    return Path(pkg_resources.resource_filename(__name__, "data/membrane.top"))


@pytest.fixture
def init_pymol(structure):
    cmd.load(structure)


@pytest.fixture
def system(topology):
    return Parser(topology).run()


def test_topology(system):
    molecules = set([mol for mol, _ in system.topology])
    for mol in system:
        assert mol in molecules


def test_garnish(init_pymol, topology):
    garnish(topology)
