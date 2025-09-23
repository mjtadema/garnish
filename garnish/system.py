# Copyright 2019-2019 the garnish authors. See copying.md for legal info.
import warnings


class Parent:
    """
    Base Parent class
    Containing children
    """
    def __init__(self):
        self.children = {}

    def __iter__(self):
        for child in self.children.values():
            yield child

    def __getitem__(self, item):
        return self.children[item]

    def __len__(self):
        return len(self.children.values())


class System(Parent):
    """
    System holds molecule definitions
    """
    def __init__(self):
        super().__init__()
        self.topology = []

    def __str__(self):
        return f"System with {len(self.children)} molecules"

    def extend(self, child):
        self.children[repr(child)] = child

    @property
    def bonds(self):
        return [list(molecule.bonds) for molecule in self]

    @property
    def elastics(self):
        return [list(molecule.elastics) for molecule in self]

    @property
    def molecules(self):
        return list(self.children.values())

    def order(self, name, n):
        try:
            molecule = self[name]
            self.topology.append((molecule, n))
        except KeyError:
            warnings.warn(f"molecule {name} has no definition, skipping...")


class Molecule(Parent):
    """
    Molecule holds atoms
    """
    def __init__(self, name):
        super().__init__()
        self.name = name

    def __repr__(self):
        return self.name

    def extend(self, child):
        self.children[int(child)] = child

    @property
    def bonds(self):
        for atom in self:
            _from = atom.idx
            for to in atom.bonds:
                yield _from, to

    @property
    def elastics(self):
        for atom in self:
            _from = atom.idx
            for to in atom.elastics:
                yield _from, to


class Atom:
    """
    Atoms keep track of bonded information.
    Each atom knows its own idx and holds indices of connected atoms
    """
    def __init__(self, idx, _type, *_):
        self.idx = idx
        self.bonds = []
        self.elastics = []

    def __int__(self):
        return self.idx

    def add_bond(self, to):
        self.bonds.append(to)

    def add_elastic(self, to):
        self.elastics.append(to)
