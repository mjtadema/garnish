# Copyright 2019-2019 the garnish authors. See copying.md for legal info.

from pathlib import Path
import warnings

# local imports
from .system import System, Molecule, Atom


class Parser:
    """
    Parser object to handle parsing a topology
    """
    def __init__(self, top_file):
        self.top_file = Path(top_file)
        self.system = System()
        self._active_molecule = None
        self._parsers = {
            "bonds": self._bonds,
            "moleculetype": self._moleculetype,
            "atoms": self._atoms,
            "molecules": self._molecules,
            "constraints": self._constraints
        }

    def run(self):
        """
        Initialize the toplevel parser
        :return: System
        """
        self.top(self.top_file)
        return self.system

    def top(self, top_file):
        """
        Top level parser, can be run recursively for #include files.
        Looks for headers containing topology information.
        Calls parsers for each relevant line.
        The parsers handle storing the data.
        """
        p = Path(top_file)
        if not p.exists():
            warnings.warn(f"{p.name} not found, skipping...")
            return

        header = "none"
        with open(top_file, 'r') as top_handle:
            for line in top_handle:
                # Strip whitespace always
                line = line.strip()

                if line.startswith("["):
                    header = line.strip("[ ]")
                    continue

                # Handle include files
                if line.startswith("#include"):
                    _, itp_file = line.strip().split()
                    itp_file = itp_file.strip("\"")
                    itp_file = Path(itp_file)
                    if not itp_file.name.startswith("/"):
                        itp_file = self.top_file.parent/itp_file
                    self.top(itp_file)
                    continue

                # Get rid of any comments
                line = line.rsplit(";", 1)[0]
                # ...and any define statements
                line = line.rsplit("#", 1)[0]

                # Skip empty lines
                if line == "":
                    continue

                try:
                    self._parsers[header](line)
                except KeyError:
                    # Skip unneeded headers
                    pass

    def _moleculetype(self, line):
        name, nrexcl = line.strip().split()
        molecule = Molecule(name=name)
        self.system.extend(molecule)
        self.active_molecule = molecule

    def _atoms(self, line):
        idx, _type, resi, resn, name, *_ = line.strip().split()
        idx = int(idx)
        resi = int(resi)
        atom = Atom(idx, _type, resi, resn, name)
        self.active_molecule.extend(atom)

    def _bonds(self, line):
        _from, to, bondtype, distance, fc = line.strip().split()

        _from = int(_from)
        to = int(to)
        distance = float(distance)
        fc = float(fc)

        if distance <= 0.47 and fc >= 1250:
            # Short distance, high force constant bonds are considered covalent
            self.active_molecule[_from].add_bond(to)
        else:
            # Long distance, low force constant bonds are considered elastic
            # regardless of bondtype
            self.active_molecule[_from].add_elastic(to)

    def _constraints(self, line):
        # Constraints are considered bonds
        _from, to, *_ = line.strip().split()
        _from = int(_from)
        to = int(to)
        self.active_molecule[_from].add_bond(to)

    def _molecules(self, line):
        name, n = line.strip().split()
        n = int(n)
        self.system.order(name, n)
