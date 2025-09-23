# Copyright 2019-2019 the garnish authors. See copying.md for legal info.

from pathlib import Path
import warnings
import logging

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
        self.variables = {}

    def resolve(self, string : str):
        """
        Resolve variables in the string and (safely) evaluate the expression
        :param string:
        :return: evaluated expression
        """
        code = compile(string, "<string>", "eval")
        for name in code.co_names:
            if name not in self.variables:
                # Have to do this to not allow malicious code
                raise NameError("Use of undefined names is not allowed")
        # An empty namespace and empty __builtins__ are passed for safety
        return eval(code, {"__builtins__": {}}, self.variables)

    def run(self):
        """
        Initialize the toplevel parser
        :return: System
        """
        self.top(self.top_file)
        logging.debug("Found following molecules in topology")
        logging.debug(self.system.topology)
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

                # Get rid of any comments
                line = line.split(";")[0]

                if line.startswith("["):
                    header = line.strip("[ ]")
                    continue

                # Handle include files
                if line.startswith("#include"):
                    self.parse_include(line)
                    continue

                if line.startswith("#define"):
                    self.parse_define(line)
                    continue

                # ignore any # directives we don't catch specifically
                line = line.split("#")[0]

                # Skip empty lines after filtering
                if line == "":
                    continue

                try:
                    getattr(self, "parse_"+header)(line)
                except AttributeError:
                    # Skip unneeded headers
                    pass

    def parse_include(self, line):
        _, itp_file = line.strip().split()
        itp_file = itp_file.strip("\"")
        itp_file = Path(itp_file)
        if not itp_file.name.startswith("/"):
            itp_file = self.top_file.parent / itp_file
        self.top(itp_file)

    def parse_define(self, line):
        # Define statements define a variable name
        # and an optional value
        splitline = line.strip().split()
        name, value = "", ""
        if len(splitline) == 2:
            # A variable is defined without a value
            name = splitline[1]
            value = None
        elif len(splitline) == 3:
            # A variable is defined with a value
            _, name, value = splitline
            try:
                value = float(value)
            except ValueError:
                pass
        else:
            logging.debug(f"Couldn't parse define statement:\n"+
                          line.strip())
            logging.debug(f"{name=}, {value=}")

        # Check if the variables don't overwrite existing
        # namespace
        if name in globals() or name in locals():
            raise Exception("Use of global or local names is not allowed")

        self.variables[name] = value

    def parse_moleculetype(self, line):
        name, nrexcl = line.strip().split()
        molecule = Molecule(name=name)
        self.system.extend(molecule)
        self.active_molecule = molecule

    def parse_atoms(self, line):
        idx, _type, resi, resn, name, *_ = line.strip().split()
        idx = int(idx)
        resi = int(resi)
        atom = Atom(idx, _type, resi, resn, name)
        self.active_molecule.extend(atom)

    def parse_bonds(self, line):
        _from, to, bondtype, distance, fc = line.strip().split()

        _from = int(_from)
        to = int(to)
        distance = float(distance)
        # Sometimes these are variably defined
        fc = float(self.resolve(fc))

        if distance <= 0.47 and fc >= 1250:
            # Short distance, high force constant bonds are considered covalent
            self.active_molecule[_from].add_bond(to)
        else:
            # Long distance, low force constant bonds are considered elastic
            # regardless of bondtype
            self.active_molecule[_from].add_elastic(to)

    def parse_constraints(self, line):
        # Constraints are considered bonds
        _from, to, *_ = line.strip().split()
        _from = int(_from)
        to = int(to)
        self.active_molecule[_from].add_bond(to)

    def parse_molecules(self, line):
        name, n = line.strip().split()
        n = int(n)
        self.system.order(name, n)
