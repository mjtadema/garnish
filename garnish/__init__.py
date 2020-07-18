# Copyright 2019-2019 the garnish authors. See copying.md for legal info.
from pathlib import Path
import sys

from .parser import Parser

# Check if importing from pymol
try:
    # check if module was called by pymol
    main_module_name = Path(sys.modules['__main__'].__file__).parent.name
except AttributeError:
    # importing from an interpreter like ipython raises this error
    main_module_name = None

print(main_module_name)
if main_module_name == 'pymol':
    from . import extensions
