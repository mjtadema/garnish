from skewer.utils import clean_path
import sys
import os

# check if module was called by pymol, if yes run skewer.main
main_initfile = sys.modules['__main__'].__file__

# get the name of the parent module
main_modulename = os.path.basename(clean_path(os.path.join(main_initfile, os.pardir)))

if main_modulename == 'pymol':
    from skewer.skewer import load
    load()
