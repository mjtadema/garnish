
# coding: utf-8

# In[38]:


import subprocess, shlex


# In[39]:


from pathlib import Path


# In[40]:


import re,io


# In[ ]:


import networkx as nx


# In[41]:


tpr = Path(sys.argv[1])
assert tpr.is_file()


# In[63]:



gmxdump = "/usr/bin/gmx dump -s "+str(tpr.absolute())
gmxdump = subprocess.Popen(shlex.split(gmxdump), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

p_grep = re.compile(".*\#atoms|.*\#beads.*|.*moltype.*|.*\#molecules.*|.*\(BONDS\).*|.*\(CONSTR\).*|.*\(HARMONIC\).*")

regexp_all = {
    'molid': re.compile("^\s+moltype\s+=\s+(\d+)"),
    'occs': re.compile("^\s+\#molecules\s+=\s+(\d+)"),
    'bead_nr': re.compile("^\s+\#atoms_mol\s+=\s+(\d+)"),
    'total_beads': re.compile("^\s+\#atoms+\s+=+\s+(\d+)"),
}
regexp_bonds = {
    'bonds': re.compile("^\s+\d+\s\w+=\d+\s\(BONDS\)\s+(\d+)\s+(\d+)"),
    'constr': re.compile("^\s+\d+\s\w+=\d+\s\\(CONSTR\)\s+(\d+)\s+(\d+)"),
    'harmonic': re.compile("^\s+\d+\s\w+=\d+\s\\(HARMONIC\)\s+(\d+)\s+(\d+)")
}
regexp_is_mol = re.compile("^\s+moltype\s+\((\d+)\):")

regex_data = {
    k: []
    for k in regexp_all.keys()
}

molecules = {}

reading_header = True
for line in io.TextIOWrapper(gmxdump.stdout, encoding="utf-8"):
    if p_grep.match(line):
        if reading_header:
            for k, p in regexp_all.items():
                if p.match(line):
                    regex_data[k] = p.findall(line)[0]
            if regexp_is_mol.match(line):
                reading_header = False
                molid = regexp_is_mol.findall(line)[0]
                bonds = {
                    k:[]
                    for k in regexp_bonds
                }
        else:
            if not regexp_is_mol.match(line):
                for k, p in regexp_bonds.items():
                    if p.match(line):
                        bond = p.findall(line)[0]
                        bonds[k].append(bond)
            else:
                molecules[molid] = bonds
                molid = regexp_is_mol.findall(line)[0]
                bonds = {
                    k:[]
                    for k in regexp_bonds
                }
                
for molid, molecule in molecules.items():
    for bondtype, bonds in molecule.items():
        g = nx.Graph()
        g.add_edges_from(bonds)
        molecules[molid][bondtype] = g
                

