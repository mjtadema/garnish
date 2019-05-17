#!/usr/bin/env python
# coding: utf-8

# In[121]:


import re, subprocess, shlex, io


# In[96]:


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


# In[95]:


regex_data = {
    k: []
    for k in regexp_all.keys()
}


# In[ ]:


molecules = {}


# In[101]:


with open("./test/gmxdump", 'r') as f:
    for line in f:
        for k, p in regexp_all.items():
            if p.match(line):
                regex_data[k] = p.findall(line)[0]
        if regexp_is_mol.match(line):
            break
    molid = regexp_is_mol.findall(line)[0]
    bonds = {
        k:[]
        for k in regexp_bonds
    }
    for line in f:
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

