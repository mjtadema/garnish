# Structure of system dictionary

This is the output of `parse_tpr` and `parse_top`:

``` 
system = {
    blocks: {
        block_id: {
            'moltype': str
            'n_molecules': int
        }
        ...
        ...
    }
    topology: {
        moltype: {
            'connectivity': {
                'bonds': [(atom1, atom2), (...), ...]
                'constr': [(atom1, atom2), (...), ...]
                'harmonic': [(atom1, atom2), (...), ...]
            }
            'n_atoms': int
            'backbone': [atom1, atom2, atom3, ...]
        }
        ...
        ...
    }
}
```

# Structure of graphs dictionary

This is the output of `make_graphs`:

``` 
bond_graphs = {
    'molecule_id': {
        'bonds': <graph>
        'constr': <graph>
        'harmonic': <graph>
    }
    ...
    ...
}
```
