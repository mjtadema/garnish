# Structure of system dictionary

This is the structure of a system dictionary (e.g.: `parse_tpr` and `parse_top`):

``` 
system = {
    blocks: {
        block_id: {
            'moltype': str,
            'n_molecules': int,
        },
        ...,
        ...,
    },
    topology: {,
        moltype: {
            'connectivity': {
                'bonds': [(atom1, atom2), (...), ...],
                'constr': [(atom1, atom2), (...), ...],
                'harmonic': [(atom1, atom2), (...), ...],
                'vsiten': [(atom1, atom2), (...), ...],
            },
            'n_atoms': int,
            'atomtypes': {atom1: type1, atom2: type2, ...},
            'backbone': [atom1, atom2, atom3, ...],
        },
        ...,
        ...,
    },
},
```

# Structure of system graph

This is the structure of `System.graph`:

``` 
System.graph:
    - nodes:
        (atom_id, {'moltype': A, 'block': B, 'atomtype': C})
        ...
        ...
    - edges:
        (atom1_id, atom2_id, {'type': 'bonds'|'constr'|'harmonic'}
        ...
        ...
```
