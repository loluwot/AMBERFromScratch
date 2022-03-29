from classes import *
# from geom import *
import numpy as np
from collections import Counter, defaultdict
mol = Molecule.read_xyz('benzoic_acid.xyz')
mol.assign_bo()
mol.ring_typing()
print(mol.resonance_states)
print(mol.fingerprint(radius=4))
D = defaultdict(lambda: [])
for i, v in enumerate(mol.node_messages):
    D[v].append(mol.atoms[i].basic_type)
print(D)


# for atom in mol.atoms:
#     print(atom.basic_type)
# mol.draw()
# print(mol.find_all([i for i in range(len(mol.atoms))], 'C[R5,R6,AR1.AR2.AR3]'))
# mol.find_env([])