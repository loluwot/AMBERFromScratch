from classes import *
from geom import *
import numpy as np
mol = Molecule.read_xyz('benzoic_acid.xyz')
mol.assign_bo()
mol.ring_typing()
for atom in mol.atoms:
    print(atom.ring_types)
#mol.draw()

print(mol.find_all([i for i in range(len(mol.atoms))], 'C[R5,R6,AR1.AR2.AR3]'))
mol.find_env([])