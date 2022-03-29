from classes import *
import numpy as np
A = Atom('C', np.array([0,0,0]))
B = Atom('C', np.array([1,0,0]))
C = Atom('C', np.array([-1,1,0]))

print(f'{A.angle(B, C)} deg')