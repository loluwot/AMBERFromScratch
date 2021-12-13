import numpy as np
import math
from constants import *
from vpython import *
from networkx.algorithms import isomorphism
import pickle
from ortools.sat.python import cp_model
from collections import defaultdict
import string


RAD2DEG = 360/2/math.pi
def normalize(vec):
    return vec/np.linalg.norm(vec)

STANDARD_UNIT = 'pm'
scene = canvas()
class hashabledict(dict):
    def __hash__(self):
        return hash(tuple(sorted(self.items())))
class Atom:
    def __init__(self, name, pos_vec, units='ang', vel_vec=np.array([0,0,0])):
        conv_factor = UNITS[units]/UNITS[STANDARD_UNIT]
        name = name.upper()
        self.name = name
        self.pos_vec = pos_vec*conv_factor
        self.vel_vec = vel_vec*conv_factor
        self.basic_type = None
        self.type = None
        self.radius = ATOMS[name]['radius']*UNITS['pm']/UNITS[STANDARD_UNIT]
        self.mass = ATOMS[name]['amu']
        self.color = ATOMS[name]['color']
        self.ring_types = set()
        self.nhs = 0
        self.hav = 0

    def dist(self, other):
        return np.linalg.norm(self.pos_vec - other.pos_vec)

    def angle(self, other1, other2):
        d23 = np.linalg.norm(other2.pos_vec - self.pos_vec)
        d21 = np.linalg.norm(other1.pos_vec - self.pos_vec)
        d13 = np.linalg.norm(other2.pos_vec - other1.pos_vec)
        return np.arccos((d23**2 + d21**2 - d13**2)/(2*d23*d21))*RAD2DEG

    def __str__(self):
        return self.name
    def copy(self):
        return pickle.loads(pickle.dumps(self))

class Molecule:
    def __init__(self, atoms, adj_list=None):
        self.atoms = atoms
        self.atoms_to_id = {}
        for i, atom in enumerate(atoms):
            self.atoms_to_id[atom] = i  
        self.adj_list = [[] for __ in range(len(atoms))]
        if adj_list is not None:
            self.adj_list = adj_list
        else:
            for i, atom1 in enumerate(atoms):
                for j, atom2 in enumerate(atoms[:i]):
                    if atom1.dist(atom2) <= K*(atom1.radius + atom2.radius):
                        self.adj_list[i].append(j)
                        self.adj_list[j].append(i)
        self.nx_graph = nx.Graph()
        self.nx_graph.add_nodes_from([(i, {'name':self.atoms[i].name, 'n':len(self.adj_list[i])}) for i in range(len(self.atoms))])
        for i, neighbours in enumerate(adj_list):
            self.nx_graph.add_edges_from([(i, x) for x in neighbours])
        self.bos = {}
        self.resonance_states = []

    def bond_length(self, atom1, atom2):
        if isinstance(atom1, int):
            atom1 = self.atoms[atom1]
        if isinstance(atom2, int):
            atom2 = self.atoms[atom2]
        return atom1.dist(atom2)

    def bond_angle(self, atom1, atom2, atom3):
        #atom1 -- atom2 -- atom3
        if isinstance(atom1, int):
            atom1 = self.atoms[atom1]
        if isinstance(atom2, int):
            atom2 = self.atoms[atom2]
        if isinstance(atom3, int):
            atom3 = self.atoms[atom3]

        return atom2.angle(atom1, atom3)

    def torsion_angle(self, atom1, atom2, atom3, atom4):
        #atom1 -- atom2 -- atom3 -- atom4
        if isinstance(atom1, int):
            atom1 = self.atoms[atom1]
        if isinstance(atom2, int):
            atom2 = self.atoms[atom2]
        if isinstance(atom3, int):
            atom3 = self.atoms[atom3]
        if isinstance(atom4, int):
            atom4 = self.atoms[atom4]

        r21 = normalize(atom1.pos_vec - atom2.pos_vec)
        r23 = normalize(atom3.pos_vec - atom2.pos_vec)
        r34 = normalize(atom4.pos_vec - atom3.pos_vec)
        r2123 = normalize(np.cross(r21, r23))
        r2334 = normalize(-np.cross(r23, r34))
        sign = 1.0 if np.dot(r2123, r34) <= 0.0 else -1.0
        return np.arccos(np.dot(r2123,r2334))*sign*RAD2DEG

    def oop_angle(self, atom1, atom2, atom3, atom4):
        #atom2 -- atom4 -- atom1
        #           |
        #         atom3
        #oop of r41
        if isinstance(atom1, int):
            atom1 = self.atoms[atom1]
        if isinstance(atom2, int):
            atom2 = self.atoms[atom2]
        if isinstance(atom3, int):
            atom3 = self.atoms[atom3]
        if isinstance(atom4, int):
            atom4 = self.atoms[atom4]
        r41 = normalize(atom1.pos_vec - atom4.pos_vec)
        r42 = normalize(atom2.pos_vec - atom4.pos_vec)
        r43 = normalize(atom3.pos_vec - atom4.pos_vec)
        r4243 = normalize(np.cross(r42, r43))
        np.arcsin(np.dot(r4243, r41)/np.sin(self.bond_angle(atom2, atom4, atom3)))*RAD2DEG

    @classmethod
    def read_xyz(cls, filename, unit='ang'):
        #print(unit)
        conv_factor = UNITS[unit]/UNITS[STANDARD_UNIT]
        f = open(filename)
        atoms = []
        try:
            n_atoms = int(f.readline())
        except:
            raise Exception('Invalid format.')
        adj_list = [[] for __ in range(n_atoms)]
        f.readline() #skipping comments
        for i, l in enumerate(f):
            strings = []
            temp = ''
            last = ''
            for c in l:
                if c != ' ':
                    if last == ' ':
                        if len(temp) > 0:
                            strings.append(temp)
                        temp = c
                    else:
                        temp += c
                last = c
            strings.append(temp)
            a, x, y, z = strings
            pos_vec = None
            try:
                pos_vec = np.array([float(x)*conv_factor, float(y)*conv_factor, float(z)*conv_factor])
            except:
                raise Exception('Invalid format.')
            atoms.append(Atom(a, pos_vec, units=STANDARD_UNIT))
            for j, atom in enumerate(atoms[:-1]):
                if atom.dist(atoms[-1]) <= K*(atom.radius + atoms[-1].radius):
                    adj_list[i].append(j)
                    adj_list[j].append(i)
        return Molecule(atoms, adj_list)

    def __str__(self):
        temp = ''.join([str(atom) for atom in self.atoms])
        return temp + ' ' + str(self.adj_list)

    def draw(self, bond_radius=30):
        
        for index in range(len(self.resonance_states)):
            for atom in self.atoms:
                sphere(canvas=scene, pos=vector(*atom.pos_vec)+ vector(0,1,0)*index*700, radius=atom.radius/2, color=atom.color)
                #print(atom.ring_types)
            for i, bonds in enumerate(self.adj_list):
                for neighbour in bonds:
                    if neighbour < i:
                        a12 = self.atoms[neighbour].pos_vec - self.atoms[i].pos_vec
                        pos = self.atoms[i].pos_vec
                        #bo = self.bos[neighbour, i]
                        #print(self.resonance_states)
                        bo = self.resonance_states[index][neighbour, i]
                        #print(bo)
                        perp = hat(vector(*a12).cross(scene.forward))
                        increment = bond_radius*2/(bo+1)
                        for ii in range(bo):
                            cylinder(canvas=scene, pos=vector(*pos)-perp*(bond_radius-increment*(ii+1))+ vector(0,1,0)*index*700, axis=vector(*a12), radius=10)
                        #cylinder(canvas=scene, pos=vector(*pos), axis=vector(*a12), radius=10)
    def basic_typing(self):
        #print(self.nx_graph.nodes.data())
        for name, graph, id1 in BASIC_TYPES:
            GM = isomorphism.GraphMatcher(self.nx_graph, graph, node_match=isomorphism.categorical_node_match(['name', 'n'],[None, None]))
            for s in GM.subgraph_isomorphisms_iter():
                #print(name)
                self.atoms[list(s.keys())[id1]].basic_type = name
        for i, atom in enumerate(self.atoms):
            if atom.basic_type is None:
                atom.basic_type = '{}X{}'.format(atom.name, len(self.adj_list[i]))
                #print(atom.basic_type)
        #type n of hydrogens
        for i in range(len(self.atoms)):
            hs = 0
            for neighbour in self.adj_list[i]:
                if self.atoms[neighbour].name == 'H':
                    hs += 1
            self.atoms[i].nhs = hs

        
    def copy(self):
        return pickle.loads(pickle.dumps(self))

    def generate_vs(self, limit=1000):
        cur_states = []
        arr = [[] for _ in range(10000)]
        init_state = [0 for _ in range(len(self.atoms))]
        arr[0].append((init_state, 0))
        index = 0
        while len(cur_states) < limit:
            cur_states.extend(arr[index])
            for possible, a in arr[index]:
                for i, value in enumerate(possible):
                    #print(tup)
                    #print(value)
                    temp = possible[::]
                    if value+1 >= len(APS[self.atoms[i].basic_type]):
                        continue
                    temp[i] = value+1
                    #print(index - APS[self.atoms[i].basic_type][value][1] + APS[self.atoms[i].basic_type][value+1][1])
                    new_index = index - APS[self.atoms[i].basic_type][value][1] + APS[self.atoms[i].basic_type][value+1][1]
                    arr[new_index].append((temp, new_index))
                    #cur_states.append(temp)
            index += 1
            #print(len(cur_states))
        return cur_states
    def resolve_vs(self, vs):
        model = cp_model.CpModel()
        ac_vs = [APS[self.atoms[i].basic_type][x][0] for i, x in enumerate(vs)]
        #print(ac_vs)
        conn = [len(self.adj_list[i]) for i in range(len(self.atoms))]
        variables = {}
        for i in range(len(self.atoms)):
            for neighbour in self.adj_list[i]:
                if (min(neighbour, i), max(neighbour, i)) not in variables.keys():
                    variables[min(neighbour, i), max(neighbour, i)] = model.NewIntVar(1, 3, '{},{}'.format(min(neighbour, i), max(neighbour, i)))
        
        for i in range(len(self.atoms)):
            s = 0
            for neighbour in self.adj_list[i]:
                tup = (min(neighbour, i), max(neighbour, i))
                s += variables[tup]
            model.Add(s == ac_vs[i])
        
        class VarArraySolutionPrinter(cp_model.CpSolverSolutionCallback):
            """Print intermediate solutions."""
            def __init__(self, variables, limit):
                cp_model.CpSolverSolutionCallback.__init__(self)
                self.__variables = variables
                self.solutions = []
            def on_solution_callback(self):
                new_solution = {}
                for key, value in self.__variables.items():
                    new_solution[key] = self.Value(value)
                self.solutions.append(new_solution)
            def get_solutions(self):
                #print(self.solutions)
                return self.solutions
        
        
        print('A')
        solver = cp_model.CpSolver()
        callbacks = VarArraySolutionPrinter(variables, 2000)
        print('A')
        status = solver.SearchForAllSolutions(model, callbacks)
        print('A')
        print(solver.StatusName(status))
        # if status == cp_model.OPTIMAL:
        #     # for key, value in variables.items():
        #     #     #print('{} BO {}'.format(key, solver.Value(value)))
        #     #     variables[key] = solver.Value(value)
        #     # #print(self.adj_list)
        return callbacks.solutions
        # else:
        #     #print(':(')
        #     return -1
    def assign_bo(self):
        self.basic_typing()
        vss = self.generate_vs()
        min_index = 100000
        for vs, index in vss:
            #print(index)
            if index <= min_index:
                #print('----------')
                results = self.resolve_vs(vs)
                #print(results)
                for result in results:
                    if result != -1:
                        self.bos = result
                        min_index = index
                    self.resonance_states.append(result)
            else:
                break

    def ring_typing(self):
        ring_type = defaultdict(lambda : 4)
        rings = []
        atoms_in_ring = set()
        for ring in nx.cycle_basis(self.nx_graph):
            for node in ring:
                atoms_in_ring.add(node)

        for ring in nx.cycle_basis(self.nx_graph):
            ring_key = tuple(sorted(ring))
            rings.append(ring_key)
            for resonance in self.resonance_states:
                isPlanar = True
                externaldoubles = []
                bonds = []
                for i, node in enumerate(ring):
                    for neighbour in self.adj_list[node]:
                        if neighbour not in ring:
                            tup = (min(node, neighbour), max(node, neighbour))
                            if resonance[tup] == 2 and neighbour not in atoms_in_ring:
                                externaldoubles.append(node)
                                break
                    tup = (min(node, ring[(i+1) % len(ring)]), max(node, ring[(i+1) % len(ring)]))
                    bonds.append(resonance[tup])
                    if self.atoms[node  ].basic_type not in PLANAR_TYPES: 
                        isPlanar = False
                    self.atoms[node].ring_types.add('R{}'.format(len(ring)))
                external_bonds = bonds[::]
                for double in externaldoubles:
                    external_bonds.insert(double, 2)
                alternating = True
                hasExternalDouble = True
                hasConSingle = False
                allSingle = True
                #print(allSingle)
                last = bonds[-1]
                for bond in bonds:
                    if not bond == 1:
                        allSingle = False
                    if bond == last:
                        alternating = False
                    if bond == 1 and last == 1 and not hasConSingle:
                        hasConSingle = True
                    else:
                        hasConSingle = False
                    last = bond
                
                last = external_bonds[-1]
                for bond in external_bonds:
                    if bond == last:
                        hasExternalDouble = False
                if alternating:
                    ring_type[ring_key] = 1
                
                if ring_type[ring_key] == 4:
                    if hasExternalDouble:
                        ring_type[ring_key] = 3
                    elif hasConSingle:
                        ring_type[ring_key] = 2
                    elif allSingle:
                        ring_type[ring_key] = 5
                

                print(bonds)
                print(alternating)
                print(hasExternalDouble)
                
        for ring in rings:
            for node in ring:
                self.atoms[node].ring_types.add('AR{}'.format(ring_type[ring]))
            #return

    def find_ringprim(self, atoms, prim):
        filtered = []
        for atom in atoms:
            if prim == 'NR':
                if len(self.atoms[atom].ring_types) == 0:
                    filtered.append(atom)
            else:
                if prim in self.atoms[atom].ring_types:
                        filtered.append(atom)
        return filtered

    def find_ring(self, atoms, ring):
        unions = ring.split('.')
        total = set()
        for union in unions:
            
            intersects = union.split(',')
            temp = set(self.find_ringprim(atoms, intersects[0]))
            #print(temp)
            for intersect in intersects:
                temp = temp.intersection(set(self.find_ringprim(atoms, intersect)))
                #print(temp)
            total = total.union(temp)
        return list(total)
    def find_all(self, atoms, atom_str):
        vals = ['', '', '']
        mode = 0
        cur_str = ''
        for c in atom_str:
            if mode == 0:
                if c not in string.digits:
                    if c == '[':
                        vals[mode] = cur_str
                        cur_str = ''
                        mode = 2
                    else:
                        cur_str += c
                else:
                    vals[mode] = cur_str
                    cur_str = c
                    mode = 1
            elif mode == 1:
                if c in string.digits:
                    cur_str += c
                elif c == '[':
                    vals[mode] = int(cur_str)
                    cur_str = ''
                    mode = 2
                else:
                    raise ValueError
            elif mode == 2:
                if c != ']':
                    cur_str += c
                else:
                    vals[mode] = cur_str
                    break
        #print(vals)
        filtering = atoms[::]
        for i, val in enumerate(vals):
            temp = []
            if i == 0:
                possible_atoms = []
                if val == 'XX':
                    possible_atoms = ['C', 'N', 'O', 'S', 'P']
                elif val == 'XA':
                    possible_atoms = ['O', 'S']
                elif val == 'XB':
                    possible_atoms = ['N', 'P']
                elif val == 'XD':
                    possible_atoms = ['S', 'P']
                else:
                    possible_atoms = [val]
                for atom in filtering:
                    if self.atoms[atom].name in possible_atoms:
                        temp.append(atom)
            elif i == 1:
                if val == '':
                    continue
                else:
                    for atom in filtering:
                        if len(self.adj_list[atom]) == val:
                            temp.append(atom)
            elif i == 2:
                if val == '':
                    continue
                else:
                    temp = self.find_ring(filtering, val)
            #print(filtering)
            filtering = temp[::]
        return filtering


    def find_env(self, atoms, env):
        print(env)
        atom = ''
        index = 0
        for i, c in enumerate(env):
            if c != '(':
                atom += c
            else:
                index = i
                break
        if ')' in env:
            env = env[index+1:-1]
        else:
            return self.find_all(atoms, env)
        print(env)
        bracket = 0
        atom_strs = []
        last_str = ''
        for c in env:
            if c == '(':
                bracket += 1
            elif c == ')':
                bracket -= 1
            elif c == ',' and bracket == 0:
                atom_strs.append(last_str)
                last_str = ''
            else:
                last_str += c
        atom_strs.append(last_str)
        print(atom_strs)
        filtered_atoms = []
        possible_atoms = self.find_all(atoms, atom)
        for at in possible_atoms:
            possible = True
            for atstr in atom_strs:
                if len(self.find_env(self.adj_list[at], atstr)) < 1:
                    possible = False
                    break
            if possible:
                filtered_atoms.append(at)
        return filtered_atoms

    

