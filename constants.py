from vpython import *
import networkx as nx
ATOMS = {}
ATOMS.update(dict.fromkeys(['C', 6], {'radius':70, 'amu':12.01, 'color':vector(0.3, 0.3, 0.3)}))
ATOMS.update(dict.fromkeys(['H', 1], {'radius':25, 'amu':1.008, 'color':color.white}))
ATOMS.update(dict.fromkeys(['O', 8], {'radius':60, 'amu':16.0, 'color':color.red}))
ATOMS.update(dict.fromkeys(['N', 7], {'radius':65, 'amu':14.01, 'color':color.blue}))

K = 1.2
class Node:
    def __init__(self, neighbours):
        self.neighbours = neighbours
        
UNITS = {'pm':1, 'A':100, 'ang':100, 'nm':1000}

CNR_graph = nx.Graph()
COO_graph = nx.Graph()
NNR_graph = nx.Graph()
NITRO_graph = nx.Graph()
P1O_graph = nx.Graph()

CNR_graph.add_nodes_from([(0, {'name':'C', 'n':1}), (1, {'name':'N', 'n':2})])
CNR_graph.add_edge(0, 1)
COO_graph.add_nodes_from([(0, {'name':'C', 'n':3}), (1, {'name':'O', 'n':1}), (2, {'name':'O', 'n':2})])
COO_graph.add_edges_from([(0, 1), (0, 2)])
NNR_graph.add_nodes_from([(0, {'name':'N', 'n':1}), (1, {'name':'N', 'n':2})])
NNR_graph.add_edge(0, 1)
NITRO_graph.add_nodes_from([(0, {'name':'N', 'n':3}), (1, {'name':'O','n':1}), (2, {'name':'O','n':1})])
NITRO_graph.add_edges_from([(0,1), (0, 2)])
P1O_graph.add_nodes_from([(0, {'name':'C', 'n':2}), (1, {'name':'N', 'n':3}), (2, {'name':'O', 'n':1}), (3, {'name':'C', 'n':2})])

BASIC_TYPES = [('CNR', CNR_graph, 0), ('COO', COO_graph, 0), ('NNR1', NNR_graph, 0), ('NNR2', NNR_graph, 1), ('NITRO', NITRO_graph, 0), ('P1ON', P1O_graph, 1), ('P1OO', P1O_graph, 2)]


APS = {}
APS.update(dict.fromkeys(['HX1', 'FX1', 'CLX1', 'BRX1', 'IX1'], [(1, 0), (0, 64), (1, 64)]))
APS['CNR'] = [(3, 0), (4, 1), (5, 32)]
APS['CX1'] = [(4, 0), (3, 1), (5, 32)]
APS['COO'] = [(4, 0), (5, 32), (6, 32)]
APS.update(dict.fromkeys(['CX2', 'CX3', 'CX4'], [(4, 0), (3, 32), (5, 32), (2, 64), (6, 64)]))
APS['NNR1'] = [(2, 0), (3, 0)]
APS['NNR2'] = [(4, 0), (3, 1)]
APS['NX1'] = [(3,0), (2,3), (4,32)]
APS['NX2'] = [(3, 0), (4, 2), (2, 4)]
APS['NITRO'] = [(5, 0), (4, 32), (6, 32), (3, 64)]
APS['P1ON'] = [(4, 0), (3, 1)]
APS['NX3'] = [(3, 0), (4, 1), (5, 2), (2, 32)]
APS['NX4'] = [(3, 0), (4, 64), (2, 64)]
APS['P1OO'] = [(1, 0), (0, 1)]
APS['OX1'] = [(2, 0), (1, 1), (3, 64)]
APS['OX2'] = [(2 , 0), (1, 32), (3, 64)]

default_fingerprint = {'fingerprint_length': 1024}
FINGERPRINT_HASH_CONSTANT = 10000
PLANAR_TYPES = ['CX3', 'NX2', 'NX3', 'OX2']
