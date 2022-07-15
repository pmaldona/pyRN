import networkx as nx

def build_reaction_network(mr, mp):
    G = nx.MultiDiGraph()

    print(mr)
    print(mp)