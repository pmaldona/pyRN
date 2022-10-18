import json
from typing import Set
import networkx as nx

def isSubset(set, subset) -> bool:
    for element in subset:
        if element not in set:
            return False
    return True

def isEqual(A: list, B: list) -> bool:
    if len(A) is not len(B):
        return False
    for a in A:
        if a not in B:
            return False
    return True


def getOrgsSmallerContaining(organizations, sm, speciesList):
    result = list()
    for i in range(len(organizations)):
        possibleSolutions = []
        root = sm[i]
        
        if i+1 < len(sm):
            for j in range(i+1,len(sm)):
                org = sm[j]
                if isSubset(org, root):
                    possibleSolutions.append(org)

            solutions = []
            for j in range(len(possibleSolutions)):
                _root = possibleSolutions[j]
                remove = False
                for k in range(len(possibleSolutions)):
                    _org = possibleSolutions[k]
                    if isEqual(_root, _org):
                        continue
                    if isSubset(_root, _org):
                        remove = True
                        break
                
                if remove == False:
                    solutions.append(str(_root))
                    continue
            
            if len(solutions) > 0:
                result.append((root, solutions))
    print("RESULT:")    
    print(result)
    return result


def build_hasse(sm, RN):

    ssm_ssms = RN.SsmStrSsmListSpArray
    species_list = RN.MpDf.index.values

    self_maintained = []
    for org in sm:
        self_maintained.append(org[1]['sp'])
    
    orgs = list()
    for org in sm:
        _org = list()
        for species in species_list:
            if species in org[1]['sp']:
                _org.append(1)
            else:
                _org.append(0)
        orgs.append(_org)

    G = nx.Graph()
    node_id_count = 0
    self_maintained.sort(key=len)
    length_dict = {}
    for semi_self in ssm_ssms:
        if length_dict.get(str(len(semi_self)), False):
            length_dict[str(len(semi_self))] = length_dict[str(len(semi_self))] + 1
        else :
            length_dict[str(len(semi_self))] = 1

    c = 0
    species_count = 0
    reaction_count = 0
    for _set in self_maintained:
        _len = len(_set)
        species_count = species_count + _len
        reactions = set()
        for spcs in _set:
            row = RN.MrDf.loc[spcs]
            for i in range(len(row)):
                if row[i] >= 1:
                    reactions.add(i)
            row = RN.MpDf.loc[spcs]
            for i in range(len(row)):
                if row[i] >= 1:
                    reactions.add(i)
        reaction_count = reaction_count + len(reactions)

        count = length_dict[str(_len)]
        _x = 0
        if count > 1:
            _x = -(75*(count-1)/2) + c*(75*(count-1))
            c += 1
        else:
            c = 0
        _y = -75*(_len-1)
        G.add_node(node_id_count, group = 1, label="O"+str(node_id_count), size=20, title=str(_set), x = _x, y = _y, fixed = json.loads('{ "x":false, "y":true}'))
        node_id_count += 1
        print(str(_set))
        print(node_id_count)
        print(G.nodes)
    
    avg_species = species_count/len(sm)
    avg_reactions = reaction_count/len(sm)

    edges = getOrgsSmallerContaining(orgs, self_maintained, species_list)

    for edgeList in edges:
        
        _from = -1
        for index in list(G.nodes):
            if str(edgeList[0]) == G.nodes[index]['title']:
                _from = index
                break
        for edge in edgeList[1:]:
            _to = -1
            for index in list(G.nodes()):
                print(edge[0])
                if str(edge[0]) == G.nodes[index]['title']:
                    _to = index
                    break
                
            G.add_edge(_from, _to, color="gray", smooth = False)


    nodes = list(G.nodes.items())
    print(nodes)
    nodes.reverse()
    _range = list(range(len(nodes)))
    _range.reverse()
    
    return G, avg_species, avg_reactions