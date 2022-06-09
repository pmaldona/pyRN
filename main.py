import json
import eel
import os
from CRNS import CRNS
import pandas as pd
from pyvis.network import Network
import tkinter as tk
from tkinter import filedialog
import networkx as nx
import platform

eel.init('static')

eel.say_hello_js('Python World!')   # Call a Javascript function

@eel.expose
def pick_file():
    global RN
    root = tk.Tk()
    root.withdraw()

    path = filedialog.askopenfilename()
    root.destroy()
    root = None
    if os.path.exists(path):
        #RN = RN.from_txt(path)
        RN = CRNS.form_txt(path)
        create_reaction_network(path)
        return path
    else:
        return 'Not valid folder'

def build_networkX_graph(path: str):
    # Generates 24 reaction columns even though there are only 12 reactions
    # My Interpretation: r1, p1, r2, p2, ... but why then save matrix twice?
    file = open(path)
    G = nx.MultiDiGraph()
    node_id_count = 0
    reactions = []
    for line in file.readlines():
        reaction = line.replace(" ","").replace("\t", "").replace("\n", "").replace(";", "").split(":")
        
        reaction_body = reaction.pop()
        reaction_node = reaction.pop()
        reactions.append(reaction_body)
        G.add_node(node_id_count, group = 2, label=reaction_node, size=20, title=reaction_node)
        node_id_count += 1

    #print(RN.mr)
    #print(RN.mp)

    df_mr = pd.DataFrame(data=RN.mr)
    reaction_size = df_mr.columns.size
    #print("reaction count: %s"%reaction_size)
    #for i in range(reaction_size):
    #    r_net.add_node(i, title='R%s'%i, color='red')

    #print(RN.reac)
    #print('RN.prod: len = %s'%len(RN.prod))
    #print('RN.prod[0]: len = %s'%len(RN.prod[0]))
    #print(RN.prod)

    species = df_mr.index.values
    #print("species count: %s"%len(species))
    species_id_map = {}
    for i in range(len(species)):
        G.add_node(node_id_count, group = 1, label=species[i], size=20, title=species[i])
        species_id_map[species[i]] = node_id_count
        node_id_count += 1
    #    stoich = list(df_mr.loc[species[i]])
    #    for j in range(reaction_size):
    #        if stoich[j] != 0.0:
    #            if j%2 == 0:
    #                r_net.add_edge(reaction_size+i, j/2)
    #            else:
    #                r_net.add_edge(reaction_size+i, (j-1)/2)

    reaction_id = 0
    for reaction in reactions:
        reactants, products = reaction.split('>')
        reversible = reactants[-1:] == "="
        reactants = reactants[0:-1]
        reactants = reactants.split("+")
        products = products.split("+")
        for reactant in reactants:
            num = reactant[0] if reactant[0].isdigit() else 1
            reac = reactant[1:] if reactant[0].isdigit() else reactant
            G.add_edge(species_id_map[reac], reaction_id, label=f'{num}', color="gray", font= {
				"size": 10,
				"align": "top"
				}, smooth = False, title=num)
            if reversible:
                G.add_edge(reaction_id, species_id_map[reac],label=f'{num}', color="gray", font= {
				"size": 10,
				"align": "top"
				}, smooth = False, title=num)
        
        for product in products:
            num = product[0] if product[0].isdigit() else 1
            prod = product[1:] if product[0].isdigit() else product
            G.add_edge(species_id_map[prod], reaction_id, label=f'{num}', color="gray", font= {
				"size": 10,
				"align": "top"
				}, smooth = False, title=num)
            if reversible:
                G.add_edge(reaction_id, species_id_map[prod],label=f'{num}', color="gray", font= {
				"size": 10,
				"align": "top"
				}, smooth = False, title=num)
        #string[0].isdigit()
        #print(reactants)
        #print(reversible)
        #print(products)

        reaction_id += 1

    return G
         

def create_reaction_network(path: str):
    # reaction network object
    r_net = Network('500px', '500px', directed =True)
    r_net.toggle_physics(False)
    nx_graph = build_networkX_graph(path)
    r_net.from_nx(nx_graph)
    #print(r_net)
               
    r_net.show('./static/templates/reaction_networks/graph.html')

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

            result.append((root, solutions))
        

    return result


def build_hasse(sm, species_list):

    orgs = list()
    for org in sm:
        _org = list()
        for species in species_list:
            if species in org:
                _org.append(1)
            else:
                _org.append(0)
        orgs.append(_org)

    G = nx.Graph()
    node_id_count = 0
    sm.sort(key=len)
    length_dict = {}
    for semi_self in RN.ssms:
        if length_dict.get(str(len(semi_self)), False):
            length_dict[str(len(semi_self))] = length_dict[str(len(semi_self))] + 1
        else :
            length_dict[str(len(semi_self))] = 1

    for set in sm:
        _len = len(set)
        count = length_dict[str(_len)]
        _x = 75*(count-1)-(75*(count-1)/2)
        _y = -75*(_len-1)
        G.add_node(node_id_count, group = 1, label="O"+str(node_id_count), size=20, title=str(set), x = _x, y = _y, fixed = json.loads('{ "x":false, "y":true}'))
        node_id_count += 1
    
    edges = getOrgsSmallerContaining(orgs, sm, species_list)
    
    for edgeList in edges:
        _from = -1
        for index in list(G.nodes()):
            if str(edgeList[0]) == G.nodes[index]['title']:
                _from = index
                break
        for edge in edgeList[1:]:
            _to = -1
            for index in list(G.nodes()):
                if str(edge[0]) == G.nodes[index]['title']:
                    _to = index
                    break
                
            G.add_edge(_from, _to, color="gray", smooth = True)


    nodes = list(G.nodes.items())
    nodes.reverse()
    _range = list(range(len(nodes)))
    _range.reverse()
    
    return G

@eel.expose
def calculate_orgs():
    if RN:
        RN.gen_basics()
        #RN.gen_syn_str()
        RN.gen_ssm_str()
        #print(RN.ssms)
        self_maintained = []
        
        for semi_self in RN.ssms:
            if RN.is_sm(semi_self):
                self_maintained.append(semi_self)
        
        hasse = Network('500px', '500px')
        hasse.toggle_physics(False)
        nx_graph = build_hasse(self_maintained, RN.sp)
        hasse.from_nx(nx_graph)
        #print(r_net)
               
        hasse.show('./static/templates/hasse/graph.html')
    return True

print(platform.system())

if platform.system() == "Darwin":
    eel.browsers.set_path('electron', 'node_modules/electron/dist/Electron.app/Contents/MacOS/Electron')
elif platform.system() == "Windows":
    eel.browsers.set_path('electron', 'node_modules/electron/dist/Electron.exe')
else:
    eel.browsers.set_path('electron', 'node_modules/electron/dist/Electron')
eel.start('templates/index.html', mode='electron', size=(800, 700), jinja_templates='templates', args=['--disable-http-cache'])