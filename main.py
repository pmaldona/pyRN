import eel
import matplotlib.pyplot as plt
import matplotlib as mpl
import base64
from io import BytesIO
import numpy as np
from pyvis.network import Network
from pyRN import pyRN
from network.network import build_reaction_network

from organisations.hasse import build_hasse

print("[eel]: Init");

eel.init('static')

@eel.expose                         # Expose this function to Javascript
def say_hello_py(x):
    print('Hello from %s' % x)

@eel.expose
def openFile(path):
    global RN
    print(path)
    extension = str(path).split('.')[-1]
    if extension == 'txt':
        RN = pyRN.from_txt(path)
        return True
    elif extension == 'sbml' or extension == 'xml':
        RN = pyRN.from_sbml(path)
        return True
    else:
        return False

@eel.expose
def gen_network():
    if RN:
        network = RN.display_RN()
        nodes, edges, heading, height, width, options = network.get_network_data()
        species_count = len(RN.sp)
        reaction_count = len(RN.mr.columns)
        return {"nodes": nodes, "edges": edges, "heading": heading, "height": height, "width": width, "options": options, "species_count": species_count, "reaction_count": reaction_count}
        #return {"network": str(network)}

@eel.expose
def gen_synergetic():
    if RN:
        RN.gen_syn_str()
        network = RN.display_str(RN.syn_str)
        nodes, edges, heading, height, width, options = network.get_network_data()
        basic = 0
        closed = 0
        basic_closed = 0
        ssm = 0
        orgs = 0
        union = 0
        syn_union = 0
        conn_orgs = 0
        for node in nodes:
            if node['shape'] == 'dot':
                basic = basic + 1
            elif node['shape'] == 'square':
                closed = closed + 1
            if node['color'] == 'red':
                basic_closed = basic_closed + 1
            elif node['color'] == 'blue':
                ssm = ssm + 1
            elif node['color'] == 'green':
                orgs = orgs + 1
        for edge in edges:
            if edge['color'] == 'green':
                union = union + 1
            elif edge['color'] == 'blue':
                syn_union = syn_union + 1
        return {"nodes": nodes, "edges": edges, "heading": heading, "height": height, "width": width, "options": options, "basic": basic, "closed": closed, "basic_closed":basic_closed, "ssm":ssm, "orgs": orgs, "union": union, "syn_union":syn_union, "conn_orgs":conn_orgs}
        #return {"network": str(network)}

@eel.expose
def gen_protosynergetic():
    if RN:
        RN.gen_atoms()
        RN.gen_mgen()
        RN.all_syn()
        network = RN.display_syn()
        nodes, edges, heading, height, width, options = network.get_network_data()
        transitions = len(RN.syn)
        generators = len(nodes) - transitions
        return {"nodes": nodes, "edges": edges, "heading": heading, "height": height, "width": width, "options": options, "generators": generators, "transitions": transitions}
        #return {"network": str(network)}

@eel.expose
def calculate_orgs():
    if RN:
        RN.gen_atoms()
        #RN.gen_syn_str()
        RN.gen_ssm_str()
        #print(RN.ssms)
        self_maintained = []
        
        for semi_self in RN.ssm_ssms:
            if RN.is_sm(semi_self):
                self_maintained.append(semi_self)
        
        hasse = Network('500px', '500px')
        hasse.toggle_physics(False)
        nx_graph, species_count, reaction_count = build_hasse(self_maintained, RN)
        hasse.from_nx(nx_graph)
        #print(r_net)

        nodes, edges, heading, height, width, options = hasse.get_network_data()

        return {"nodes": nodes, "edges": edges, "heading": heading, "height": height, "width": width, "options": options, "species_count": species_count, "reaction_count": reaction_count}
    return None

@eel.expose
def plot_basics_sp():
    # Data for plotting
    t = np.arange(0.0, 2.0, 0.01)
    s = 1 + np.sin(2 * np.pi * t)

    fig, ax = plt.subplots()
    ax.plot(t, s)

    ax.set(xlabel='time (s)', ylabel='voltage (mV)',
        title='About as simple as it gets, folks')
    ax.grid()

    tmpfile = BytesIO()
    fig.savefig(tmpfile, format='png')
    encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

    return encoded


@eel.expose
def plot_basics_r():
    RN.plot_basic_r_presence()

@eel.expose
def plot_stoichiometry():
    sp=RN.sp
    sp_i=range(len(sp))
    r=RN.mp.columns
    r_i=range(len(r))
    S=RN.mp.iloc[sp_i,r_i]-RN.mr.iloc[sp_i,r_i]
    fig, ax = plt.subplots(figsize = (10, 5))
    
    #  Ploting
    
    ax.matshow(S, cmap=mpl.colormaps['viridis'])
    
    
    for i in range(len(sp_i)):
        for j in range(len(r_i)):
            c = S.iloc[i,j]               
            ax.text(j, i, str(c), va='center', ha='center')
    
    sp_ticks=np.array(list(map(lambda x: x+ 0.5,range(len(sp_i)))))
    r_ticks=np.array(list(map(lambda x: x+ 0.5,range(len(r_i)))))
    
    ax.set_yticks(sp_ticks,sp)
    ax.set_xticks(r_ticks,r)

    tmpfile = BytesIO()
    fig.savefig(tmpfile, format='png')
    encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

    return encoded
    #plt.show()

@eel.expose
def plot_concentrations(ti=0, tf=50, steps=100, cutoff=0.1, i_sp=None, rt=None):
    if RN:
        RN.set_model_ma(i_sp,rt)
        RN.run_model(ti,tf,steps,cutoff)
        print(RN.con)
        print(RN.rate)
        print(RN.abst)
        print(RN.a_sp)
        print(RN.a_r)
        fig, ax = plt.subplots(1, 1, figsize = (10, 5))
        ax.plot(RN.con) #cmap=mpl.colormaps['viridis'])
        ax.title.set_text("Concentrations")
        tmpfile = BytesIO()
        fig.savefig(tmpfile, format='png')
        encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

        return encoded

@eel.expose
def plot_rates(ti=0, tf=50, steps=100, cutoff=0.1, i_sp=None, rt=None):
    if RN:
        RN.set_model_ma(i_sp,rt)
        RN.run_model(ti,tf,steps,cutoff)
        print(RN.con)
        print(RN.rate)
        print(RN.abst)
        print(RN.a_sp)
        print(RN.a_r)
        fig, ax = plt.subplots(1, 1, figsize = (10, 5))
        ax.plot(RN.rate) #cmap=mpl.colormaps['viridis'])
        ax.title.set_text("Rates")
        tmpfile = BytesIO()
        fig.savefig(tmpfile, format='png')
        encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

        return encoded

@eel.expose
def random_network(): 
    global RN
    RN=pyRN.rand_gen_whit_inflow()
    return True

#say_hello_py('Python World!')
#eel.say_hello_js('Python World!')   # Call a Javascript function

print("[eel]: Start");

eel.start('index.html', mode=None)