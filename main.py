import eel
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
        print(RN.prod)
        build_reaction_network(RN.mr, RN.mp)

@eel.expose
def calculate_orgs():
    if RN:
        RN.gen_basics()
        #RN.gen_syn_str()
        RN.gen_ssm_str()
        #print(RN.ssms)
        self_maintained = []
        
        for semi_self in RN.ssm_ssms:
            if RN.is_sm(semi_self):
                self_maintained.append(semi_self)
        
        hasse = Network('500px', '500px')
        hasse.toggle_physics(False)
        nx_graph = build_hasse(self_maintained, RN.ssm_ssms, RN.sp)
        hasse.from_nx(nx_graph)
        #print(r_net)
               
        nodes, edges, heading, height, width, options = hasse.get_network_data()
        return {"nodes": nodes, "edges": edges, "heading": heading, "height": height, "width": width, "options": options}
    return None

say_hello_py('Python World!')
eel.say_hello_js('Python World!')   # Call a Javascript function

print("[eel]: Start");

eel.start('index.html', mode=None)