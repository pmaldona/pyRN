from enum import Enum
from typing import Union
from pandas import DataFrame
from pyRN import pyRN
from pyvis.network import Network
import networkx as nx
from bitarray import bitarray as bt
from .WalkTypes import WalkTypes

from ..misc.CaptureOut import Capturing

class State:
    def __init__(self) -> None:
        self.set_new_rn(None, None)
    
    def set_new_rn(self, path: Union[str, None], RN: Union[pyRN, None]):
        self.file_path: Union[str, None] = path
        self.reaction_network: Union[pyRN, None] = RN
        self.constants_file: Union[str, None] = None
        self.rates_file: Union[str, None] = None
        self.statistics = None
        self.organizations: Union[list, None] = None
        self.network_graph: Union[Network, None] = None
        self.protosynergies_graph: Union[Network, None] = None
        self.hasse_graph: Union[Network, None] = None
        self.synergies_graph: Union[Network, None] = None
        self.model: bool = False
        self.simple_rw: bool = False
        self.random_walk_type: WalkTypes = WalkTypes.SIMPLE_RANDOM_WALK
        self.simulation_parameters: list = []
    
    def save_parameters(self, timeStart, timeFinal, steps, cutoff, w, l, d, nmin, n, trys, save, fname, convPert, keys) -> bool:
        self.simulation_parameters = [timeStart, timeFinal, steps, cutoff, w, l, d, nmin, n, trys, save, fname, convPert, keys]
        return True

    
    def get_parameters(self) -> list:
        return self.simulation_parameters


    def save_network(self, path: str) -> bool:
        if self.reaction_network != None:
            try:
                self.reaction_network.saveToText(path)
                if self.model == True:
                    self.reaction_network.saveInitCondToText()
            except:
                return False
            else:
                return True
        else:
            return False
        
    def open_file(self, file_path: str) -> bool:
        RN = None
        extension = file_path.split('.')[-1]
        if extension == 'txt':
            try:
                RN = pyRN.setFromText(file_path)
            except:
                pass
            
        elif extension == 'sbml' or extension == 'xml':
            try:
                RN = pyRN.setFromSbml(file_path)
            except:
                pass
        else:
            return False
        if RN == None:
            return False
        else:
            self.set_new_rn(file_path, RN)
            return True
    
    def generate_random_network(self, random_species=12, random_vector=None) -> bool:
        RN = None
        RN = pyRN.setSimpleRandomgenerate(Ns=random_species,rv=random_vector)
        if RN == None:
            return False
        else:
            self.set_new_rn("RandomNetwork", RN)
            return True
        
    def add_species(self, Nse, p, extra, m, l) -> bool:
        if Nse == 0:
            Nse = None
        if extra == 0:
            extra = None
        if self.reaction_network != None:
            try:
                self.reaction_network.setExtraRandomgenerated(p=p,m=m,Nse=Nse,extra=extra,l=l)
            except:
                return False
            else:
                return True
        else:
            return False
    
    def add_inflow(self, extra) -> bool:
        if self.reaction_network != None:
            try:
                self.reaction_network.setExtraRandomgeneratedInflow(extra=extra)
            except:
                return False
            else:
                return True
        else:
            return False

    def add_outflow(self, extra) -> bool:
        if self.reaction_network != None:
            try:
                self.reaction_network.setExtraRandomgeneratedOutflow(extra=extra)
            except:
                return False
            else:
                return True
        else:
            return False
    
    def generate_organizations(self) -> bool:
        if self.reaction_network != None:
            try:
                self.reaction_network.setGenerators()
                self.reaction_network.setSynStr()
                self.organizations = self.reaction_network.SynStrOrgListBtArray
                net=self.reaction_network.getHasseNxFromBtList(self.reaction_network.SynStrOrgListBtArray,setlabel="L")
                rn: pyRN = self.reaction_network
                net = nx.relabel_nodes(net, lambda x: str(rn.getIndArrayFromBt(bt(x))))
                nt = Network('500px', '500px',directed=False,notebook=False)
                nt.toggle_physics(False)
                nt.from_nx(net)
                self.hasse_graph = nt
            except:
                return False
            else:
                return True
        else:
            return False

    def generate_network_complete(self) -> bool:
        if self.reaction_network != None:
            try:
                network = self.reaction_network.getRnDisplayPv()
                self.network_graph = network
            except:
                return False
            else:
                return True
        else:
            return False

    def generate_synergetic_structures_complete(self) -> bool:
        if self.reaction_network != None:
            try:
                if self.organizations == None:
                    self.generate_organizations()
                self.synergies_graph = self.reaction_network.getStrDisplayPv(self.reaction_network.SynStrNx)
            except:
                return False
            else:
                return True
        else:
            return False
    
    def generate_protosynergetic_structures_complete(self) -> bool:
        if self.reaction_network != None:
            try:
                if self.organizations == None:
                    self.generate_organizations()
                self.reaction_network.setMgen()
                self.reaction_network.setSyn()
                self.protosynergies_graph = self.reaction_network.displaySynPv()
            except:
                return False
            else:
                return True
        else:
            return False

    def init_model(self, SpConFileNameStr=None, KConstFileNameStr=None, ti=0, tf=50, steps=100, cutoff=0.1) -> bool:
        if self.reaction_network != None:
            try:
                self.reaction_network.setMakModelFromFile(SpConFileNameStr ,KConstFileNameStr ,cutoff)
                self.reaction_network.runMakModel(ti,tf,steps,cutoff)
                self.model = True
                self.constants_file = SpConFileNameStr
                self.rates_file = KConstFileNameStr
            except:
                return False
            else:
                return True
        else:
            return False
    
    def set_random_walk_type(self, random_walk_type: WalkTypes) -> bool:
        self.random_walk_type = random_walk_type
        return True

    def init_simple_random_walk(self, w=range(10),l=10,d=1,nmin=3, fname="rand_walk.json") -> bool:
        if self.reaction_network != None:
            if self.organizations == None:
                if self.generate_organizations() == False:
                    return False
            try:
                self.reaction_network.setRwSimple(None,w,l,d,nmin,fname=fname)
                self.simple_rw = True
            except:
                return False
            else:
                return True
        else:
            return False
    
    def init_mak_random_walk(self, w=range(10),l=10,n=500, trys=10, save=True, fname="rand_mak_walk.json") -> bool:
        if self.reaction_network != None:
            if self.organizations == None:
                if self.generate_organizations() == False:
                    return False
            
            try:
                self.reaction_network.setMakRw(None, w=w,l=l,n=n,trys=trys,sim_save=save,fname=fname)
                self.simple_rw = True
            except:
                return False
            else:
                return True
        else:
            return False
        
    def get_current_state(self) -> dict:
        return {
            "file_path": self.file_path,
            "rates_path": self.rates_file,
            "concentrations_path": self.constants_file,
            "statistics": self.statistics,
            "RN": self.reaction_network,
            "network": self.network_graph,
            "protosynergies": self.protosynergies_graph,
            "hasse": self.hasse_graph,
            "synergies": self.synergies_graph
        }
    
    def get_walk_type(self) -> WalkTypes:
        return self.random_walk_type

    def get_reactions(self) -> list:
        if self.reaction_network:
            with Capturing() as output:
                self.reaction_network.printRp()
                return output
        else:
            return list()

    def get_network(self) -> Union[Network, None]:
        if self.network_graph != None:
            return self.network_graph
        elif self.network_graph == None and self.reaction_network != None:
            if self.generate_network_complete():
                return self.network_graph
            else:
                return None
        else:
            return None
    
    def get_synergetic_structure(self) -> Union[Network, None]:
        if self.synergies_graph != None:
            return self.synergies_graph
        elif self.synergies_graph == None and self.reaction_network != None:
            if self.generate_synergetic_structures_complete():
                return self.synergies_graph
            else:
                return None
        else:
            return None
    
    def get_protosynergetic_structure(self) -> Union[Network, None]:
        if self.protosynergies_graph != None:
            return self.protosynergies_graph
        elif self.protosynergies_graph == None and self.reaction_network != None:
            if self.generate_protosynergetic_structures_complete():
                return self.protosynergies_graph
            else:
                return None
        else:
            return None
    
    def get_hasse_structure(self) -> Union[Network, None]:
        if self.hasse_graph != None:
            return self.hasse_graph
        elif self.hasse_graph == None and self.reaction_network != None:
            if self.generate_organizations():
                return self.hasse_graph
            else:
                return None
        else:
            return None
    
    def get_concentrations_run(self, SpConFileNameStr=None, KConstFileNameStr=None, ti=0, tf=50, steps=100, cutoff=0.1) -> Union[DataFrame, None]:
        if self.model == True and self.reaction_network != None:
            if self.constants_file != SpConFileNameStr or self.rates_file == KConstFileNameStr:
                self.init_model(SpConFileNameStr, KConstFileNameStr, ti, tf, steps, cutoff)
            return self.reaction_network.SpConDf
        elif self.model == False and self.reaction_network != None:
            self.init_model(SpConFileNameStr, KConstFileNameStr, ti, tf, steps, cutoff)
            return self.reaction_network.SpConDf
        else:
            return None
    
    def get_rates_run(self, SpConFileNameStr=None, KConstFileNameStr=None, ti=0, tf=50, steps=100, cutoff=0.1) -> Union[DataFrame, None]:
        if self.model == True and self.reaction_network != None:
            if self.constants_file != SpConFileNameStr or self.rates_file == KConstFileNameStr:
                self.init_model(SpConFileNameStr, KConstFileNameStr, ti, tf, steps, cutoff)
            return self.reaction_network.RpRateDf
        elif self.model == False and self.reaction_network != None:
            self.init_model(SpConFileNameStr, KConstFileNameStr, ti, tf, steps, cutoff)
            return self.reaction_network.RpRateDf
        else:
            return None
        
    def get_simple_rw(self, walk_range=range(10),l=10,d=1,nmin=3, fname="rand_walk.json"):
        if self.reaction_network == None:
            return None
        else:
            try:
                self.reaction_network.RwDict[self.random_walk_type.value].keys()
                print("Have run already")
            except:
                print("NO run")
            # if self.reaction_network.RwDict[self.random_walk_type.value].keys() == None:
                if self.init_simple_random_walk(w=walk_range,l=l,d=d,nmin=nmin, fname=fname) == False:
                    return None
            print(self.reaction_network.RwDict[self.random_walk_type.value].keys())
            print(walk_range)
            return self.reaction_network.RwDict[self.random_walk_type.value].keys()
    
    def get_mak_rw(self, walk_range=range(10),l=10,n=500, trys=10, save=True, fname="rand_mak_walk.json"):
        if self.reaction_network == None:
            return None
        else:
            if self.simple_rw == False:
                if self.init_mak_random_walk(w=walk_range,l=l,n=n,trys=trys,save=save, fname=fname) == False:
                    return None
            print(self.reaction_network.RwDict[self.random_walk_type.value].keys())
            return self.reaction_network.RwDict[self.random_walk_type.value].keys()
            

