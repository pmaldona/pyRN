from typing import Union
from pyRN import pyRN
from ..App.App import state
from pyvis.network import Network

def get_network_from_set(species_set) -> Union[Network, None]:
    network = state.get_network()
    if network == None:
        return None
    else:
        print(network.get_nodes())
        return Network()
