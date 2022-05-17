from .CRNS import CRNS as crns

def from_txt(path):
  reaction_network = crns.form_txt(path)
  return reaction_network
