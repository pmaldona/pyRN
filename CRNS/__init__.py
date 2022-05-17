from .CRNS import CRNS as crns

def new_from_txt(path):
  reaction_network = crns.form_txt(path)
  return reaction_network
