import matplotlib.pyplot as plt
import numpy as np

def plot_raw(subplot, abstractions, title=''):
    m = np.transpose(np.array(abstractions))
    subplot.imshow(m,vmin=0, vmax=1)
    subplot.set_xlabel('Iterations')
    subplot.set_ylabel('Species')
    subplot.set_title(title+'\nRaw Plot')
