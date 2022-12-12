import base64
from io import BytesIO
from typing import Any, Union, Callable
import matplotlib.pyplot as plt
from pandas import DataFrame

def plot_stoich(plot: Union[Any, None]) -> Union[str, None]:
    if plot != None:
        #fig, ax = plt.subplots(1, 1, figsize = (10, 5))
        #ax.plot(plot)
        tmpfile = BytesIO()
        # plot.savefig('foo.png')
        plot.savefig(tmpfile, format='png')
        encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')
        
        return encoded
    else:
        return None
    
def plot_df(plot: Union[DataFrame, None], title="") -> Union[str, None]:
    if plot != None:
        fig, ax = plt.subplots(1, 1, figsize = (10, 5))
        ax.plot(plot)
        ax.set_title(title)
        tmpfile = BytesIO()
        fig.savefig(tmpfile, format='png')
        encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

        return encoded
    else:
        return None

def plot_abstractions_callable(callable: Callable, walk_type, walk_index, abstractions_type='c', title='', show_indices=True, index_spacing=1, legend=True) -> str:
    fig, ax = plt.subplots(1, 1, figsize = (10, 5))
    callable(ax, walk_type, walk_index, abstractions_type, title=title, show_indices=show_indices, index_spacing=index_spacing, legend=legend)
    tmpfile = BytesIO()
    fig.savefig(tmpfile, format='png')
    encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

    return encoded

def plot_random_walk_callable(callable: Callable, walk_type, walk_index, abstractions_type='c', title='') -> str:
    fig, ax = plt.subplots(1, 1, figsize = (10, 5))
    callable(ax, walk_type, walk_index, abstractions_type, title)
    tmpfile = BytesIO()
    fig.savefig(tmpfile, format='png')
    encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

    return encoded

def plot_trajectory_callable(callable: Callable, walk_type, conv_pert: bool, walk_index, abstractions_type='c', title='', loga = False) -> str:
    fig, ax = plt.subplots(1, 1, figsize = (10, 5))
    if conv_pert:
        callable(ax, walk_type, walk_index=walk_index, title=title, loga=loga)
    else:
        callable(ax, walk_type, walk_indexes=None, title=title, loga=loga)
    tmpfile = BytesIO()
    fig.savefig(tmpfile, format='png')
    encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

    return encoded

def plot_markov_callable(callable: Callable, walk_type, walk_index = None, abstractions_type='c', save=False, file_path=None) -> str:
    fig, ax = plt.subplots(1, 1, figsize = (18, 7))
    callable(ax, walk_type, walk_indexes=walk_index, abstraction_type=abstractions_type, save=save, file_path=file_path)
    tmpfile = BytesIO()
    fig.savefig(tmpfile, format='png')
    encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

    return encoded

def plot_hist_simple_rw(callable: Callable, walk_type, walk_index = None, abstractions_type='c') -> str:
    fig, ax = plt.subplots(1, 1, figsize = (18, 7))
    callable(ax, walk_type, walk_indexes=walk_index, abstraction_type=abstractions_type)
    tmpfile = BytesIO()
    fig.savefig(tmpfile, format='png')
    encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

    return encoded