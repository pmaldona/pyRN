from itertools  import chain
from matplotlib import ticker
from numpy      import array, polyfit, transpose
from pandas     import concat
from random     import choice

import pickle
import os
import re
import seaborn

### Reoccuring helper functions ###

def generate_random_colors(n: int)->list:
    '''
    Parameters:
        n (int)
    Returns:
        colors (list), a list of random hex-codes for colors
    '''
    return ['#' + ''.join([choice('0123456789ABCDEF') for j in range(6)]) for i in range(n)]

### 1 Single networks, varying maximum perturbation size ###

### 1.1 Loading of dataframes ###

def get_file_paths(path, extension=None):
    """
    Returns a list of file paths in the directory located at `path` and its subdirectories.
    
    Parameters:
        path (string): A string representing the path to the directory to search.
        extension: An optional string representing the file extension to search for. If not provided, all files will be returned.
    Return:
        A list of strings representing the file paths.
    """
    file_paths = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if extension is None or file.endswith(extension):
                file_paths.append(os.path.join(root, file))
    return file_paths

def depkl(path):
    '''
    Load python object from given path using pickle
    
    Parameters:
        path   (string), filepath
        
    Returns:
        object (python object)
    '''
    with open(path, 'rb') as f:
        object = pickle.load(f)
    return object

def sort_strings_by_number(strings: list)->list:
    '''
    Parameters:
        strings (list), a list of strings
    Returns:
        strings (list), the same list but sorted by the last number occuring in the string
    Ensures that abstraction dataframes and transition dataframes will be sorted by maximum perturbation sizes
    even if they are >9
    '''
    return sorted(strings, key=lambda x: int(re.findall("\\d+", x)[0]))

def get_abstractions_dfs_file_paths(path, min_max_pert_size=1, max_max_pert_size=0):
    '''
    Parameters:
        path (str),          path to a folder with abstractions_dfs for a single reaction networks
        max_pert_size (int), limit for maximum pert size (no limit if 0)
    Returns:
        abstractions_dfs (list), a list of the file paths of to abstractions dataframes stored at path where the maximum perturbation
                                 size is smaller than max_pert_size if max_pert_size != 0
    '''
    abstractions_dfs_file_paths = [file_path for file_path in get_file_paths(path, extension=None) if 'abstractions_df' in file_path]
    abstractions_dfs_file_paths = sort_strings_by_number(abstractions_dfs_file_paths)
    if not max_max_pert_size:
        return abstractions_dfs_file_paths
    return abstractions_dfs_file_paths[min_max_pert_size-1:max_max_pert_size]

def get_transtions_dfs_file_paths(path, min_max_pert_size=1, max_max_pert_size=0):
    '''
    Parameters:
        path (str),          path to a folder with transitions_dfs for a single reaction networks
        max_pert_size (int), limit for maximum pert size (no limit if 0)
    Returns:
        abstractions_dfs (list), a list of the file paths of transitions dataframes stored at path where the maximum perturbation
                                 size is smaller than max_pert_size if max_pert_size != 0
    '''
    transtions_dfs_file_paths = [file_path for file_path in get_file_paths(path, extension=None) if 'transitions_df' in file_path]
    transtions_dfs_file_paths = sort_strings_by_number(transtions_dfs_file_paths)
    if not max_max_pert_size:
        return transtions_dfs_file_paths
    return transtions_dfs_file_paths[min_max_pert_size-1:max_max_pert_size]

def get_abstractions_dfs(path, min_max_pert_size=1, max_max_pert_size=0):
    '''
    Parameters:
        path (str),          path to a folder with abstractions_dfs for a single reaction networks
        max_pert_size (int), limit for maximum pert size (no limit if 0)
    Returns:
        abstractions_dfs (list), a list of the abstractions dataframes stored at path where the maximum perturbation
                                 size is smaller than max_pert_size if max_pert_size != 0
    '''
    abstractions_dfs_file_paths = get_abstractions_dfs_file_paths(path, min_max_pert_size, max_max_pert_size)
    return [depkl(path) for path in abstractions_dfs_file_paths]

def get_transitions_dfs(path: str, min_pert_size=1, max_pert_size=0)->list:
    '''
    Parameters:
        path (str),          path to a folder with transitions_dfs for a single reaction networks
        max_pert_size (int), limit for maximum pert size (no limit if 0)
    Returns:
        abstractions_dfs (list), a list of the transitions dataframes stored at path where the maximum perturbation
                                 size is smaller than max_pert_size if max_pert_size != 0
    '''
    transitions_dfs_file_paths   = get_transtions_dfs_file_paths(path, min_pert_size, max_pert_size)
    return [depkl(path) for path in transitions_dfs_file_paths]

### 1.2. abstractions_df plots ###

def cumulate(lists):
    '''
    Parameter:
        lists (list of lists)

    Returns
        lists (list of lists)
    '''
    for i in range(1, len(lists)):
        lists[i] = [lists[i][j]+lists[i-1][j] for j in range(len(lists[i]))]
    return lists

def plot_series_in_cumulative_plot(subplot, steps, series, color, label):
    subplot.plot(steps, series, c='black')
    subplot.fill_between(steps, series, [0 for step in series], color=color, label=label)

def cumulative_plot(subplot, steps, series, colors, labels, legend=False):
    for i, s in enumerate(reversed(series)):
        plot_series_in_cumulative_plot(subplot, steps, s, color=colors[len(series)-i-1], label=labels[len(series)-i-1])

def plot_global_resilience(subplot, abstractions_dfs, colors, legend=False)->None:
    '''
    Parameters:
        subplot (matplotlib.axes._subplots.AxesSubplot)
        abstractions_dfs (list), a list of abstractions dataframes with increasing maximum perturbation size
        colors (list), a list of strings representing colors
        legend (boolean), display a legend?
    x-axis: maximum perturbation size
    y-axis: cumulated global resilience
    '''
    global_resiliences = list(transpose([abstractions_df['global_resilience'] for abstractions_df in abstractions_dfs]))
    global_resiliences = cumulate(global_resiliences)
    labels=abstractions_dfs[0]['abstraction']
    cumulative_plot(subplot, range(1, len(global_resiliences[0])+1), global_resiliences, colors=colors, labels=labels)
    subplot.set_xlabel('Maximum perturbation size', fontsize=20)
    subplot.set_ylabel('Global resilience', fontsize=20)
    subplot.set_title('Global resilience', fontsize=20)
    subplot.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    subplot.tick_params(labelsize=20)
    if legend:
        subplot.legend()

def plot_local_resilience(subplot, abstractions_dfs, colors, legend=False):
    '''
    Parameters:
        subplot (matplotlib.axes._subplots.AxesSubplot)
        abstractions_dfs (list), a list of abstractions dataframes with increasing maximum perturbation size
        colors (list), a list of strings representing colors
    x-axis: maximum perturbation size
    y-axis: local resilience
    '''
    local_resiliences = list(transpose([abstractions_df['local_resilience'] for abstractions_df in abstractions_dfs]))
    #colors = ['#ffff00','#ff00ff','#00ffff','#ff0000','#00ff00','#0000ff','#000000']
    labels=abstractions_dfs[0]['abstraction']
    for i, local_resilience in enumerate(local_resiliences):
        subplot.scatter(range(1, len(local_resilience)+1),
                        local_resilience,
                        label=list(labels[i]),
                        color=colors[i])
        subplot.plot(range(1, len(local_resilience)+1),
                     local_resilience,
                     color=colors[i])
    subplot.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    subplot.tick_params(labelsize=20)
    subplot.set_xlabel('Maximum perturbation size', fontsize=20)
    subplot.set_ylabel('Local resilience', fontsize=20)
    subplot.set_ylim(ymin=0)
    subplot.set_title('Local resilience', fontsize=20)
    if legend:
        subplot.legend()

def plot_local_resilience_vs_global(subplot, abstractions_dfs, colors, legend=False):
    '''
    Parameters:
        subplot (matplotlib.axes._subplots.AxesSubplot)
        abstractions_dfs (list), a list of abstractions dataframes with increasing maximum perturbation size
        colors (list), a list of strings representing colors
    x-axis: local resilience
    y-axis: global resilience
    '''
    local_resiliences = list(transpose([abstractions_df['local_resilience'] for abstractions_df in abstractions_dfs]))
    global_resiliences = list(transpose([abstractions_df['global_resilience'] for abstractions_df in abstractions_dfs]))
    labels = abstractions_dfs[0]['abstraction']
    for i in range(len(local_resiliences)):
        subplot.scatter(local_resiliences[i], global_resiliences[i], label=list(labels[i]), color=colors[i], s=[x*10 for x in range(len(local_resiliences[i]))])
        subplot.plot(local_resiliences[i], global_resiliences[i], color=colors[i])
    subplot.tick_params(labelsize=20)
    subplot.set_xlabel('Local resiliences', fontsize=20)
    subplot.set_ylabel('Global resiliences', fontsize=20)
    subplot.set_title('Local resilience vs global resilience', fontsize=20)
    if legend:
        subplot.legend()

def plot_org_size_vs_resilience(subplot, abstractions_dfs, colors, type='local'):
    '''
    Parameters:
        subplot (matplotlib.axes._subplots.AxesSubplot)
        abstractions_dfs (list),   a list of abstractions dataframes with increasing maximum perturbation size
        colors           (list),   a list of strings representing colors
        type             (string), 'local' or 'global'
    x-axis: organization size
    y-axis: 
    '''
    subplot.set_ylim(ymin=0, ymax=1)
    subplot.tick_params(labelsize=20)
    subplot.set_xlabel('number of species in organization', fontsize=20)
    subplot.set_xticks(list(range(0,len(abstractions_dfs[0]['abstraction'].iloc[0]), 5)))
    subplot.set_ylabel('local resilience', fontsize=20)
    for i, a_df in enumerate(abstractions_dfs):
        subplot.scatter(a_df[['number of species']], a_df[[type+'_resilience']], color=colors[i], label=f'max_pert-size: {i+1}', s=100)
        # Linear regression
        coefficients = polyfit(a_df['number of species'], a_df[type+'_resilience'], 1)
        subplot.plot(a_df['number of species'], coefficients[0]*array(a_df['number of species']) + coefficients[1], color=colors[i])
    subplot.legend(fontsize=20)

### 1.3 Transitions ###

def plot_number_of_transitions(subplot, transitions_dfs):
    '''
    Parameters:
        subplot (matplotlib.axes._subplots.AxesSubplot),
        transitions_dfs (list), list of transitions dataframes
    x-axis: maximum perturbation size
    y-axis: transitions with probability > 0
    '''
    number_of_transitions = [transitions_df.loc[transitions_df['probability']>0].shape[0] for transitions_df in transitions_dfs]
    subplot.scatter(range(1,len(number_of_transitions)+1), number_of_transitions, s=100)
    subplot.plot(range(1,len(number_of_transitions)+1), number_of_transitions)
    subplot.set_xlabel('Maximum perturbation size', fontsize=20)
    subplot.set_ylabel('Number of transitions', fontsize=20)
    subplot.set_ylim(ymin=0)
    subplot.set_title('Number of possible transitions', fontsize=20)
    subplot.grid(False)

def plot_transition_probabilities(subplot, transitions_dfs, colors):
    '''
    Parameters:
        subplot (matplotlib.axes._subplots.AxesSubplot),
        transitions_dfs (list), list of transitions dataframes
        colors (list), list of strings representing colors
    x-axis: maximum perturbation size
    y-axis: transitions probability
    '''
    probabilities = transpose([transitions_df['probability'] for transitions_df in transitions_dfs])
    maximum_probabilities = [transitions_df['probability'].max() for transitions_df in transitions_dfs]
    minimal_probabilities = [transitions_df['probability'].min() for transitions_df in transitions_dfs]
    subplot.fill_between(range(1,len(minimal_probabilities)+1), minimal_probabilities, maximum_probabilities, color='#0000ff', alpha=0.1)
    for i,p in enumerate(probabilities):
        subplot.plot(range(1,len(p)+1),p, color=colors[i])
    subplot.set_xlabel('Maximum perturbation size', fontsize=20)
    subplot.set_ylabel('Transition probabilities', fontsize=20)
    subplot.set_title('Transition probabilities', fontsize=20)
    subplot.set_ylim(ymin=0)

def plot_transition_size_vs_probabiliy(subplot, transitions_dfs, colors):
    '''
    Parameters:
        subplot (matplotlib.axes._subplots.AxesSubplot),
        transitions_dfs (list), list of transitions dataframes
        colors (list), list of strings representing colors
    x-axis: maximum perturbation size
    y-axis: transitions size (net change of species)
    '''
    sizes         = transpose([transitions_df['size_difference'] for transitions_df in transitions_dfs])
    probabilities = transpose([transitions_df['probability'] for transitions_df in transitions_dfs])
    for i in range(len(sizes)):
        subplot.scatter(probabilities[i], sizes[i], color=colors[i], s=[x*20 for x in range(len(sizes[i]))])
        subplot.plot(probabilities[i], sizes[i], color=colors[i])
    subplot.set_xlabel('Transition probabilities', fontsize=20)
    subplot.set_ylabel('Difference in number of species', fontsize=20)
    subplot.set_title('Transition probabilities by size', fontsize=20)

### 2. Multiple networks ###

def get_abstractions_dfs_for_multiple_networks(path, min_max_pert_size=1, max_max_pert_size=0):
    '''
    Parameters:
    '''
    return list(chain(*[get_abstractions_dfs(os.path.join(path, folder), min_max_pert_size, max_max_pert_size) for folder in os.listdir(path)]))

def binning_abstraction_dfs_by_size(abstractions_dfs: list, bins: list)->list:
    '''
    Parameters:
        abstractions_dfs (list), list of abstractions DataFrames
        bins             (list), list of tuples, each indicating the minimal and maximal size for the bins
    Returns:
        abstractions_dfs grouped by their size according to bins
    '''
    return [[abstractions_df for abstractions_df in abstractions_dfs if bin[0]<=len(abstractions_df['abstraction'])<=bins[1]] for bin in bins]

def binning_abstraction_dfs_by_column_value(abstractions_dfs: list, bins: list, column)->list:
    '''
    Parameters:
        abstractions_dfs (list), list of abstractions DataFrames
        bins             (list), list of tuples, each indicating the minimal and maximal size for the bins
        column           (str),  name of the column
    Returns:
        abstractions_dfs grouped by their column value according to bins
    '''
    return [[abstractions_df for abstractions_df in abstractions_dfs if bin[0]<=abstractions_df[column]<=bins[1]] for bin in bins]

def pairwise_plot(abstractions_dfs, column_names):
    '''
    Parameters:
        abstractions_dfs (list), list of abstractions DataFrames
        column_names     (list), list of the column names
    Returns:
        A nice seaborn pairplot
    '''
    plot_df = concat(abstractions_dfs, ignore_index=True)[column_names]
    g = seaborn.pairplot(plot_df, height=3)
    g.map_upper(seaborn.kdeplot, fill=True)
    g.map_upper(seaborn.regplot, line_kws={'color':'red'}, ci=None)
    g.set(xlim=(0,None), ylim=(0,None))
    g.fig.suptitle('Pairwise plot for the variables in the abstractions_dfs\nwith linear regression and density estimation in the right upper half', fontsize=20)
    g.fig.subplots_adjust(top=.925)
    return g

def correlations_heatmap(abstractions_dfs, column_names):
    '''
    Parameters:
        abstractions_dfs (list), list of abstractions DataFrames
        column_names     (list), list of the column names
    Returns:
        A heatmap of correlations
    '''
    plot_df = concat(abstractions_dfs, ignore_index=True)[column_names]
    seaborn.set(rc={'figure.figsize':(10,10)})
    h = seaborn.heatmap(plot_df.corr(), annot=True, vmin=0, vmax=1, cmap='Reds', linewidth=.5, square=True)
    h.set_title('Heatmap of correlations', fontsize=20)
    h.set_xticklabels(h.get_xticklabels(), fontsize=15, rotation=45)
    h.set_yticklabels(h.get_yticklabels(), fontsize=15, rotation=45)
    return h

def plot_correlation_curves(abstractions_dfs: list, column_names: list, axes, network_sizes):
    '''
    Parameters:
        abstractions_dfs (list), list of lists of abstractions DataFrames grouped by some property
        column_names     (list), list of the column names
        axes             
    Returns:
        A heatmap of correlations
    '''
    for i, column_name_1 in enumerate(column_names):
        for j, column_name_2 in enumerate(column_names[i+1:]):
            correlations = [abstractions_dfs[k][column_name_1].corr(abstractions_dfs[k][column_name_2]) for k in range(len(abstractions_dfs))]
            axes.plot(network_sizes, correlations, label=f'{column_name_1} , {column_name_2}', lw=3)
    
    axes.set_ylabel('Network sizes', fontsize=20)
    axes.set_ylabel('Correlations',  fontsize=20)
    axes.set_ylim(ymin=-1, ymax=1)
    axes.tick_params(labelsize=20)
    axes.set_xticks(network_sizes)
    axes.legend(fontsize=15)