from pyRN.SEA import sos
import hasseNetworkx

from bitarray import frozenbitarray
import matplotlib.pyplot
import matplotlib.lines
import networkx
import numpy
import pickle

def is_subset_of(abstraction_1, abstraction_2):
    ''' 
    Evaluates if abstraction_1 is a subset of abstraction_2

        Parameters:
            abstraction_1 (string), 
            abstraction_1 (string),
            Both binary strings with each digit indicating either
            the presence (1) or the absence of a species (0)

        Returns:
        (boolean), True if abstraction_1 is a subset of abstraction_2

    Examples:
    is_subset_of('0110', '0010') -> True
    is_subset_of('1100', '0010') -> False
    '''
    set_1 = list(frozenbitarray(abstraction_1))
    set_2 = list(frozenbitarray(abstraction_2))
    return sos.is_subset_of(set_1, set_2)

def set_size(abstraction):
    '''
    Determines the number of species of an abstraction

        Parameter:
            abstraction (string),
            A binary string with each digit indicating either
            the presence (1) or the absence of a species (0)

        Returns:
            (int), the number of species present in abstraction
    '''
    return sos.n_elements(list(frozenbitarray(abstraction)))

def create_subset_graph(abstractions_df):
    '''
    Creates a graph where the nodes are the abstractions stored in abstractions_df
    and the edges are the subset relationship between said abstractions.
    Note that the Graph object uses a binary string representation of the abstractions.

        Parameters:
           abstractions_df (pandas.DataFrame)
           transitions_df  (pandas.DataFrame)

        Returns:
            subset_graph (networkx.DiGraph)
    '''
    nodes = [abstraction.to01() for abstraction in abstractions_df['abstraction']]
    edges = [(node_1, node_2) for node_1 in nodes for node_2 in nodes if ((node_1 != node_2) and is_subset_of(node_1, node_2))]
    # Defining the Graph object
    subset_graph = networkx.DiGraph()
    subset_graph.add_nodes_from(nodes)
    subset_graph.add_edges_from(edges)
    # Editing of the graph object (Removing edges implied by transitivity)
    hasseNetworkx.transitivity_elimination(subset_graph)
    return subset_graph

def plot_subset_relationships(subplot, subset_graph, node_positions, node_size):
    '''
    A function for the plotting of the subset_graph

        Parameters:
           subplot        (matplotlib.pyplot subplot),
           subset_graph   (networkx.DiGraph),
           node_positions (dictionary),
           node_size      (scalar or array)

        Returns:
            (matplotlib.lines.Line2D), a handle to use in the legend to indicate the subset relation_ships    
    '''
    networkx.draw_networkx_edges(subset_graph,
                                 pos=node_positions,
                                 edge_color="#0000ff",
                                 node_size=node_size,
                                 arrowstyle='-',
                                 ax=subplot)
    return matplotlib.lines.Line2D(range(3), range(3), color="white", marker="_", mec="#0000ff", markersize=30)

def create_transitions_graph(abstractions_df, transitions_df):
    '''
    Creates a graph where the nodes are the abstractions stored in abstractions_df
    and the edges are the subset relationship between said abstractions.
    Note that the Graph object uses a binary string representation of the abstractions.

        Parameters:
           abstractions_df (pandas.DataFrame),
           transitions_df  (pandas.DataFrame)

        Returns:
            transitions_graph (networkx.Graph object)
    '''
    nodes = [abstraction.to01() for abstraction in abstractions_df['abstraction']]
    edges = [(transition[0].to01(),transition[1].to01()) for transition in zip(transitions_df['initial_state'], transitions_df['convergent_state'])]
    transitions_graph = networkx.DiGraph()
    transitions_graph.add_nodes_from(nodes)
    transitions_graph.add_edges_from(edges)
    return transitions_graph

def plot_transitions(subplot, transitions_graph, transitions_df, node_positions, node_size):
    '''
    A function for the plotting of the transitions_graph

        Parameters:
           subplot        (matplotlib.pyplot subplot),
           subset_graph   (networkx.DiGraph),
           node_positions (dictionary),
           node_size      (scalar or array)

        Returns:
            (matplotlib.lines.Line2D), a handle to use in the legend to indicate the subset relation_ships    
    '''
    probabilities = [transitions_df.loc[(transitions_df['initial_state']==frozenbitarray(transition[0])) &
                                        (transitions_df['convergent_state']==frozenbitarray(transition[1]))]['probability']
                     for transition in transitions_graph.edges()]
    networkx.draw_networkx_edges(transitions_graph,
                            pos = node_positions,
                            edge_color="#ff0000",
                            node_size=node_size,
                            connectionstyle='arc3, rad = 0.2',
                            arrowstyle='->',
                            width=numpy.array(probabilities)*5,
                            ax=subplot)
    return matplotlib.lines.Line2D(range(3), range(3), color="white", marker="_", mec="#ff0000", markersize=30)

def get_node_labels(abstractions_df, index_labels):
    '''
    Returns the labels for the nodes of the markov plot

        Parameters:
            abstractions_df (pandas.DataFrame), dataframe with information on the organizations
            index_labes     (boolean),          if True: labels will correspond to the index of the organizations in the dataframe
                                                else:    labels will correspond to the bitstring representation of the organizations
        
        Returns:
            (dictionary)
    '''
    if index_labels:
        return {abstraction.to01(): int(abstractions_df.index[abstractions_df['abstraction']==abstraction][0]) for abstraction in abstractions_df['abstraction']}
    else:
        return {abstraction.to01(): abstraction.to01() for abstraction in abstractions_df['abstraction']}
    
def plot_markov(subplot, abstractions_df, transitions_df, shift_x=False, node_size='global_resilience', index_labels=False):
    '''
    A function for the plotting of the complete markov_graph

        Parameters:
           subplot,
           abstractions_df (pandas.DataFrame),
           transitions_df  (pandas.DataFrame),
           shift_x         (boolean)
           node_size       (string), key of the abstractions_df (default='global_resilience)
           index_labels    (boolean), if True: node labels will correspond to the index of the organizations in the dataframe
                                      else:    labels will correspond to the bitstring representation of the organizations (default)

        Returns:
            (matplotlib.lines.Line2D), a handle to use in the legend to indicate the subset relation_ships    
    '''
    # Create Graphs
    subset_graph      = create_subset_graph(abstractions_df)
    transitions_graph = create_transitions_graph(abstractions_df, transitions_df)
    # Prepare plotting
    node_positions = hasseNetworkx.layout(subset_graph, set_size, shift_x=shift_x)
    node_sizes = [float(abstractions_df.loc[abstractions_df['abstraction']==frozenbitarray(node)][node_size].iloc[0]) for node in subset_graph.nodes()]
    node_sizes = [node_size*5000/max(node_sizes) for node_size in node_sizes]
    node_labels = get_node_labels(abstractions_df, index_labels)
    textbox = dict(boxstyle='round', edgecolor='black', facecolor='none', alpha=1, lw=2)
    # Plot
    layer_1 = plot_subset_relationships(subplot, subset_graph, node_positions, node_size=node_sizes)
    layer_2 = plot_transitions(subplot, transitions_graph, transitions_df, node_positions, node_size=node_sizes)
    layer_3 = matplotlib.lines.Line2D(range(5), range(5), color="white", marker='o', markerfacecolor="#00ff00")
    networkx.draw_networkx_nodes(subset_graph, node_positions, node_size=node_sizes, node_color='#00ff00')
    networkx.draw_networkx_labels(subset_graph, node_positions, node_labels, font_size=15, bbox=textbox)
    matplotlib.pyplot.legend((layer_1,layer_2, layer_3),("Subset relationship","Transition (width~probability)", node_size), fontsize=15)
    subplot.tick_params(left=True, labelleft=True)
    subplot.set_ylabel('Number of species', fontsize=15)