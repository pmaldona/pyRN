import matplotlib.pyplot
import networkx
import pandas
import sos

def nodes(G, abstractions_df):
    '''
    In:  G, Empyty directed multigraph networkx object
         abstractions_df (pandas data frame)
    Out: G with nodes from abstractions_df added
    '''

    abstractions = [sos.from_string(a) for a in list(abstractions_df['a'])]

    bins = [[] for i in range(len(abstractions[0])+1)]
    for i,abstraction in enumerate(abstractions):
        bins[sos.n_elements(abstraction)].append(abstraction)

    max_bin_size = 1

    for bin in bins:
        if len(bin) > max_bin_size:
            max_bin_size = len(bin)

    shift = 1
    for i,bin in enumerate(bins):
        if (i>0) and (len(bins[i])==len(bins[i-1])):
            shift = 0.5
        for j,abstraction in enumerate(bin):       
            G.add_node(list(abstractions_df.index[abstractions_df['a']==str(abstraction)])[0],
                pos=(max_bin_size-(max_bin_size/(len(bin)+1)*(j+shift)),sos.n_elements(abstraction)))
        shift = 1
    return G

def proper_subset_relationship_matrix(abstractions):

    n = len(abstractions)

    M = [[0 for j in range(n)] for i in range(n)]

    for i in range(n):
        for j in range(i+1,n):
            if sos.is_subset_of(abstractions[i],abstractions[j]):
                M[i][j] = 1
    
    return M

def transitivity_elimination(relationship_matrix):

    n = len(relationship_matrix)

    for i in range(n):
        for j in range(n-1,0,-1):
            if relationship_matrix[i][j] == 1:
                for i_ in range(n):
                    if relationship_matrix[i_][j] == 1 and relationship_matrix[i][i_] == 1:
                        relationship_matrix[i][j] = 0

    return relationship_matrix

def subset_relationships(G, abstractions_df):

    abstractions = [sos.from_string(a) for a in list(abstractions_df['a'])]

    proper_subset_matrix = proper_subset_relationship_matrix(abstractions)

    proper_subset_matrix = transitivity_elimination(proper_subset_matrix)

    n = len(abstractions)
    for i in range(n):
        for j in range(i+1, n):
            if proper_subset_matrix[i][j] == 1:
                G.add_edge(i, j, type='subset relationship')
    
    return G

def transition_relations(G, transitions_df):

    for i in range(transitions_df.shape[0]):
        G.add_edge(transitions_df.loc[i,'a1'], transitions_df.loc[i,'a2'], type='transition', probability=transitions_df.loc[i, 'probability'])
    
    return G

def create_graph(abstractions_df, transitions_df):

    # Initialize Graph object
    G = networkx.MultiDiGraph()

    # Initialize nodes and position them on the grid in a hasse-diagramm
    G = nodes(G, abstractions_df)

    # Add edges for subset relations
    G = subset_relationships(G, abstractions_df)

    # Add edges for transitions
    G = transition_relations(G, transitions_df)

    return G

def draw_nodes(G, abstractions_df, ax):

    print('Node legend:')
    print(abstractions_df['a'])

    pos=networkx.get_node_attributes(G,'pos')
    networkx.draw_networkx(G, pos, edgelist=[],node_size=0, ax=ax)

    nodes = list(G.nodes())

    maintainabilities = [abstractions_df.loc[n, 'maintainability'] for n in nodes]
    reachabilities    = [abstractions_df.loc[n, 'reachability'] for n in nodes]
    scaled_maintainability  = [maintainabilities[i]/max(maintainabilities)*2500 for i in range(len(maintainabilities))]
    scaled_reachability     = [reachabilities[i]/max(reachabilities)*2500 for i in range(len(reachabilities))]
    scaled_reachability     = [2500-r for r in scaled_reachability]

    mb = [nodes[i] for i in range(len(nodes)) if (scaled_maintainability[i]>=scaled_reachability[i])]
    rb = [nodes[i] for i in range(len(nodes)) if (scaled_maintainability[i]<=scaled_reachability[i])]
    ms = [nodes[i] for i in range(len(nodes)) if (scaled_maintainability[i]<=scaled_reachability[i])]
    rs = [nodes[i] for i in range(len(nodes)) if (scaled_maintainability[i]>=scaled_reachability[i])]

    mbs = [scaled_maintainability[i] for i in range(len(nodes)) if (scaled_maintainability[i]>=scaled_reachability[i])]
    rbs = [scaled_reachability[i] for i in range(len(nodes)) if (scaled_maintainability[i]<=scaled_reachability[i])]
    mss = [scaled_maintainability[i] for i in range(len(nodes)) if (scaled_maintainability[i]<=scaled_reachability[i])]
    rss = [scaled_reachability[i] for i in range(len(nodes)) if (scaled_maintainability[i]>=scaled_reachability[i])]

    networkx.draw_networkx_nodes(G, pos, nodelist=mb, node_size=mbs, node_color="#b000b0", alpha=0.8, ax=ax)
    networkx.draw_networkx_nodes(G, pos, nodelist=rb, node_size=rbs, node_color="#00ff00", alpha=0.8, ax=ax)
    networkx.draw_networkx_nodes(G, pos, nodelist=rs, node_size=rss, node_color="#00ff00", alpha=0.8, ax=ax)
    networkx.draw_networkx_nodes(G, pos, nodelist=ms, node_size=mss, node_color="#b000b0", alpha=0.8, ax=ax)

    line1 = matplotlib.lines.Line2D(range(5), range(5), color="white", marker='o', markerfacecolor="#00ff00")
    line2 = matplotlib.lines.Line2D(range(5), range(5), color="white", marker='o', markerfacecolor="#b000b0")

    return (line1, line2)

def draw_subset_relationships(G, ax):

    pos=networkx.get_node_attributes(G,'pos')
    subset_relationships = [e for e in list(G.edges(data=True)) if e[2].get('type')=='subset relationship']
    networkx.draw_networkx_edges(G, pos, edgelist=subset_relationships, edge_color="#0000ff", node_size=0, label='subset_relations', arrowstyle='-', ax=ax)
    return matplotlib.lines.Line2D(range(3), range(3), color="white", marker="_", mec="#0000ff", markersize=20)

def draw_transition_probabilities(G, ax):

    pos=networkx.get_node_attributes(G,'pos')
    transitions = [e for e in list(G.edges(data=True)) if e[2].get('type')=='transition']
    probabilities = [e[2].get('probability')*5 for e in transitions]
    networkx.draw_networkx_edges(G, pos, edgelist=transitions, edge_color="#ff0000", node_size=2500, connectionstyle='arc3, rad = 0.05', width=probabilities, arrowstyle='->', ax=ax)
    return matplotlib.lines.Line2D(range(3), range(3), color="white", marker="_", mec="#ff0000", markersize=20)

def draw_graph(G, abstractions_df):
    
    fig, ax = matplotlib.pyplot.subplots()

    # Draw nodes with markov properties
    layer1, layer2 = draw_nodes(G, abstractions_df)

    # Draw subset relationships
    layer3 = draw_subset_relationships(G)

    # Draw transitions with probability
    layer4 = draw_transition_probabilities(G)
    
    # Display polt with matplolib.pyplot
    ax.tick_params(left=True, labelleft=True)
    ax.set_ylabel('Number of species')
    matplotlib.pyplot.legend((layer1,layer2,layer3, layer4),("Reachability","Maintainability","Subset relationship","Transition (width~probability)"))
    matplotlib.pyplot.show()
