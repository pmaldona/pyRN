from itertools  import combinations, product
from matplotlib import pyplot, gridspec
from networkx   import *
from numpy      import array, matmul, mean, round
from pandas     import DataFrame, options
from random     import choice

from pyRN.SEA.markov import *

def transition_graph_from_dataframes(abstractions_df, transitions_df):
    '''
    Parameters:
        abstractions_df (pandas.DataFrame), 
        transitions_df  (pandas.DataFrame)
    Returns:
        G (networkx.DiGraph), a directed graph of the transitions
    '''
    T = transition_matrix_from_dataframes(abstractions_df, transitions_df)
    G = DiGraph()
    G.add_edges_from([(j,i) for i in range(abstractions_df.shape[0]) for j in range(abstractions_df.shape[0]) if T[i][j]])
    return G

####################################################################################################
##############################      GLOBAL COLLECTIVE RESILIENCE      ##############################
####################################################################################################

def global_collective_resilience(abstractions_df, organization_indices):
    '''
    Parameters:
        abstractions_df,     (pandas.DataFrame)
        organization_indices (list), a list of indices refering to organizations in the abstractions_df
    Returns:
        The global collective resilience of the set of organizations corresponding to the indices
        in the organization_indices list as the sum of individual global resiliences of the organizations
    '''
    return sum([abstractions_df.iloc[i]['global_resilience'] for i in organization_indices])

####################################################################################################
##############################      LOCAL COLLECTIVE RESILIENCE       ##############################
####################################################################################################

##############################   CONDITIONAL PROBABILITY ESTIMATION    #############################

#                                           ESTIMATOR 1                                            #

def conditional_probability_estimator_1(abstractions_df, transitions_df, organization_indices):
    '''
    Parameters:
        abstractions_df,
        transitions_df,
        organization_indices
    Returns:
        estimation (list), a probability vector 
    Estimation of p(o|O) for all organizations.
    Assumes uniform distribution over all organizations in the set
    '''
    try:
        return [1/len(organization_indices) if index in organization_indices else 0 for index in range(abstractions_df.shape[0])]
    except:
        return [0 for index in range(abstractions_df.shape[0])]

#                                           ESTIMATOR 2                                            #
    
def readjust_probabilities(transitions_df):
    '''
    In a transitions_df probabilities are scaled in a way such that
    For every organization the sum of transitions starting from it is 1
    '''
    organizations = list(set(transitions_df['initial_state']))
    for start, start_organization in enumerate(organizations):
        try: 
            scaling_factor = 1/sum(transitions_df.loc[(transitions_df['initial_state']==start_organization)]['probability'])
        except:
            scaling_factor = 0
        for end, end_organization in enumerate(organizations):
            transitions_df.loc[(transitions_df['initial_state']==start_organization) & (transitions_df['convergent_state']==end_organization), 'probability'] *= scaling_factor
    return transitions_df

def block_leaving_transitions(abstractions_df, transitions_df, organization_indices):
    '''
    In a transitions_df transition probabilities of organizations with indices from
    organization_indices (refering to the abstractions_df) are set to 0
    '''
    for start, start_organization in enumerate(list(abstractions_df['abstraction'])):
        # Setting probabilities of transitions leaving the set to 0
        for end, end_organization in enumerate(list(abstractions_df['abstraction'])):
            if (start in organization_indices) and (not end in organization_indices):
                transitions_df.loc[(transitions_df['initial_state']==start_organization) & (transitions_df['convergent_state']==end_organization), 'probability'] = 0
    return transitions_df

def make_set_unleavable(abstractions_df, transitions_df, organization_indices):
    '''
    Returns: new_transitions_df
    Transitions starting from organizations with indices from organization_indices
    can only end in organizations with indices from organization_indices.
    Probabilities will be adapted accordingly.
    '''
    new_transitions_df = transitions_df.copy(deep=True)
    new_transitions_df = block_leaving_transitions(abstractions_df, new_transitions_df, organization_indices)
    new_transitions_df = readjust_probabilities(new_transitions_df)
    return new_transitions_df

def conditional_probability_estimator_2(abstractions_df, transitions_df, organization_indices):
    '''
    Parameters:
        abstractions_df,
        transitions_df,
        organization_indices
    Returns:
        estimation (list), a probability vector 
    Estimation of p(o|O) for all organizations.
    Creates a changed version of the original Markov-model that
    excludes the complement of the set of organization and returns
    its stationary distribution
    '''

    # Adapted stationary transition matrix
    T = make_set_unleavable(abstractions_df, transitions_df, organization_indices)
    T = T_inf(transition_matrix_from_dataframes(abstractions_df, T), max_steps=100)
    # Initial uniform distribution
    p = [1/len(organization_indices) if i in organization_indices else 0 for i in range(abstractions_df.shape[0])]
    # Stationary distribution
    s = matmul(T, p)
    s = [round(x, 10) for x in s]

    return s

#                                           ESTIMATOR 3                                            #

def conditional_probability_estimator_3(abstractions_df, transitions_df, organization_indices):
    '''
    Parameters:
        abstractions_df,
        transitions_df,
        organization_indices
    Returns:
        estimation (list), a probability vector 
    Estimation of p(o|O) for all organizations.
    Rescales the global resiliences of the organizations to 1
    '''
    s = sum([gr for i, gr in enumerate(abstractions_df['global_resilience'].tolist()) if i in organization_indices])
    if s:
        return [r/s if i in organization_indices else 0 for i, r in enumerate(abstractions_df['global_resilience'])]
    else:
        return [0 for index in range(abstractions_df.shape[0])]
    

############################   LOCAL COLLECTIVE RESILIENCE ESTIMATION    ###########################
    

def local_collective_resilience(abstractions_df, transitions_df, organization_indices, estimator):
    '''
    Parameters:
        abstractions_df      (pandas.DataFrame),
        transitions_df       (pandas.DataFrame),
        organization_indices (list)
        estimator            (int)
    Returns:
        local_collective_resilience (float)
    Computes an estimated value of the local resilience using the chosen estimator
    '''
    if estimator == 1:
        s = conditional_probability_estimator_1(abstractions_df, transitions_df, organization_indices)
    if estimator == 2:
        s = conditional_probability_estimator_2(abstractions_df, transitions_df, organization_indices)
    if estimator == 3:
        s = conditional_probability_estimator_3(abstractions_df, transitions_df, organization_indices)

    lcr = 0
    for i, org_1 in enumerate(organization_indices):
        for j, org_2 in enumerate(organization_indices):
            start = abstractions_df.iloc[org_1]['abstraction']
            end   = abstractions_df.iloc[org_2]['abstraction']
            p     = transitions_df.loc[(transitions_df['initial_state']==start) & (transitions_df['convergent_state']==end)]['probability'].item()*s[organization_indices[i]]
            lcr  += p
    return round(lcr,10)

####################################################################################################
##############################        TWO STATE MARKOV_MODELS         ##############################
####################################################################################################

def estimated_transitions_matrix(abstractions_df, transitions_df, organization_indices, estimator):
    '''
    Parameters:
        abstractions_df      (pandas.DataFrame),
        transitions_df       (pandas.DataFrame),
        organization_indices (list),
        estimator            (int),
    Returns:
        T (2D-array)
        The transitions matrix of the Markov-model for O and Ō
        based on estimate local collective resiliences
    '''

    # Determine Ō
    complement_indices = [i for i in range(abstractions_df.shape[0]) if i not in organization_indices]

    # Calculate estimate for localresilience of O and Ō

    r_O = local_collective_resilience(abstractions_df, transitions_df, organization_indices, estimator)
    r_Ō = local_collective_resilience(abstractions_df, transitions_df, complement_indices, estimator)

    return [[r_O,1-r_Ō],[1-r_O,r_Ō]]

def estimate_global_collective_resilience(abstractions_df, transitions_df, organization_indices, estimator):
    '''
    Parameters:
        abstractions_df      (pandas.DataFrame),
        transitions_df       (pandas.DataFrame),
        organization_indices (list),
        estimator            (int), 1,2 or 3
    Returns:
        T (2D-array)
        Global collective resilience based on an
        estimate local collective resiliences by tomas or simons approach
    '''
    T = estimated_transitions_matrix(abstractions_df, transitions_df, organization_indices, estimator)
    T = T_inf(T)
    return matmul(T,[1,0])[0]

def estimate_global_collective_resilience(abstractions_df, transitions_df, organization_indices, type):
    '''
    Parameters:
        abstractions_df      (pandas.DataFrame),
        transitions_df       (pandas.DataFrame),
        organization_indices (list),
        type                 (str), 'tomas' or 'simon'
    Returns:
        T (2D-array)
        Global collective resilience based on an
        estimate local collective resiliences by tomas or simons approach
    '''
    T = estimated_transitions_matrix(abstractions_df, transitions_df, organization_indices, type)
    T = T_inf(T)
    return matmul(T,[1,0])[0]

####################################################################################################
####################### IDENTIFICATION OF INTERESTING SETS OF ORGANIZATIONS ########################
####################################################################################################

def find_index_of_oganization_with_minimal_resilience(organization_index, abstractions_df):
    '''
    Parameters:
        organization_index (list), a list of allowed organization indices
        abstractions_df    (pandas.DataFrame)
    Returns:
        Of all indices in organization_index one is returned where the global resilience is minimal
    '''
    minimal_resilience = min(abstractions_df.iloc[organization_index]['global_resilience'])
    return choice([i for i in organization_index if abstractions_df.iloc[i]['global_resilience']==minimal_resilience])

def algorithm(abstractions_df, transitions_df):
    '''
    Parameters:
        abstractions_df (pandas.DataFrame),
        transitions_df  (pandas.DataFrame)
    Returns:
        global_resilience_course (list)
        sequence_of_removals     (list)
    Starting from the full set of organizations this algorithm succesively
    organizations such that
    a) the removal of the organization leaves the remaining set of organizations connected
    b) the removal causes the minimal possible loss of global collective resilience 
    For each step the global collective resilience of the current subset and the index of
    the currently removed organization are documented in the lists that will be returned
    '''

    # Creating a directed Graph of the transitions
    G = transition_graph_from_dataframes(abstractions_df, transitions_df)
    
    # Genereating a networkx.DiGraph object for the markov-model
    organization_index = list(range(len(abstractions_df)))

    global_resilience_course = []
    sequence_of_removals     = ['-']

    while organization_index:

        # Evaluate global collective resiliece
        global_resilience_course.append(global_collective_resilience(abstractions_df, organization_index))

        # Check which organizations can be removed without disrupting the connectedness of the set
        try:
            removable_index = [i for i in organization_index if is_weakly_connected(G.subgraph([i_ for i_ in organization_index if i_!=i]))]
            if not removable_index:
                break
        except:
            break

        # Remove an organization with minimal global resilience
        organization_index_to_remove = find_index_of_oganization_with_minimal_resilience(organization_index, abstractions_df)
        organization_index.remove(organization_index_to_remove)
        sequence_of_removals.append(organization_index_to_remove)

    return global_resilience_course+[0], sequence_of_removals+[organization_index[0]]
        
def plot_for_algorithm(abstractions_df, global_resilience_course, sequence_of_removals, axes):

    # Plotting curve of collective resilience
    axes.plot(range(len(global_resilience_course)), global_resilience_course, lw=4, label='Collective global resilience')

    # Plotting additional information on removed organizations
    for key in ['number of species', 'complexity', 'local_resilience']:
        values            = [abstractions_df.iloc[i][key] for i in sequence_of_removals[1:]]
        normalized_values = [value/max(values) for value in values]
        axes.plot(list(range(1,len(normalized_values)+1)), normalized_values, label=f'{key} of removed organization', linestyle='dashed')

    # Labeling etc
    axes.tick_params(labelsize=20)
    axes.set_xlabel('Number of excluded organizations\nIndex of last removed organization', fontsize=20)
    axes.set_xticks(list(range(len(sequence_of_removals))), [f'{i}\n{r}' for i, r in enumerate(sequence_of_removals)])
    axes.legend(fontsize=15)
    axes.grid()

def plot_removal_sequence(sequence_of_removals, global_resilience_course, axes):
    axes.tick_params(labelsize=20)
    axes.set_xlabel('index of organization', fontsize=20)
    axes.set_ylabel('iteration', fontsize=20)
    organizations =  list(range(len(sequence_of_removals)-1))
    M = [[1 for organization_1 in organizations] for organization_1 in organizations]
    for i, organization in enumerate(sequence_of_removals[1:]):
        M[i] = [global_resilience_course[i] if j in organizations else 0 for j in range(len(sequence_of_removals)-1)]
        organizations.remove(organization)
    axes.imshow(M)