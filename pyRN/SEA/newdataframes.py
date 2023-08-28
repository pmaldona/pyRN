'''
This scipt provides methods to calculate the dataframes
with information about abstractions and transitions 
including a Markov-model for the system evolution.

The Markov-model can be calculated in various ways based
on different parameterizations of the perturbations.

TABLE OF CONTENTS:                                                                                          | Methods
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Dataframe initialialization             |                   |                                               | initialize_abstractions_df()
                                        |                   |                                               | initialize_transitions_df()
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Helper functions                        |                   |                                               | find_all_perturbations()
                                        |                   |                                               | probability_perturbation_leads_to_convergent_state()
                                        |                   |                                               | determine_perturbation_size()
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Calculating transition probabilities    | 1. for conn=False | Parameterization by size of the perturbation  | additive_perturbations()
                                        |                   |                                               | eliminative_perturbations(initial_state)
                                        |                   |                                               | number_of_possible_perturbations_by_size_matrix_using_combinatorics()
                                        |                   |                                               | number_of_possible_perturbations_by_size_matrix_using_dataframe()
                                        |                   |                                               | number_of_observed_perturbations_by_size_matrix()
                                        |                   |                                               | fix_perturbation_probability_vector()
                                        |                   |                                               | additive_probabilities()
                                        |                   |                                               | eliminative_probabilities(initial_state, alpha)
                                        |                   |                                               | perturbation_size_probabilities_matrix()
                                        |                   |                                               | transition_probability_1_1()
                                        |                   |                                               | add_probabilities_to_transitions_df_1_1()
                                        |                   |-------------------------------------------------------------------------------------------------------------------------
                                        |                   | Parameterization by species specific rates    | perturbation_probability_1_2()
                                        |                   |                                               | fix_perturbation_probabilities()  (in Progress)
                                        |                   |                                               | add_perturbation_probabilities_to_SimpleTransSpDf()
                                        |                   |                                               | transition_probability_1_2()
                                        |---------------------------------------------------------------------------------------------------------------------------------------------
                                        | 2. for conn=True  | Uniform distribution                          | transition_probability_2_1()
                                        |                   |                                               | add_probabilities_to_transitions_df_2_1()
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Fixing numerical errors from the        |                   |                                               | fix_transition_probabilities()
transition probability calculations     |                   |                                               | fix_transition_probabilities_for_all_initial_states())
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Adding additional information           |                   |                                               | add_number_of_species()
                                        |                   |                                               | add_complexities()
                                        |                   |                                               |
'''

from pyRN.SEA import sos

import math
import numpy
import pandas
from   bitarray      import bitarray
from   bitarray      import frozenbitarray
from   scipy.special import binom
import pandas as pd
from bitarray import frozenbitarray as fbt
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
'''
Dataframe initialization
'''
###############################################################################################################################################################################
def TransDfFromRw(RwDict,Rw_type="simple",init_label='i',pert_label='pr',pert_state_label='p',conver_label='c',w=None):
    
    if w is None:
        w=RwDict[Rw_type].keys()
        
    TransDf=pd.DataFrame(pd.DataFrame(columns=["initial state",'perturbation',"final perturbed state","convergent state"]))
        
    for i in w:
        init_states=list(map(lambda x: fbt(RwDict[Rw_type][i][init_label][x].tolist()),
                             RwDict[Rw_type][i][init_label].columns))
        pertrurbations=list(map(lambda x: fbt(RwDict[Rw_type][i][pert_label][x].tolist()),
                             RwDict[Rw_type][i][pert_label].columns))
        pert_states=list(map(lambda x: fbt(RwDict[Rw_type][i][pert_state_label][x].tolist()),
                             RwDict[Rw_type][i][pert_state_label].columns))
        conv_states=list(map(lambda x: fbt(RwDict[Rw_type][i][conver_label][x].tolist()),
                             RwDict[Rw_type][i][conver_label].columns))
        tmp_df=pd.DataFrame(list(map(list, zip(*[init_states,pertrurbations,pert_states,conv_states]))),
                            columns=["initial state",'perturbation',"final perturbed state","convergent state"])
        TransDf=pd.concat([TransDf,tmp_df],axis=0,ignore_index=True)
        
    return TransDf

def initialize_abstractions_df(SimpleTransSpDf):
    '''
    Returns a Pandas dataframe with the abstractions
    in the SimpleTransSpDf as frozenbitarrays.

        Parameters:
            SimpleTransSpDf (pandas.DataFrame)

        Returns:
            abstractions_df (pandas.DataFrame)
    '''
    abstractions    = list(SimpleTransSpDf['initial state'])       # Get a list of all initial abstractions from the SimpleTransSpDf dataframe
    abstractions    = list(set(abstractions))                         # Getting a list of unique abstractions
    abstractions    = [list(a) for a in abstractions]                 # Convert to list
    abstractions.sort(key=sos.n_elements)                             # Sort by size
    abstractions    = [frozenbitarray(a) for a in abstractions]       # Reconverting
    abstractions_df = pandas.DataFrame({'abstraction': abstractions}) # Create the dataframe
    return abstractions_df

def initialize_transitions_df(abstractions_df):
    '''
    Returns a Pandas dataframe with the transitions
    in the SimpleTransSpDf with abstractions as
    frozenbitarrays.

        Parameters:
            abstractions_df, (pandas.DataFrame)

        Returns:
            transitions_df, (pandas.DataFrame)
    '''
    transitions = [(a1,a2) for a1 in abstractions_df['abstraction'] for a2 in abstractions_df['abstraction']] # Get all pairs of abstractions from the abstractions_df dataframe
    a_1 = [t[0] for t in transitions]                                                                         # Get a list of all initial abstractions
    a_2 = [t[1] for t in transitions]                                                                         # Get a list of all convergent abstractions
    transitions_df = pandas.DataFrame({'initial_state': a_1, 'convergent_state': a_2})                        # Create the dataframe
    return transitions_df

###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
'''
Helper functions that are used in multiple other functions
'''
###############################################################################################################################################################################

def find_all_perturbations(SimpleTransSpDf, initial_state, convergent_state):
    '''
    Identifies all perturbations that cause a transition
    from the initial_state to the convergent_state and
    returns them as a list of bitarrays

        Parameters:
            SimpleTransSpDf  (pandas.DataFrame)
            initial_state    (some binary representation)
            convergent_state (some binary representation)'
            Allowed binary representations:
                '01',
                bitarray('01'),
                frozenbitarray('01'),
                [0,1],
                [False, True]

        Returns:
            perturbations (list)
    '''
    perturbations = SimpleTransSpDf.loc[(SimpleTransSpDf['initial state']==bitarray(initial_state)) &
                                           (SimpleTransSpDf['convergent state']==bitarray(convergent_state))]['perturbation']
    perturbations = list(set(perturbations)) # Filter for unique perturbations
    return perturbations

def probability_perturbation_leads_to_convergent_state(SimpleTransSpDf, initial_state, convergent_state, perturbation):
    '''
    Determines the frequency at which a perturbation of an
    initial_state leads to a convergent_state.
    (If a perturbation causes the system to transition to
    a closed set that is not an organization itself the 
    system can transition to various organizations below,
    with the same probability for each of them)

        Parameters:
            SimpleTransSpDf  (pandas.DataFrame)
            initial_state    (some binary representation)
            convergent_state (some binary representation)'
            Allowed binary representations:
                '01',
                bitarray('01'),
                frozenbitarray('01'),
                [0,1],
                [False, True]

        Returns:
            perturbations (list)
    '''
    events = SimpleTransSpDf.loc[(SimpleTransSpDf['initial state']==bitarray(initial_state)) &
                                    (SimpleTransSpDf['perturbation']==bitarray(perturbation))]
    n = events.shape[0]
    events = events.loc[events['convergent state']==bitarray(convergent_state)]
    k = events.shape[0]
    return k/n

def determine_perturbation_size(state, perturbation):
    '''
    In: state,
        perturbation
    Out: perturbation_size (tuple), (#added_species, #removed_species)
    '''
    state        = list(state)
    perturbation = list(perturbation)
    added        = sum([1 for i in range(len(state)) if state[i]==0 and perturbation[i]==1])
    removed      = sum([1 for i in range(len(state)) if state[i]==1 and perturbation[i]==1])
    return (added, removed)

###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
'''
1  Creating a markov model for both connected and not connected additive perturbations (for conn=False)
'''
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
'''
1.1 Parameterization by size of the perturbation

The idea is to operate with 3 matrices:<br>

- number_of_possible_perturbations_per_size_matrix
- number_of_observed_perturbations_per_size_matrix
- perturbation_size_probability_matrix

number_of_possible_perturbations_per_size_matrix[i][j]:
the number of possible perturbations that can add i new species and remove j existing species from the system<br>

number_of_observed_perturbations_per_size_matrix[i][j]:
the number of observed perturbations that cause a specific transition by adding i new species and removing j existing species from the system<br>

perturbation_size_probability_matrix[i][j]:
the probability that a perturbation orrurs which adds i new species and removes j existing species from the system

All of the matrices have the same size (ns+1)^2. Their values of the depend on the initial_state before the transition.
'''
###############################################################################################################################################################################

###############################################################################################################################################################################
'''
1.1 number_of_possible_perturbations_per_size_matrix
'''
###############################################################################################################################################################################

def additive_perturbations(initial_state):
    '''
        Parameters:
            initial_state         (some binary representation)

        Returns:
            number of ways in which a maximum of max_perturbation_size species
            can be added to the initial_state
    '''
    n = len(initial_state)                  # Overall number of species
    k = sos.n_elements(list(initial_state)) # Number of species present
    return [binom(n-k,l) for l in range(len(initial_state)+1)]

def eliminative_perturbations(initial_state):
    '''
        Parameters:
            initial_state         (some binary representation)

        Returns:
            number of ways in which a maximum of max_perturbation_size species
            can be removed from the initial_state
    '''
    k = sos.n_elements(list(initial_state)) # Number of species present
    return [binom(k,l) for l in range(len(initial_state)+1)]

def number_of_possible_perturbations_by_size_matrix_using_combinatorics(initial_state):
    '''
    Returns a matrix number_of_possible_perturbations_by_size_matrix
    with number_of_possible_perturbations_by_size_matrix[i][j]
    the number of possible perturbations that can add i new species
    and remove j existing species from the system. No more than
    max_perturbation_size species can be added or removed from the
    initial_state.
    The calculation is done using combinatorics.

        Parameters:
            initial_state         (some binary representation)
            max_perturbation_size (int), maximum number of species to be added

        Returns:
            number_of_possible_perturbations_by_size_matrix
    '''
    a = additive_perturbations(initial_state)
    e = eliminative_perturbations(initial_state)
    return [[int(a[i]*e[j]) for j in range(len(initial_state)+1)] for i in range(len(initial_state)+1)]
    
def number_of_possible_perturbations_by_size_matrix_using_dataframe(SimpleTransSpDf, initial_state):
    '''
    Returns a matrix number_of_possible_perturbations_by_size_matrix
    with number_of_possible_perturbations_by_size_matrix[i][j]
    the number of possible perturbations that can add i new species
    and remove j existing species from the system. No more than
    max_perturbation_size species can be added or removed from the
    initial_state.
    The calculation is done by counting the perturbations in the 
    SimpleTransSpDf.

        Parameters:
            SimpleTransSpDf       (pandas.DataFrame)
            initial_state         (some binary representation)
            max_perturbation_size (int), maximum number of species to be added

        Returns:
            number_of_possible_perturbations_by_size_matrix
    '''
    t = SimpleTransSpDf.loc[SimpleTransSpDf['initial state']==bitarray(initial_state)]
    perturbations = list(set(t['perturbation'])) # Filter unique perturbations
    possible_transitions = [[0 for j in range(len(initial_state)+1)] for i in range(len(initial_state)+1)]
    for perturbation in perturbations:
        perturbation_size = determine_perturbation_size(initial_state, perturbation)
        possible_transitions[perturbation_size[0]][perturbation_size[1]] += 1
    return possible_transitions

###############################################################################################################################################################################
'''
1.2 number_of_observed_perturbations_per_size_matrix
'''
###############################################################################################################################################################################

def number_of_observed_perturbations_by_size_matrix(SimpleTransSpDf, initial_state, convergent_state, max_perturbation_size):
    '''
    Returns a matrix number_of_observed_perturbations_by_size_matrix
    with number_of_observed_perturbations_by_size_matrix[i][j]
    the number of observed perturbations that add i new species
    and remove j existing species from the system and lead to a
    specific transition. No more than max_perturbation_size species
    can be added or removed from the initial_state.
    The calculation is done by counting the perturbations in the 
    SimpleTransSpDf.

        Parameters:
            SimpleTransSpDf (pandas.DataFrame)
            initial_state
            convergent_state
            max_perturbation_size (int)

        Returns:
            number_of_observed_perturbations_by_size_matrix
    '''
    perturbations = find_all_perturbations(SimpleTransSpDf, initial_state, convergent_state)
    observed_perturbations = [[0 for j in range(max_perturbation_size+1)] for i in range(max_perturbation_size+1)]
    for perturbation in perturbations:
        perturbation_size = determine_perturbation_size(initial_state, perturbation)
        observed_perturbations[perturbation_size[0]][perturbation_size[1]] += probability_perturbation_leads_to_convergent_state(SimpleTransSpDf,
                                                                                                                                 initial_state,
                                                                                                                                 convergent_state,
                                                                                                                                 perturbation)
    return observed_perturbations

###############################################################################################################################################################################
'''
1.3 perturbation_size_probabilities_matrix

Each species has a probability alpha to be affected by perturbation.
This means that there is a probability of alpha that a species gets 
removed if it is present or added if it is not present.
'''
###############################################################################################################################################################################

def fix_perturbation_probability_vector(probability_vector):
    '''
    If the sum of probabilities of the input probability_vector is not 1
    all of the values that are not 0 are increased equally so probabilities
    add up to 1.

        Parameters:
            probability_vector (list)

        Returns:
            probability_vector (list)
    '''
    error = 1-numpy.sum(probability_vector)
    compensation = error / len([p for p in probability_vector if p>0])
    probability_vector = [p if p==0 else p+compensation for p in probability_vector]
    return probability_vector

def additive_probabilities(initial_state, alpha):
    n = len(initial_state)
    state       = list(initial_state)
    full_state  = [1 for i in range(n)]
    if state == full_state:
        return [1] + [0 for i in range(n)]
    non_present = sos.difference(full_state, state)
    p_additive  = [0 for i in range(n+1)]
    for i in range(sos.n_elements(non_present)+1):
        p_additive[i] = binom(n,i)*(alpha**i)*((1-alpha)**(n-i))
    return fix_perturbation_probability_vector(p_additive)

def eliminative_probabilities(state, alpha):
    state       = list(state)
    empty_state = [0 for i in range(len(state))]
    if state == empty_state:
        return [1] + empty_state
    p_eliminative  = [0 for i in range(len(state)+1)]
    n = sos.n_elements(state)
    for i in range(sos.n_elements(state)+1):
        p_eliminative[i] = binom(n,i)*(alpha**i)*((1-alpha)**(n-i))
    return fix_perturbation_probability_vector(p_eliminative)

def perturbation_size_probabilities_matrix(initial_state, alpha):
    '''
    Returns the perturbation_size_probabilities_matrix
    perturbation_size_probabilities_matrix[i][j] is the
    probability that a perturbation occurs that adds i
    new species  and removes j species.  Every species
    has the same probability to be affected in a
    perturbation event.

        Parameters:
            initial_state
            alpha (float), the probability for each  species to be affected by a perturbation

        Returns:
            perturbation_size_probabilities_matrix
    '''
    n = len(initial_state)+1
    a = additive_probabilities(initial_state, alpha)
    e = eliminative_probabilities(initial_state, alpha)
    M = [[a[i]*e[j] for j in range(n)] for i in range(n)]
    return M

###############################################################################################################################################################################
'''
1.4 Calculating transition probabilities and add them to the dataframe

With the three matrices we can now calculate transition probabilities as
Σ_{i,j}(probabilities[i][j]*observed[i][j]/possible[i][j])
'''
###############################################################################################################################################################################

def transition_probability_1_1(SimpleTransSpDf, initial_state, convergent_state, alpha):
    '''
    Returns the probability of a transition given the system
    is in the initial_state. Probabilities are calculated as
    Σ_{i,j}(probabilities[i][j]*observed[i][j]/possible[i][j])

        Parameters:
            SimpleTransSpDf (pandas.DataFrame)
            initial_state
            convergent_state
            alpha (float), the probability for each  species to be affected by a perturbation

        Returns:
            transition probability
    '''
    max_perturbation_size = len(initial_state)
    # Perturbations
    possible      = number_of_possible_perturbations_by_size_matrix_using_dataframe(SimpleTransSpDf,  initial_state)
    observed      = number_of_observed_perturbations_by_size_matrix(SimpleTransSpDf, initial_state, convergent_state, max_perturbation_size)
    probabilities = perturbation_size_probabilities_matrix(initial_state, alpha)

    M = [[0 for j in range(max_perturbation_size+1)] for i in range(max_perturbation_size+1)]
    for i in range(len(M)):
        for j in range(len(M)):
            if possible[i][j] > 0:
                M[i][j] = probabilities[i][j]*observed[i][j]/possible[i][j]
    return numpy.sum(M)

def add_probabilities_to_transitions_df_1_1(SimpleTransSpDf, transitions_df, alpha):
    '''
    Adds a column with transition probabilities
    to an existing transitions_df

        Parameters:
            SimpleTransSpDf (pandas.DataFrame)
            transitions_df  (pandas.DataFrame)
            alpha           (float), the probability for each  species to be affected by a perturbation

        Returns:
    '''
    probabilities = []
    for row in transitions_df.iterrows():
        transition = (bitarray(row[1][0]), bitarray(row[1][1]))
        probabilities.append(transition_probability_1_1(SimpleTransSpDf, transition[0], transition[1], alpha))
    transitions_df['probability']=probabilities

###############################################################################################################################################################################
###############################################################################################################################################################################

def perturbation_probability_1_2(initial_state, perturbation, appearing_rates, disappearing_rates):
    '''
    Returns the probability of a perturbation for an initial_state

    Parameter:
        initial_state,
        perturbation,
        appearing_rates, 
        disappearing_rates

    Returns:
        probability
    '''
    probability = 1
    for i in range(len(initial_state)):
        if initial_state[i]:
            if perturbation[i]:
                probability = probability*disappearing_rates[i]
            else:
                probability = probability*(1-disappearing_rates[i])
        else:
            if perturbation[i]:
                probability = probability*appearing_rates[i]
            else:
                probability = probability*(1-appearing_rates[i])
    return probability

def fix_perturbation_probabilities(SimpleTransSpDf):
    '''
    If (Σ P(perturbation|initial_state)) < 1 the probabilities will be increased equally
    '''
    initial_states = SimpleTransSpDf['initial state'].unique().tolist()
    for initial_state in initial_states:
        df = SimpleTransSpDf.loc[(SimpleTransSpDf['initial state']==initial_state) & (SimpleTransSpDf['perturbation_probability']>0)]
        error = 1 - (sum(df['perturbation_probability'].tolist()))/df.shape[0]
        df['perturbation_probability'] += error
         
def add_perturbation_probabilities_to_SimpleTransSpDf(SimpleTransSpDf, appearing_rates, disappearing_rates):
    '''
    Adds a column 'perturbation_probability' to the SimpleTransSpDf pandas.DataFrame

    Parameter:
        SimpleTransSpDf    (pandas.DataFrame)
        appearing_rates    (list), list of floats with values v, 0≤v≤1
        disappearing_rates (list), list of floats with values v, 0≤v≤1

    Returns:
    '''
    perturbation_probabilities = [0 for i in range(SimpleTransSpDf.shape[0])]
    for i, row in SimpleTransSpDf.iterrows():
        perturbation_probabilities[i] = perturbation_probability_1_2(row[0], row[1], appearing_rates, disappearing_rates)
    SimpleTransSpDf['perturbation_probability'] = perturbation_probabilities

def transition_probability_1_2(SimpleTransSpDf, initial_state, convergent_state, appearing_rates, disappearing_rates):
    perturbations= SimpleTransSpDf.loc[(SimpleTransSpDf['initial state']==initial_state) & (SimpleTransSpDf['convergent state']==convergent_state)]['perturbation'].unique()
    transition_probability = 0
    for perturbation in perturbations.tolist():
        a = probability_perturbation_leads_to_convergent_state(SimpleTransSpDf, initial_state, convergent_state, perturbation)
        b = float(SimpleTransSpDf.loc[(SimpleTransSpDf['initial state']==initial_state) & (SimpleTransSpDf['perturbation']==perturbation)]['perturbation_probability'])
        if not math.isnan(a*b):
            transition_probability += a*b
    return transition_probability

###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################

def transition_probability_2_1(SimpleTransSpDf, initial_state, convergent_state):
    '''
    Calculate transition probabilities
    '''
    perturbations_of_initial_state = (SimpleTransSpDf.loc[SimpleTransSpDf['initial state']==initial_state]['perturbation'])
    perturbations_of_initial_state = list(set(perturbations_of_initial_state.tolist()))
    o = 0
    n = len(perturbations_of_initial_state)
    transition_probability = 0
    for perturbation in perturbations_of_initial_state:
        df = SimpleTransSpDf.loc[(SimpleTransSpDf['initial state']==initial_state) &
                                    (SimpleTransSpDf['perturbation']==perturbation) &
                                    (SimpleTransSpDf['convergent state']==convergent_state)]
        if not df.empty:
            p = probability_perturbation_leads_to_convergent_state(SimpleTransSpDf, initial_state, convergent_state, perturbation)
            
            transition_probability += p/n
            o += p
    values = {'possible': n,
              'observed': o,
              'probability': transition_probability}                                  
    return values

def add_probabilities_to_transitions_df_2_1(SimpleTransSpDf, transitions_df):
    '''
    Adds a column with transition probabilities
    to an existing transitions_df

        Parameters:
            SimpleTransSpDf (pandas.DataFrame)
            transitions_df  (pandas.DataFrame)

        Returns:
    '''
    possible    = []
    observed    = []
    probability = []
    for row in transitions_df.iterrows():
        transition = (bitarray(row[1][0]), bitarray(row[1][1]))
        values = transition_probability_2_1(SimpleTransSpDf, transition[0], transition[1])
        possible.append(values['possible'])
        observed.append(values['observed'])
        probability.append(values['probability'])
    transitions_df['perturbations_of_initial_state'] = possible
    transitions_df['perturbations_causing_transition'] = observed
    transitions_df['probability'] = probability

###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
'''
Fixing numerical errors from the transition probability calculations
'''
###############################################################################################################################################################################

def fix_transition_probabilities(transitions_df, initial_state):
    '''
    If (Σ P(perturbation|initial_state)) < 1 the highest probability is increased
    to prevent future error propagation / error amplification

    Parameters:
        transitions_df,
        initial_state
    '''
    
    # Determine what and where to update  probabilities
    transitions_from_initial_state         = transitions_df.query('initial_state==@initial_state')
    probabilities                          = transitions_from_initial_state['probability'].tolist()
    if any(numpy.array(probabilities)!=0): 
        error                                  = 1-sum(probabilities)
        maximimum_probability                  = max(probabilities)
        transition_with_max_probability        = transitions_from_initial_state.query('probability==@maximimum_probability')
        convergent_states_with_max_probability = transition_with_max_probability['convergent_state'].tolist()
        new_probability                        = maximimum_probability + (error/len(convergent_states_with_max_probability))

        # Update
        mask = (transitions_df[['initial_state', 'probability']] == [initial_state, maximimum_probability]).all(axis=1)
        for convergent_state in convergent_states_with_max_probability:
            transitions_df.loc[mask, 'probability'] = new_probability
    else:
        transitions_df.loc[(transitions_df['initial_state']==initial_state) &\
                       (transitions_df['convergent_state']==initial_state),'probability']=1
            
def fix_transition_probabilities_for_all_initial_states(transitions_df):
    '''
    For all initial states:
    If (Σ P(perturbation|initial_state)) < 1 the highest probability is increased
    to prevent future error propagation / error amplification

    Parameters:
        transitions_df
    '''
    initial_states = transitions_df['initial_state'].unique()
    for initial_state in initial_states:
        fix_transition_probabilities(transitions_df, initial_state)

###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
'''
Adding additional information to the dataframes
'''
###############################################################################################################################################################################

def add_number_of_species(abstractions_df):
    '''
    Adds a new column with the number of species to 
    an existing abstractions_df

        Parameters:
            abstractions_df (pandas.DataFrame)
    '''
    abstractions = abstractions_df['abstraction'].tolist()
    n_species    = [sos.n_elements(abstraction) for abstraction in abstractions]
    abstractions_df['number of species'] = n_species

def add_complexities(RN, abstractions_df, elem_type="generators"):
    '''
    Adds a new column with the complexities to 
    an existing abstractions_df

        Parameters:
            RN              (pyRN object), a reaction network,
            abstractions_df (pandas.DataFrame),
            elem_type       (string) optional,
                            default: generators, options: generators, basics
    '''
    abstractions = [bitarray(abstraction) for abstraction in abstractions_df['abstraction']]
    complexities = [RN.getComplexityFloat(abstraction,abst_type="species",elem_type="generators") for abstraction in abstractions]
    abstractions_df['complexity'] = complexities
    return abstractions_df

def add_size_difference(transitions_df):
    '''
    Adds a new column with the size (number of species) difference to 
    an existing transitions_df

        Parameters:
            transitions_df (pandas.DataFrame)
    '''
    org_size = lambda organization: sos.n_elements(list(organization))
    size_differences = [org_size(transitions_df.iloc[i]['convergent_state'])-org_size(transitions_df.iloc[i]['initial_state']) for i in range(transitions_df.shape[0])]
    transitions_df['size_difference'] = size_differences
