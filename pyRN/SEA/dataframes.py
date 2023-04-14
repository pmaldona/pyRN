from . import sos
from . import dataRetrivalService as drs
import itertools
import pandas
from . import markov

def initialize_abstractions_df(rndws):
    '''
    In:  rndws (list of random walks,
                where each random walk is a list of abstractions and
                where each abstraction is a list of 0 and 1)
    Out: pandas dataframe with only one column: the different abstractions 
    '''
    abstractions = list(itertools.chain(*rndws))                # Creates a joined list of all abstractions
    abstractions = [str(a) for a in abstractions]               # Converting abstractions to a string representation for elimination of duplicates
    abstractions = list(set(abstractions))                      # Filters all duplicates
    abstractions = [sos.from_string(a) for a in abstractions]   # Converting abstractions to binary array representation to count number of elements
    abstractions.sort(key=sos.n_elements)                       # Sorting abstractions by their number of species       
    abstractions = [str(a) for a in abstractions]               # Converting abstractions to a string representation for storing in the dataframe

    return pandas.DataFrame({'abstraction': abstractions})

def abstractions_list(abstractions_df):
    '''
    In:  abstractions_df
    Out: list of the binary array representations of the abstractions in abstractions_df
    '''
    abstractions = list(abstractions_df['abstraction'])
    abstractions = [sos.from_string(a) for a in abstractions]

    return abstractions

def abstractions_df_add_species_number(abstractions_df):
    '''
    In:  abstractions_df (pandas dataframe)
    Out: abstractions_df
    Adds a column with the number of species to the dataframe
    '''
    abstractions = abstractions_list(abstractions_df)
    n_species    = [sos.n_elements(a) for a in abstractions]

    return abstractions_df.assign(n_species=n_species)

def abstractions_df_add_initial_distribution(abstractions_df, rndws):
    '''
    In:  abstractions_df (pandas dataframe)
    Out: abstractions_df
    Adds a column with the probability of each abstraction to be the initial one of an random walk
    '''
    abstractions = abstractions_list(abstractions_df)
    init_distrib = [0 for i in range(len(abstractions))]
    for rndw in rndws:
        for i,abstraction in enumerate(abstractions):
            if rndw[0] == abstraction:
                init_distrib[i] += 1
    init_distrib = [x/len(rndws) for x in init_distrib]

    return abstractions_df.assign(initial_distribution=init_distrib)

def abstractions_df_add_complexities(abstractions_df, rndws, complexities):
    '''
    In:  abstractions_df (pandas dataframe)
         rndws           (list of random walks,
                          where each random walk is a list of abstractions and
                          where each abstraction is a list of 0 and 1)
         complexities    (2D-array with complexities[i][j] being the
                          complexity of the system with abstraction rndws)
    Out: abstractions_df
    Adds a column with the complexities of each abstraction
    '''
    abstractions = abstractions_list(abstractions_df)
    df_complexities = [0 for a in abstractions]
    # Search for the positions of the abstractions from the df in the random walk
    # and gets the corresponding complexity value
    for df_index, df_abstraction in enumerate(abstractions):
        finded=False
        for rndw_index,rndw in enumerate(rndws):
            for a_index, a in enumerate(rndw):
                if a  == df_abstraction:
                    df_complexities[df_index] = complexities[rndw_index][a_index]
                    finded = True
                    break
            if finded:
                break
        
    return abstractions_df.assign(complexity=df_complexities)

def initialize_transitions_df(rndws):
    '''
    In:  rndws (list of random walks,
                where each random walk is a list of abstractions and
                where each abstraction is a list of 0 and 1)
    Out: pandas dataframe with only one column: the different abstractions 
    '''
    transitions    = []
    counts         = []
    for rndw in rndws:
        for i in range(len(rndw)-1):
            initial_state    = str(rndw[i])
            convergent_state = str(rndw[i+1])
            if (initial_state,) not in transitions:
                transitions.append((initial_state,convergent_state))
                counts.append(1)
            else:
                counts[transitions.index((initial_state,convergent_state))] +=1
    
    transitions_df = pandas.DataFrame({'initial_state': [t[0] for t in transitions], 'convergent_state': [t[1] for t in transitions], 'counts': counts})
    return transitions_df

def transitions_df_add_set_changes(transitions_df):
    '''
    In:  transitions_df (pandas dataframe)
    Out: transitions_df
    Adds a column with the number of species to the dataframe
    '''
    nrows = transitions_df.shape[0]
    transitions_df.assign(change_local=[0 for i in range(nrows)], change_global=[0 for i in range(nrows)])
    for i in range(nrows):
        a1 = sos.from_string(transitions_df.loc[i, 'initial_state']) # Abstraction before transition
        a2 = sos.from_string(transitions_df.loc[i, 'convergent_state']) # Abstraction after transition
        transitions_df.loc[i, 'change_local']  = sos.change_local(a1,a2)
        transitions_df.loc[i, 'change_global'] = sos.normalized_hamming_distance(a1,a2)

    return transitions_df

def transitions_df_add_complexity_changes(transitions_df, abstractions_df):
    '''
    In:  transitions_df (pandas dataframe)
    Out: transitions_df
    Adds a column with the number of species to the dataframe
    '''
    nrows = transitions_df.shape[0]
    transitions_df.assign(change_complexity=[0 for i in range(nrows)])
    for i in range(nrows):
        a1 = transitions_df.loc[i, 'initial_state']    # Abstraction before transition
        a2 = transitions_df.loc[i, 'convergent_state'] # Abstraction after transition
        c1 = abstractions_df.loc[abstractions_df['abstraction'] == a1, 'complexity']
        c2 = abstractions_df.loc[abstractions_df['abstraction'] == a2, 'complexity']
        transitions_df.loc[i, 'change_complexity']  = c1.values[0]-c2.values[0]

    return transitions_df

def transitions_df_add_abstraction_indexes(transitions_df, abstractions_df):
    '''
    In:  transitions_df  (pandas dataframe),
         abstractions_df (pandas dataframe)
    Out: transitions_df
    Adds a column with the indices of the abstractions in the abstractions_df to the transitons_df
    '''
    a1 = []
    a2 = []
    for i in range(transitions_df.shape[0]):
        a1.append(list(abstractions_df.index[abstractions_df['abstraction']==transitions_df.loc[i, 'initial_state']])[0])
        a2.append(list(abstractions_df.index[abstractions_df['abstraction']==transitions_df.loc[i, 'convergent_state']])[0])
    transitions_df.insert(0, 'initial_state', a2)
    transitions_df.insert(0, 'convergent_state', a1)

    return transitions_df

def transitions_df_add_transitions_probabilities(df):
    '''
    In:  transitions_df (pandas dataframe), with at least three columns: a_1, a_2 and count
    Out: transitions_df (pandas dataframe)
    Adds a column to the transitions_df with the transition probabilities in the Markov-model
    '''
    df.assign(probability=[0 for i in range(df.shape[0])])
    starts = df['initial_state'].unique()
    for start in starts:
        transitions_from_start = df.loc[df['a_1']==start]
        n = transitions_from_start['counts'].sum()
        ends = transitions_from_start['a_2']
        for end in ends:
            df.loc[(df['initial_state']==start) & (df['convergent_state']==end),'probability'] = df.loc[(df['initial_state']==start) & (df['convergent_state']==end),'counts']/n
    return df

def dataframes(path, abstraction_type):
    '''
    In: path             (string), path to .json-file with random walks,
        abstraction_type (string), identifier of the type of abstraction we are interested in
    '''
    abst = drs.get_RNDWs(path, abstraction_type)     # Reads abstractions from .json-file
    cps = drs.get_CPs(path, 'c'+abstraction_type)

    transitions_df  = initialize_transitions_df(abst)
    abstractions_df = initialize_abstractions_df(abst)

    abstractions_df = abstractions_df_add_initial_distribution(abstractions_df, abst)
    abstractions_df = abstractions_df_add_species_number(abstractions_df)
    abstractions_df = abstractions_df_add_complexities(abstractions_df, abst, cps)

    transitions_df = transitions_df_add_abstraction_indexes(transitions_df, abstractions_df)
    transitions_df = transitions_df_add_transitions_probabilities(transitions_df)
    transitions_df = transitions_df_add_set_changes(transitions_df)
    transitions_df_add_complexity_changes(transitions_df, abstractions_df)

    abstractions_df = markov.add_markov_properties_to_dataframe(abstractions_df, transitions_df)

    return abstractions_df, transitions_df

def dataframesFromLists(abstraction_list,complexity_list):
    '''
    Parameters
    ----------
    abstraction_list : List of lists
        Lists of all abstraction form a random walk.
    complexity_list : TYPE
        Lists of all complexities form a random walk correlative to the abstraction_list.

    Returns
    -------
    abstractions_df : Pandas DataFrame
        Abstraction dataframe whit resilient data.
    transitions_df : Pandas DataFrame
        Transition dataframe whit resilient data.

    '''
    abst = abstraction_list     # Reads abstractions from .json-file
    cps = complexity_list

    transitions_df  = initialize_transitions_df(abst)
    abstractions_df = initialize_abstractions_df(abst)

    abstractions_df = abstractions_df_add_initial_distribution(abstractions_df, abst)
    abstractions_df = abstractions_df_add_species_number(abstractions_df)
    abstractions_df = abstractions_df_add_complexities(abstractions_df, abst, cps)

    transitions_df = transitions_df_add_abstraction_indexes(transitions_df, abstractions_df)
    transitions_df = transitions_df_add_transitions_probabilities(transitions_df)
    transitions_df = transitions_df_add_set_changes(transitions_df)
    transitions_df_add_complexity_changes(transitions_df, abstractions_df)

    abstractions_df = add_markov_properties_to_dataframe(abstractions_df, transitions_df)

    return abstractions_df, transitions_df


    def dataframesFromAbst(abstrac_list):
        '''
        Parameters
        ----------
        abstlist : list of list of 0,1 vectors
            List of all random walks abstractions.

        Returns
        -------
        Dataframe of abstractions and transition whit resilience infromation.
        '''
        
        transitions_df  = initialize_transitions_df(abst)
        abstractions_df = initialize_abstractions_df(abst)

        abstractions_df = abstractions_df_add_initial_distribution(abstractions_df, abst)
        abstractions_df = abstractions_df_add_species_number(abstractions_df)
        abstractions_df = abstractions_df_add_complexities(abstractions_df, abst, cps)

        transitions_df = transitions_df_add_abstraction_indexes(transitions_df, abstractions_df)
        transitions_df = transitions_df_add_transitions_probabilities(transitions_df)
        transitions_df = transitions_df_add_set_changes(transitions_df)
        transitions_df_add_complexity_changes(transitions_df, abstractions_df)

        abstractions_df = markov.add_markov_properties_to_dataframe(abstractions_df, transitions_df)
        
        

# a1 = [[0,0,1,0],[1,1,1,0]]
# a2 = [[1,0,1,0],[1,0,1,0]]
# a3 = [[0,0,1,1],[1,1,0,0]]
# a=[a1,a2,a3]


# import json
# path="../Research/RGRW/Nsp10_Nsteps10_rep2_net1.json"
# with open(path, 'r') as f:
#         rndws = json.loads(f.read())

# abst = drs.get_RNDWs(path, 'a')     # Reads abstractions from .json-file
# cps = drs.get_CPs(path, 'ca')
# transitions_df  = initialize_transitions_df(abst)
# abstractions_df = initialize_abstractions_df(abst)

# abstractions_df = abstractions_df_add_initial_distribution(abstractions_df, abst)
# abstractions_df = abstractions_df_add_species_number(abstractions_df)
# abstractions_df = abstractions_df_add_complexities(abstractions_df, abst, cps)

# transitions_df = transitions_df_add_abstraction_indexes(transitions_df, abstractions_df)
# transitions_df = transitions_df_add_transitions_probabilities(transitions_df)
# transitions_df = transitions_df_add_set_changes(transitions_df)
# transitions_df_add_complexity_changes(transitions_df, abstractions_df)

# abstractions_df = add_markov_properties_to_abstractions_df(abstractions_df, transitions_df)
# # data=dataframes(path,'a')[0]

