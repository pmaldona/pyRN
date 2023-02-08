import numpy
import copy

def T_inf(T, max_steps=500):
    '''
    In:  T     (2D-array), the transition matrix for a finite number of steps
    Out: T_inf (2D-array), the transition matrix for infinitive steps
    Calculates the transition matrix for infinitive steps by multiplying 
    transition matrix for a finite number of steps with itself until
    it doesn't change significantly anymore
    (Waiting for it to not chang at all could lead to bad results
    due to numerical errors)
    '''
    for step in range(max_steps):
        T_inf = [[T[i][j] for j in range(len(T))] for i in range(len(T))]
        T     = numpy.matmul(T,T)
        if (numpy.all(numpy.round(T,10)==numpy.round(T_inf,10))):
            return T_inf

def pn(p0,T,n):
    Tn = T
    for i in range(n):
        Tn = numpy.matmul(Tn,T)
        pn = numpy.matmul(Tn,p0)
    return pn

def reachability(p0,gs,T,t,max=100):
    '''
    In:  p0 initial distribution
         gs goal state
         T  transition matrix
         t  probability threshold
    Out: number of steps needed to have a probability of t to get to gs starting from p0 
    '''

    T_ = copy.copy(T)

    # Changing the transition matrix so gs cannot be left again
    
    T_[:,gs] = [0 for i in range(len(T_))]
    T_[gs][gs] = 1

    n      = 0     # number of steps
    pn     = p0    # probability distribution after n steps 
    Tn     = T_    # transition matrix for n steps

    while pn[gs]<t and n<max+1:
        n = n+1
        pn = numpy.matmul(Tn,pn)
        Tn = numpy.matmul(Tn,T_)
    return(n)

def add_transition_probbilities_to_dataframe(df):
    '''
    In:  df (pandas dataframe), with at least three columns: initial_state, convergent_state and count
    Out: df (pandas dataframe)
    Adds a column to thdf with the transition probabilities in the markov-model
    '''
    df = df.assign(probability=[0 for i in range(df.shape[0])])
    starts = df['initial_state'].unique()
    for start in starts:
        transitions_from_start = df.loc[df['initial_state']==start]
        n = transitions_from_start['counts'].sum()
        ends = transitions_from_start['convergent_state']
        for end in ends:
            df.loc[(df['initial_state']==start) & (df['convergent_state']==end),'probability'] = df.loc[(df['initial_state']==start) & (df['convergent_state']==end),'counts']/n
    return df

def transition_matrix_from_dataframes(abstractions_df, transitions_df):
    '''
    In:  df (pandas dataframe), with at least three columns: initial_state, a_2 and either count or probability
    Out: Transition mastrix of Markov-model
    '''
    if ('probability' not in transitions_df.columns):
        df = add_transition_probbilities_to_dataframe(df)
    
    n = abstractions_df.shape[0]
    t_matrix = [ [0 for j in range(n)] for i in range(n)]
    for i in range(n):
        start = abstractions_df.loc[i,'abstraction']
        transitions_from_start = transitions_df.loc[transitions_df['initial_state']==start]
        for x in transitions_from_start.index.tolist():
            end = transitions_from_start.loc[x, 'a_2']
            # print(end)
            j = abstractions_df.index[abstractions_df['abstraction']==end][0]
            p = transitions_df.loc[(transitions_df['initial_state']==start) & (transitions_df['convergent_state']==end), 'probability'].tolist()[0]
            t_matrix[i][j]=p
    T = numpy.transpose(numpy.array(t_matrix))

    return T
     
def add_local_resiliences_to_dataframe(abstractions_df, transition_matrix):
    '''
    Adds local resiliences to abstractions_df

        Parameter:
            abstractions_df,
            transition_matrix (2D-Array)

        Returns:
    '''
    abstractions_df['local_resilience'] = [transition_matrix[i][i] for i in range(len(transition_matrix))]

def add_global_resiliences_to_dataframe(abstractions_df, transition_matrix):
    '''
    Adds global resiliences to abstractions_df

        Parameter:
            abstractions_df,
            transition_matrix (2D-Array)

        Returns:
    '''
    stationary_transition_matrix = T_inf(transition_matrix)
    n_states = len(transition_matrix)
    global_resiliences = [0 for i in range(n_states)]
    for i in range(n_states):
        initial_state = [(lambda x: 1 if x==i else 0)(x) for x in range(n_states)]
        global_resiliences[i] = numpy.matmul(stationary_transition_matrix, initial_state)[i]
    abstractions_df['global_resilience'] = global_resiliences

def add_reachabilities_to_dataframe(abstractions_df, transition_matrix):
    '''
    Adds local resiliences to abstractions_df

        Parameter:
            abstractions_df,
            transition_matrix (2D-Array)

        Returns:
    '''
    return

def add_markov_properties_to_dataframe(abstractions_df, transitions_df):
    '''
    Adds  markov-properties (local resiliences, local resiliences) to abstractions_df

        Parameter:
            abstractions_df,
            transition_matrix (2D-Array)

        Returns:
    '''
    transition_matrix = transition_matrix_from_dataframes(abstractions_df, transitions_df)
    add_local_resiliences_to_dataframe(abstractions_df, transition_matrix)
    add_global_resiliences_to_dataframe(abstractions_df, transition_matrix)
    return abstractions_df
