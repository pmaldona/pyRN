import numpy
import copy

def T_inf(T):
    for i in range(10):
        T = numpy.matmul(T,T)
    return T

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
    In:  df (pandas dataframe), with at least three columns: a_1, a_2 and count
    Out: df (pandas dataframe)
    Adds a column to thdf with the transitionprobabilities in the markov-model
    '''
    df = df.assign(probability=[0 for i in range(df.shape[0])])
    starts = df['a_1'].unique()
    for start in starts:
        transitions_from_start = df.loc[df['a_1']==start]
        n = transitions_from_start['counts'].sum()
        ends = transitions_from_start['a_2']
        for end in ends:
            df.loc[(df['a_1']==start) & (df['a_2']==end),'probability'] = df.loc[(df['a_1']==start) & (df['a_2']==end),'counts']/n
    return df

def transition_matrix_from_dataframes(abstractions_df, transitions_df):
    '''
    In:  df (pandas dataframe), with at least three columns: a_1, a_2 and either count or probability
    Out: Transition mastrix of Markov-model
    '''
    if ('probability' not in transitions_df.columns):
        df = add_transition_probbilities_to_dataframe(df)
    
    n = abstractions_df.shape[0]
    t_matrix = [ [0 for j in range(n)] for i in range(n)]
    for i in range(n):
        start = abstractions_df.loc[i,'a']
        transitions_from_start = transitions_df.loc[transitions_df['a_1']==start]
        for x in transitions_from_start.index.tolist():
            end = transitions_from_start.loc[x, 'a_2']
            print(end)
            j = abstractions_df.index[abstractions_df['a']==end][0]
            p = transitions_df.loc[(transitions_df['a_1']==start) & (transitions_df['a_2']==end), 'probability'].tolist()[0]
            t_matrix[i][j]=p
    T = numpy.transpose(numpy.array(t_matrix))

    return T
     
def add_markov_properties_to_dataframe(abstractions_df, transitions_df):

    T      = transition_matrix_from_dataframes(abstractions_df, transitions_df)
    T_stat = T_inf(T)

    reachabilities    = []
    maintainabilities = []
    stricts           = []

    n = abstractions_df.shape[0]

    sum_inits = abstractions_df['n_init'].sum()
    init      = [(abstractions_df.loc[i, 'n_init']/sum_inits) for i in range(n)]

    for i in range(n):
        state    = [0 for i in range(n)]
        state[i] = 1
        reachabilities.append(reachability(init, i, T, 0.2))
        maintainabilities.append(numpy.matmul(T_stat, state)[i])
        s = list(transitions_df.loc[(transitions_df['a_1']==abstractions_df.loc[i, 'a']) & (transitions_df['a_2']==abstractions_df.loc[i, 'a']), 'probability'])
        if len(s) > 0:
            stricts.append(s[0])
        else:
            stricts.append(0)
    
    abstractions_df = abstractions_df.assign(reachability=reachabilities, maintainability=maintainabilities, strict_maintainability=stricts)

    return abstractions_df
