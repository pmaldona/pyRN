import json
from . import markov
import numpy
import pandas
from . import sos
import time

def transitions_dataframe_from_SimpleRWDict(SimpleRWDict, w=None, abstraction_type='c'):
    '''
    In:  SimpleRWDict: Dictonary of simple random walk 
         w: Kyes of SimpleRWDict Dictonary considered, if w=None, it use all the 
         Dictionary
         abstraction_type: corrspond to the type of abstraction to be considered, if 
         its 'all' used both perturbation and convegence states
    Out: pandas dataframe with transitions
    '''  

    # Find all observed transitions
    transitions    = []
    counts         = []
    
    if w is None:
        w=SimpleRWDict.keys()
        
    for i in w:
        
        if abstraction_type=='all':
            
            abs_1 = SimpleRWDict[i]['p']# Converts nested dictionary to 2-dimensional array
            abs_2 = SimpleRWDict[i]['c']# Converts nested dictionary to 2-dimensional array 
            abstractions = [item for sublist in zip(abs_1, abs_2) for item in sublist]
        else:
            abstractions = SimpleRWDict[i][abstraction_type]# Converts nested dictionary to 2-dimensional array
        
        for i in range(len(abstractions)-1):
            a_1 = str([int(s) for s in abstractions[i]])
            a_2 = str([int(s) for s in abstractions[i+1]])
            if (a_1,a_2) not in transitions:
                transitions.append((a_1,a_2))
                counts.append(1)
            else:
                counts[transitions.index((a_1,a_2))] +=1

    # Create pandas dataframe and store values there
    data = {'a_1': [transition[0] for transition in transitions], 'a_2': [transition[1] for transition in transitions], 'counts': counts}
    df = pandas.DataFrame(data)

    # Calculate and add changes to dataframe
    nrows = df.shape[0]
    df = df.assign(change_local=[0 for i in range(nrows)], change_global=[0 for i in range(nrows)], complexity_change=[0 for i in range(nrows)])
    for i in range(nrows):
        df.loc[i, 'change_local']     = sos.change_local(sos.from_string(df.loc[i, 'a_1']), sos.from_string(df.loc[i, 'a_2']))
        df.loc[i, 'change_global']    = sos.normalized_hamming_distance(sos.from_string(df.loc[i, 'a_1']), sos.from_string(df.loc[i, 'a_2']))

    

    # Calculate and add transition probabilities in Markov-model
    df = markov.add_transition_probbilities_to_dataframe(df)

    return df

def abstractions_dataframe_from_SimpleRWDict_and_transitions_dataframe(SimpleRWDict,transitions_df,w=None,abstraction_type='c'):
    '''
    In:  SimpleRWDict: Dictonary of simple random walk 
         w: Kyes of SimpleRWDict Dictonary considered, if w=None, it use all the 
         Dictionary
         abstraction_type: corrspond to the type of abstraction to be considered, if 
         its 'all' used both perturbation and convegence states 
         transitions_df: Transition DataFrame obtaind form the 
         transitions_dataframe_from_SimpleRWDict function 
    Out: pandas dataframe with abstrac
    '''

    # Find unique abstractions and sort them by the number of species
    
    if w is None:
        w=SimpleRWDict.keys()
    
    abstractions = []
    for i in w:
        
        if abstraction_type=='all':
            
            abs_1 = SimpleRWDict[i]['p']# Converts nested dictionary to 2-dimensional array
            abs_2 = SimpleRWDict[i]['c']# Converts nested dictionary to 2-dimensional array 
            abstractions += [item for sublist in zip(abs_1, abs_2) for item in sublist]
        else:
            abstractions += SimpleRWDict[i][abstraction_type]# Converts nested dictionary to 2-dimensional array
    
    abstractions=list(map(lambda x: str([int(s) for s in x]),abstractions))
    abstractions = list(set(abstractions)) # Filters duplicates
    abstractions = [sos.from_string(a) for a in abstractions]
    abstractions.sort(key=sos.n_elements)
    n_species = [sos.n_elements(a) for a in abstractions]
    abstractions = [str(a) for a in abstractions]
    n = len(abstractions)

    # Determine how often abstractions occur as initial abstractions
    initials = [0 for i in range(n)]
    for i in w:
        
        if abstraction_type=='c':
            init=str([int(s) for s in SimpleRWDict[i]['c'][0]])
        else:
            init=str([int(s) for s in SimpleRWDict[i]['p'][0]])
            
        for i,a in enumerate(abstractions):
            if abstractions[i]==init:
                initials[i]+=1
    
    abstractions_df = pandas.DataFrame({'a': abstractions, '#species': n_species, 'n_init': initials})

    # Determine complexity of each abstraction
    abstractions_df = abstractions_df.assign(complexity=[0 for i in range(n)])
    for i in range(n):
        for j in w:
            
            if abstraction_type=='all':
                abs_1 = SimpleRWDict[j]['p']# Converts nested dictionary to 2-dimensional array
                abs_2 = SimpleRWDict[j]['c']# Converts nested dictionary to 2-dimensional array 
                rndw_abstractions = [item for sublist in zip(abs_1, abs_2) for item in sublist]
                
                com_1 = SimpleRWDict[j]['cp']# Converts nested dictionary to 2-dimensional array
                com_2 = SimpleRWDict[j]['cc']# Converts nested dictionary to 2-dimensional array 
                complexity = [item for sublist in zip(com_1, com_2) for item in sublist]
                
                
            else:
                 rndw_abstractions = SimpleRWDict[j][abstraction_type]# Converts nested dictionary to 2-dimensional array
                 complexity = SimpleRWDict[j]['c'+abstraction_type]
                 
            rndw_abstractions = [str([int(s) for s in a]) for a in rndw_abstractions]

            for k, abstraction in enumerate(rndw_abstractions):
                #print(rndws[j]['c'+x][k])
                if(abstractions_df.loc[i, 'a']==abstraction):
                    abstractions_df.loc[i, 'complexity']=complexity[k]


    # Indegree
    abstractions_df = abstractions_df.assign(indegree=[0 for i in range(n)])
    for i in range(n):
        start = abstractions_df.loc[i, 'a']
        transitions_to_end = transitions_df.loc[transitions_df['a_2']==start]
        abstractions_df.loc[i, 'indegree']=transitions_to_end.shape[0]

    # Outdegree
    abstractions_df = abstractions_df.assign(outdegree=[0 for i in range(n)])
    for i in range(n):
        start = abstractions_df.loc[i, 'a']
        transitions_from_start = transitions_df.loc[transitions_df['a_1']==start]
        abstractions_df.loc[i, 'outdegree']=transitions_from_start.shape[0]

    # Calculate and add markov properties
    abstractions_df = markov.add_markov_properties_to_dataframe(abstractions_df, transitions_df)

    return abstractions_df

def add_abstraction_indices_to_transitions_dataframe(abstractions_df, transitions_df):
    a1 = []
    a2 = []
    for i in range(transitions_df.shape[0]):
        a1.append(list(abstractions_df.index[abstractions_df['a']==transitions_df.loc[i, 'a_1']])[0])
        a2.append(list(abstractions_df.index[abstractions_df['a']==transitions_df.loc[i, 'a_2']])[0])
    transitions_df.insert(0, 'a2', a2)
    transitions_df.insert(0, 'a1', a1)
    return transitions_df

def add_plasticity_and_smoothness_to_transitions_dataframe(rndws, x, transitions_df, intervall_size):
    return

#add_complexity_changes(abstractions_df, transitions_df):

def dataFromSimpleRwDict(SimpleRWDict, w=None, abstraction_type='c', save=False, file_path=None):
    '''
        In:  SimpleRWDict: Dictonary of simple random walk 
             w: Kyes of SimpleRWDict Dictonary considered, if w=None, it use all the 
             Dictionary
             abstraction_type: corrspond to the type of abstraction to be considered, if 
             its 'all' used both perturbation and convegence states 
             transitions_df: Transition DataFrame obtaind form the 
             transitions_dataframe_from_SimpleRWDict function          
             save      (boolean), if True: Stores the dataframes in the same directory as the input file
             file_path (string),  path were dataframes are to be stored
             
    '''

    transitions_df  = transitions_dataframe_from_SimpleRWDict(SimpleRWDict, w=w, abstraction_type=abstraction_type)
    abstractions_df = abstractions_dataframe_from_SimpleRWDict_and_transitions_dataframe(SimpleRWDict, transitions_df, w=w, abstraction_type=abstraction_type)
    transitions_df  = add_abstraction_indices_to_transitions_dataframe(abstractions_df, transitions_df)
    #transitions_df = add_complexity_changes(abstractions_df, transitions_df)

    if(save):
        transitions_df.to_csv(file_path+'_transitions.csv')
        abstractions_df.to_csv(file_path+'_abstractions.csv')

    return (abstractions_df, transitions_df)