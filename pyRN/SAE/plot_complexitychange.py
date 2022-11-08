import json
import matplotlib as mp
import matplotlib.pyplot as plt
from . import sos
import pandas

def plot_complexitychange(subplot, abstractions, abstractions_df, title='', show_indices=True, index_spacing=1):
    '''
    In:  mandatory:
         subplot (matplotlib.pyplot subplot)
         abstractions (list), list of binary arrays that represent the abstractions in the random walk
         abstractions_df (pandas dataframe)

         optional:
         change_function (function), function that calculates the distance between pairs of abstractions
         show_indices (boolean), chooses if indices of transitions will be displayed
         index_spacing (float), modifies spacing between indices
         legend (boolean), chooses if a legend will be displayed

    '''

    complexities = []
    strict_maintainabilities = []
    maintainabilities = []
    for a in abstractions:
        complexities.append(float(abstractions_df.loc[abstractions_df['a']==str(a),'complexity']))
        strict_maintainabilities.append(float(abstractions_df.loc[abstractions_df['a']==str(a),'strict_maintainability']))
        maintainabilities.append(float(abstractions_df.loc[abstractions_df['a']==str(a),'maintainability'])*1000)  
    strict_maintainabilities.pop(0)
    maintainabilities.pop(0)

    abstraction_sizes = [sos.n_elements(a)/len(abstractions[0]) for a in abstractions]
    abstraction_sizes.pop(0)
    
        
    # Calculating the size of the sets before and after transitions
    x       = [complexities[i]   for i in range(len(complexities)-1)]
    y       = [complexities[i+1] for i in range(len(complexities)-1)]

    # Plotting a diagonale to indicate numerical stability
    min_size = min(complexities)
    max_size = max(complexities)
    subplot.plot([min_size-1,max_size+1],[min_size-1,max_size+1],color ='black', label ='numerical stability')

    # Plotting transition points
    subplot.scatter(x, y, s=maintainabilities, alpha=strict_maintainabilities, cmap='RdYlGn', c=abstraction_sizes)

    # Plotting connections
    subplot.plot(x,y)

    # Add numbers to identify order the of transitions in a way that numbers do not overlap and also add labels with the change in that transition
    if show_indices:
        x_      = [0 for xi in x]
        y_      = [0 for yi in y]
        dist    = 0.05*(max_size-min_size)*index_spacing
        for i in range(len(x)):                 # Iterate over all coordinates
            occurrences = 0                     # Number of previous transitions with the same coordinates
            for j in range(i):
                if x[i]==x[j] and y[i]==y[j]:
                    occurrences += 1
            if x[i]<y[i]:
                x_[i] = x[i]-dist*(occurrences+1)
                y_[i] = y[i]+dist*(occurrences+1)
            else:
                x_[i] = x[i]+dist*(occurrences+1)
                y_[i] = y[i]-dist*(occurrences+1)   
        for i in range(len(x_)):
            subplot.text(x_[i],y_[i],f'{i}', fontsize = 10, color='black')

    subplot.set_xlabel("Complexity before transition")
    subplot.set_ylabel("Complexity after transition")
    subplot.set_aspect('equal', 'box')
    subplot.set_title(f'{title}\nComplexity change plot')
