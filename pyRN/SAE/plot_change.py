import matplotlib.pyplot as plt
from . import sos

def plot_change(subplot, abstractions, change_function=sos.normalized_hamming_distance, title='', show_indices=True, index_spacing=1, legend=True):
    '''
    In:  mandatory:
         subplot (matplotlib.pyplot subplot)
         abstractions (list), list of binary arrays that represent the abstractions in the random walk

         optional:
         change_function (function), function that calculates the distance between pairs of abstractions
         show_indices (boolean), chooses if indices of transitions will be displayed
         index_spacing (float), modifies spacing between indices
         legend (boolean), chooses if a legend will be displayed

    '''

    # Calculating the size of the sets before and after transitions
    x       = [sos.n_elements(abstractions[i]) for i in range(len(abstractions)-1)]
    y       = [sos.n_elements(abstractions[i+1]) for i in range(len(abstractions)-1)]

    # Plotting a diagonale to indicate numerical stability
    changes = [change_function(abstractions[i], abstractions[i+1]) for i in range(len(abstractions)-1)] 
    min_size = min([sos.n_elements(a) for a in abstractions])
    max_size = max([sos.n_elements(a) for a in abstractions])
    subplot.plot([min_size-1,max_size+1],[min_size-1,max_size+1],color ='black', label ='numerical stability')

    # Divide transitions to smaller lists by the magnitude of the change
    max_change = max(changes)
    if max_change < 0.000001:
        max_change = 0.000001
    do1_3 = (max_change/3)
    do2_3 = 2*do1_3
    d1_3 = (max_change/3)
    d2_3 = 2*d1_3
    small      = [(x[i],y[i],changes[i]*125*(max_size-min_size)/max_change) for i in range(len(x)) if d1_3 >  changes[i]]
    medium     = [(x[i],y[i],changes[i]*125*(max_size-min_size)/max_change) for i in range(len(x)) if d1_3 <= changes[i] and changes[i]<d2_3]
    big        = [(x[i],y[i],changes[i]*125*(max_size-min_size)/max_change) for i in range(len(x)) if         changes[i] >=d2_3]

    # Plotting transition points
    subplot.scatter([b[0] for b in big]   , [b[1] for b in big]   , s=[b[2] for b in big]   , color = 'red')
    subplot.scatter([m[0] for m in medium], [m[1] for m in medium], s=[m[2] for m in medium], color = 'yellow')
    subplot.scatter([s[0] for s in small] , [s[1] for s in small] , s=[s[2] for s in small] , color = 'green')
    subplot.scatter([],[],s=30, label=f'small     (               change < {round((do1_3),3)})', color = 'green')
    subplot.scatter([],[],s=30, label=f'medium ({round((do1_3),3)} < change < {round((do2_3),3)})', color = 'yellow')
    subplot.scatter([],[],s=30, label=f'big         ({round((do2_3),3)} < change                    )', color = 'red')

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

    subplot.set_xlabel("Abstraction size\nbefore transition")
    subplot.set_ylabel("Abstraction size\nafter transition")
    subplot.set_aspect('equal', 'box')
    subplot.set_title(f'{title}\nChange plot')
    if legend:
        subplot.legend(bbox_to_anchor=(1.05, 1), prop={'size': 7.5})