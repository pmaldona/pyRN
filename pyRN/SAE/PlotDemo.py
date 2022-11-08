import matplotlib.pyplot as plt
import json
import plot_raw
import plot_hasse
import plot_change
import sos

########################################################################################
# EXPLANATION:                                                                         #
#                                                                                      #
# The idea was to make the plotting functions as similar to one another as possible.   #
# In order to have a more convinient way to compare different sets of data and aspects #
# of data at the same time, the user can initialize a plot with subplots and use       #
# the provided scripts to draw the plots.                                              #
#                                                                                      #
# The representation of the state abstractions are lists consisting only of 0 and 1    #
# The state abstraction {x1, x2, x5} for a system with molecules {x1, x2, x3, x4, x5}  #
# would be represented as [1,1,0,0,1].                                                 #
# The methods plot_raw.plot_raw() and plot_change.plot_change() both get a subplot and #
# a list of state abstractions.                                                        #
# The method plot_hasse.plot_hasse() also gets a subplot but it gets a list of lists   #
# abstractions. This one can have multiple trajectories in the same plot for better    #
# comparability.                                                                       #
# All methods have optional parameters.                                                #
########################################################################################



########################################################################################
# START USER-INPUT (path to .json-file)                                                #
path = 
# END USER-INPUT   (path to .json-file)                                                #
########################################################################################

# 1. Reading and converting data

# 1.1. Reading data from file-system
with open(path, 'r') as f:
    rndws = json.loads(f.read())

# 1.2. Converting data
# 1.2.1. Converting nested dictionary to 2-D arrays
data = [[list(step.values()) for step in list(json.loads(rndws[i]['a']).values())] for i in [0,1]]
# 1.2.2. Converting 2-D array from {True, False} to {0,1}
abstractions = []
for i in range(len(data)):
    a = []
    for j in range(len(data[i])):
        a.append([int(s) for s in data[i][j]])
    abstractions.append(a)

# 2. Plotting

# 2.1 Plot initialization
figure = plt.figure()
ax00 = plt.subplot2grid((3,2),(0,0))
ax01 = plt.subplot2grid((3,2),(0,1))
ax1  = plt.subplot2grid((3,2),(1,0), colspan=2)
ax20 = plt.subplot2grid((3,2),(2,0))
ax21 = plt.subplot2grid((3,2),(2,1))

# 2.2. Plotting
plot_raw.plot_raw(ax00, abstractions[0])
plot_raw.plot_raw(ax01, abstractions[1])
plot_hasse.plot_hasse(ax1, abstractions)
plot_change.plot_change(ax20, abstractions[0], legend=False, index_spacing=1.5)
plot_change.plot_change(ax21, abstractions[1], index_spacing=1.5)

# 2.3. Adjust subplots (Works only for this specific cofiguration)
plt.subplots_adjust(left=0.25, bottom=None, right=0.75, top=None, wspace=0, hspace=0.3)

# 2.4. Display plots
plt.show()