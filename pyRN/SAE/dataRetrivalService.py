import json

def get_RNDWs(path, x):
    # 1.1. Reading data from file-system
    with open(path, 'r') as f:
        rndws = json.loads(f.read())

    # 1.2. Converting data
    # 1.2.1. Converting nested dictionary to 2-D arrays
    data = [[list(step.values()) for step in list(json.loads(rndws[i][x]).values())] for i in range(len(rndws))]
    # 1.2.2. Converting 2-D array from {True, False} to {0,1}
    abstractions = []
    for i in range(len(data)):
        a = []
        for j in range(len(data[i])):
            a.append([int(s) for s in data[i][j]])
        abstractions.append(a)
    return abstractions

def get_CPs(path, x):
    # 1.1. Reading data from file-system
    with open(path, 'r') as f:
        rndws = json.loads(f.read())

    # 1.2. Converting data
    # 1.2.1. Converting nested dictionary to 2-D arrays
    complexities = [rndws[i][x] for i in range(len(rndws))]
    return complexities
