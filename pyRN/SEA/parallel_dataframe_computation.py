# pyRN imports
from pyRN import pyRN
from pyRN.SEA import markov
from pyRN.SEA import newdataframes
from pyRN.SEA import sos

# built-in imports
from contextlib import redirect_stdout
import concurrent.futures
import io
import math
import os
import pickle
import shutil
import sys
import time

import io

def move_files_to_folders(path):
    '''
    For every file in a folder it creates a folder with
    the name of the file and moves the file there
    '''
    for filename in os.listdir(path):
        file_path = os.path.join(path, filename)
        if os.path.isfile(file_path):
            folder_path = os.path.join(path, filename.replace(filename[filename.rfind('.'):], ''))
            os.makedirs(folder_path)
            shutil.move(file_path, folder_path)

def pyRN_object_from_file(path):
    '''
    Initializes and returns a pyRN-object for a reaction network file saved at path

        Parameters:
            path (string), file path to a reaction network (.txt or .xml)

        Returns:
            pyRN-object
    '''
    try:
        if path.endswith('.txt'):
            return pyRN.setFromText(path)
        if path.endswith('.xml'):
            return pyRN.setFromSbml(path)
    except:
        print(f'Could not load network from {path}')
    
def pkl(object, path):
    '''
    Store python object using pickle at given path
    
        Parameters:
            object (python object)
            path   (string), filepath
    '''
    essig = open(path, 'wb')
    pickle.dump(object, essig)
    essig.close()

def depkl(path):
    '''
    Load python object from given path using pickle
    
        Parameters:
            path   (string), filepath
        
        Returns:
            object (python object)
    '''
    essig = open(path, 'rb')
    with essig as f:
        object = pickle.load(f)
    return object

def set_and_store(path,
                  max_pert_size = math.inf,
                  conn=True,
                  closure=True,
                  include_empty_set=True):
    '''
    Initializes a pyRN object for a reaction network file saved at path,
    sets it's generators, synergetic structure and SimpleTransSpDf

        Parameters:
            path (string), file path to a reaction network (.txt or .xml)
            max_pert_size (int),
            conn (boolean)
            closure boolean,
            include_empty_set=False

        Returns:
            pyRN-object
    '''
    start = time.time()
    print(f'{path}: start set_and_store')
    RN = pyRN_object_from_file(path)
    f = io.StringIO()
    with redirect_stdout(f):
        RN.setSpConnMat()
    print(f'{path}: connectivity matrix set')
    with redirect_stdout(f):
        RN.setGenerators()
    print(f'{path}: generators set')
    with redirect_stdout(f):
        RN.setSynStr()
    print(f'{path}: synergetic structure set')
    pert_size = min(max_pert_size, len(RN.SpIdStrArray))
    RN.setSimpleTransDict(orglist=RN.SynStrOrgListBtArray,
                          pert_type="species",
                          pert_size=pert_size,
                          conn=conn,closure=closure,
                          include_empty_set=include_empty_set)
    print(f'{path}: SimpleTransSpDf set')
    pkl(RN, path.replace(path[path.rfind('.')+1:], 'pickle'))
    print(f'{path}: set_and_store completed after {round(time.time()-start,2)} seconds')

def multithread_set_and_store(file_paths, max_pert_size, threads):
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit the function with each parameter to the executor
        results = [executor.submit(set_and_store, file_path, max_pert_size) for file_path in file_paths]
        # Wait for all the results to be finished
        concurrent.futures.wait(results)

def get_file_paths(path, extension=None):
    """
    Returns a list of file paths in the directory located at `path` and its subdirectories.
    
        Parameters:
            path (string): A string representing the path to the directory to search.
            extension: An optional string representing the file extension to search for. If not provided, all files will be returned.
        Return
            A list of strings representing the file paths.
    """
    file_paths = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if extension is None or file.endswith(extension):
                file_paths.append(os.path.join(root, file))
    return file_paths

def calculate_dataframes(SimpleTransSpDf, max_perturbation_size, folder_path):
    # Initialize dataframes
    abstractions_df = newdataframes.initialize_abstractions_df(SimpleTransSpDf)
    transitions_df  = newdataframes.initialize_transitions_df(abstractions_df)
    # Calculate transition probabilities
    allowed = lambda p: sos.n_elements(list(p)) <= max_perturbation_size
    TransDf = SimpleTransSpDf.copy()
    TransDf = SimpleTransSpDf[SimpleTransSpDf['perturbation'].apply(allowed)]
    newdataframes.add_probabilities_to_transitions_df_2_1(TransDf, transitions_df)
    newdataframes.fix_transition_probabilities_for_all_initial_states(transitions_df)
    # Calculate markov properties
    abstractions_df = markov.add_markov_properties_to_dataframe(abstractions_df, transitions_df)
    # Add additional information
    newdataframes.add_number_of_species(abstractions_df)
    newdataframes.add_number_of_species(abstractions_df)
    newdataframes.add_size_difference(transitions_df)
    # Store data
    pkl(abstractions_df, f'{folder_path}\\abstractions_df_pert_size_{max_perturbation_size}.pickle')
    pkl(transitions_df, f'{folder_path}\\transitions_df_pert_size_{max_perturbation_size}.pickle')

def calculate_all_dataframes(file_path, max_perturbation_size):
    RN = depkl(file_path)
    folder_path = file_path.replace(file_path[file_path.rfind('\\')+1:], '')
    for max_perturbation_size in range(1, max_perturbation_size+1):
        print(f'{folder_path}: calculate dataframes with max_pert_size {max_perturbation_size}')
        calculate_dataframes(RN.SimpleTransSpDf, max_perturbation_size, folder_path)

def multithread_calculate_all_dataframes(file_paths, max_perturbation_size):
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit the function with each parameter to the executor
        results = [executor.submit(calculate_all_dataframes, file_path, max_perturbation_size) for file_path in file_paths]
        # Wait for all the results to be finished
        concurrent.futures.wait(results)

##################################################
#------------------------------------------------#
path = r"C:\Users\simon\Downloads\small_nets2"
path = r'C:\Users\simon\Downloads\non_connected_example'
max_perturbation_size = 2
threads = 2                
#------------------------------------------------#
##################################################

start = time.time()
# Put every reaction network in it's own folder
move_files_to_folders(path)
# Get paths to all networks
file_paths = get_file_paths(path)
# Initialize networks,
# set connectivity matrix
# set generators
# set synergetic structure
# set SimpleTransSpDf
multithread_set_and_store(file_paths, max_perturbation_size, threads)
# Get paths to all pyRN-objects
file_paths = get_file_paths(path, '.pickle')
# Load pyRN-objects, calculate and save 
# abstractions_df and transitions_df
# for all possible 
sys.stdout = sys.__stdout__
multithread_calculate_all_dataframes(file_paths, max_perturbation_size)
duration = (time.time()-start)
minutes  = int(duration/60)
seconds  = int(duration-minutes*60)
print(f'\nCompleted task in {minutes} minutes and {seconds} seconds')