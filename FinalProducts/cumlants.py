# this scripts will read the data of HDF5 file passed as terminal arguments and calculte the cumlants
# and save the results in a new HDF5 file whose name is also passed as argument.

#record starting time
from time import perf_counter
start_time = perf_counter()

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import h5py as hp5
import glob



# Pylinas core modules
import MAS_library as MASL
import smoothing_library as SL

print("Imported Libraries")


# get the arguments passed to the script
# arguments are in order: 
import sys
args = sys.argv
# Assigning the arguments to variables
# path to files
path      = args[1]
file_path = args[2]
save_path = args[3]
# Titles of the plot and verbose outputs
aux = args[4]


# #number of files
n_files = len(glob.glob(path + '*.hdf5'));
print("Processing : "+ aux)
print("Number of files: ", n_files)


#====================================================
# number of particles
n_dm = 752**3
n_gas = 752**3

## reading data from hdf5 files
#====================================================

# first getting auxilary data from one of the files
f = hp5.File(path + file_path + '0.hdf5', 'r')
a = f['Header']. attrs.get('Time') # Scale Factor .
h = f['Header']. attrs.get('HubbleParam') # h
dm_mass = f['Header'].attrs.get('MassTable')[1] # DM mass
BoxSize = f['Header'].attrs.get('BoxSize') # Box Size in cMpc/h
f.close()


# defining arrays for saving data
data   = np.empty((n_dm,3), dtype=np.float32)
dm_mass = np.ones(len(data))*dm_mass
g_data = np.empty((n_gas,3), dtype=np.float32)
g_mass = np.empty(n_gas,dtype=np.float32)
# S_data = np.array([[0,0,0]], dtype=np.float32)
# S_mass = np.array([0], dtype=np.float32)
# BH_data = np.array([[0,0,0]], dtype=np.float32)
# BH_mass = np.array([0], dtype=np.float32)


# index to keep track of the number of particles
i_dm, i_gas = 0, 0;

# reading entire dataset
for i in tqdm(range(n_files), "Reading file:"):
    fileName = path + file_path + str(i) + ".hdf5"
    f = hp5.File(fileName, 'r')
    #DM
    shape = f['PartType1/Coordinates'].shape
    data[i_dm:i_dm+shape[0], :] = f['PartType1/Coordinates'][:]
    i_dm += shape[0]
    # Gas
    shape = f['PartType0/Coordinates'].shape
    g_data[i_gas:i_gas+shape[0], :] = f['PartType0/Coordinates'][:]
    g_mass[i_gas:i_gas+shape[0]]    = f['PartType0/Mass'][:]
    i_gas += shape[0]
    # # reading Star data
    # S_data = np.concatenate((S_data, f['PartType4/Coordinates'][:]), axis=0)
    # S_mass = np.concatenate((S_mass, f['PartType4/Mass'][:]), axis=0)
    # # Reading Blackhole's data
    # BH_data = np.concatenate((BH_data, f['PartType5/Coordinates'][:]), axis=0)
    # BH_mass = np.concatenate((BH_mass, f['PartType5/BH_Mass'][:]), axis=0)
    f.close()

# trimming
# S_data = S_data[1:]
# S_mass = S_mass[1:]
# BH_data = BH_data[1:]
# BH_mass = BH_mass[1:]
print("Snapshots read")
#====================================================

#====================================================
## calulating the cumlants.
#====================================================

# defining function to calutate smoothend field
def smoothed_field(data, BoxSize, Rth, threads=28, Filter = 'Gaussian'):
    """
    Computes the smoothed field from the given positional data
    using Pylians smooting library

    Parameters:
    -----------
    data   : Positional data (numpy array)
    Rth    : Smoothing radius (float)
    Threads: Number of threads (int)
    Filter : Smoothing kernel (string) - Gaussian, Top-Hat
    """
    ## smoothing the data
    #computing the FFT of the filter
    grid = data.shape[0]
    W_k = SL.FT_filter(BoxSize, Rth, grid, Filter , threads)
 
    
    # smooth the field
    field_smoothed = SL.field_smoothing(data, W_k, threads)
    return field_smoothed


# defining the funciton to calculate the cumlants
def cumlant(Np, field, order):
    """
    Calculates the cumlant of a field.
    
    Parameters:
    -----------
    Np   : Number of particles (int)
    field: Desinty field whose cumlant will be calculated (numpy array)
    order: Order of cumlant (int)
    """
    # calculating the cumlant
    mean = np.mean(field)
    field -= mean
    pow = np.power(field, order)
    sum = np.sum(pow)
    err = np.var(pow)
    cumlant = sum/Np
    return [cumlant, err]
    
## calulating CIC fields for gas, DM and All-matter
#====================================================
grid = 752
MAS = 'CIC'
verbose = True

print("\n\nGetting position and delta array...")
# DM positions and delta
pos = np.array(data, dtype=np.float32)
delta = np.zeros((grid,grid,grid), dtype=np.float32)
# Gas positions and delta
g_pos = np.array(g_data, dtype=np.float32)
g_delta = np.zeros((grid,grid,grid), dtype=np.float32)
#Stars and BH Pos
# S_pos = np.array(S_data, dtype=np.float32)
# BH_pos= np.array(BH_data, dtype=np.float32)
# # all matter positions mass and delta
# pos_matter = np.array(np.concatenate((pos, g_pos, S_pos, BH_pos)),dtype=np.float32)
# mass_matter = np.array(np.concatenate((dm_mass, g_mass, S_mass, BH_mass)),dtype=np.float32)
# delta_matter= np.zeros((1504,1504,1504), dtype=np.float32)

print("Position and delta array assigned!!\n")
#========== construct 3D density field=================
# DM
print("Constructing 3D density field of DM")
MASL.MA(pos, delta, BoxSize, MAS, verbose=verbose)
delta /= np.mean(delta, dtype=np.float64);
delta -= 1.0

# gas
print("Constructing 3D density field of gas")
MASL.MA(g_pos, g_delta, BoxSize, MAS,W=g_mass, verbose=verbose)
g_delta /= np.mean(g_delta, dtype=np.float64);
g_delta -= 1.0

# # all matter
# print("Constructing 3D density field of all matter")
# MASL.MA(pos_matter, delta_matter, BoxSize, MAS, W=mass_matter, verbose=verbose)
# delta_matter /= np.mean(delta_matter, dtype=np.float64);
# delta_matter -= 1.0

print("3D density field constructed!!\n")

#====================================================
# Now compute values of smoothing scales 
R = np.linspace((2/752)*BoxSize, 0.1*BoxSize, 50)
# Order of cumlants
O = [2,3,4,5];

# Creating a HDF5 file for storing the data
h5f = hp5.File(save_path+'.hdf5', 'a')


# Looping over R and compute the smoothed field and cumlant
for order in O:
    #creating main group for storing the data
    group = h5f.create_group(str(order))

    # allocating array for storing the cumlants
    Cumlant_dm     = np.empty((len(R),2), dtype=np.float32)
    Cumlant_gas    = np.empty((len(R),2), dtype=np.float32)
    # Cumlant_matter = np.empty((len(R),2), dtype=np.float64)
    idx = 0
    for Rth in tqdm(R, "Smooting scales:"):
        t_inst = perf_counter()
        #calculate the smooth field
        Sf_dm = smoothed_field(delta, BoxSize, Rth)
        Sf_gas = smoothed_field(g_delta, BoxSize, Rth)
        # Sf_matter = smoothed_field(delta_matter, BoxSize, Rth)
        
        # calculating the cumlant
        Cumlant_dm[idx]    = cumlant(dm_mass.shape[0]**3, Sf_dm, order)
        Cumlant_gas[idx]   = cumlant(g_mass.shape[0]**3, Sf_gas, order)
        # Cumlant_matter[idx]= cumlant(mass_matter.shape[0], Sf_matter, order)
        #updating the index
        idx += 1

    # saving the data to the group without saving R
    group.create_dataset('DM', data=Cumlant_dm)
    group.create_dataset('Gas', data=Cumlant_gas)
    # group.create_dataset('All-Matter', data=Cumlant_matter)

h5f.close()
print("Data saved to HDF5 file")
print("Computation complete!!!")
print("\n Time taken: ", round((perf_counter()-start_time)/60), " minutes")
