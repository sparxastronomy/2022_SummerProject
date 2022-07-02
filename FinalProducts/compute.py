# python file to compute the Plots and Results 

#record starting time
from time import perf_counter
start_time = perf_counter()

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import h5py as hp5
import glob

from matplotlib import colors
from matplotlib import cm
from sctriangulate.colors import build_custom_continuous_cmap

# Pylinas core modules
import Pk_library as PKL
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
S_data = np.array([[0,0,0]], dtype=np.float32)
S_mass = np.array([0], dtype=np.float32)
BH_data = np.array([[0,0,0]], dtype=np.float32)
BH_mass = np.array([0], dtype=np.float32)
BH_Mdot = np.array([0], dtype=np.float32) #accretion rate


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
    # reading Star data
    S_data = np.concatenate((S_data, f['PartType4/Coordinates'][:]), axis=0)
    S_mass = np.concatenate((S_mass, f['PartType4/Mass'][:]), axis=0)
    # Reading Blackhole's data
    BH_data = np.concatenate((BH_data, f['PartType5/Coordinates'][:]), axis=0)
    BH_mass = np.concatenate((BH_mass, f['PartType5/BH_Mass'][:]), axis=0)
    BH_Mdot = np.concatenate((BH_Mdot, f['PartType5/BH_Mdot'][:]), axis=0)
    f.close()

# trimming
S_data = S_data[1:]
S_mass = S_mass[1:]
BH_data = BH_data[1:]
BH_mass = BH_mass[1:]
BH_Mdot = BH_Mdot[1:]
print("Snapshots read")
#====================================================

#===============================================================================================
print('Plotting density fields')
new_cmap = build_custom_continuous_cmap([9,9,14], [28, 70, 138], [227, 172, 70], [255,255,255]);
gas_cmap = build_custom_continuous_cmap([3,3,3], [88, 64, 161], [207, 126, 83], [247, 246, 181]);

## Plotting the DM density
plt.figure(figsize=(6.5,6.5), dpi=100)
ax  = plt.gca()
ax.set_aspect('equal')
plt.hist2d(data[:,0], data[:,1], norm =colors.LogNorm(), cmap=new_cmap, bins=2048);

plt.title('Z-axis projection of DM particles: '+aux)


plt.xlabel('x [cMpc/h]')
plt.ylabel('y [cMpc/h]')


plt.savefig(save_path+'_DM.jpg', dpi=400, bbox_inches='tight')

## plotting gas density field
plt.figure(figsize=(6.5,6.5), dpi=100)
ax  = plt.gca()
ax.set_aspect('equal')
plt.hist2d(g_data[:,0], g_data[:,1], norm =colors.LogNorm(), cmap=gas_cmap, bins=2048);
plt.title('Z-axis projection of Gas particles: '+aux)
plt.xlabel('x [cMpc/h]')
plt.ylabel('y [cMpc/h]')
plt.savefig(save_path+'_Gas.jpg', dpi=400, bbox_inches='tight')
#===============================================================================================


#==============================================================================================================
## Computing 2pCF for DM and Gas
#==============================================================================================================
print("\nComputing 2pCF for: "+aux)
## computing xi(r)
grid = 350
# BoxSize = 205 #already obtained from the first file
verbose = True
Np = len(data)
MAS = 'CIC'


print("Getting position and delta array")
# DM positions and delta
pos = np.array(data, dtype=np.float32)
delta = np.zeros((grid,grid,grid), dtype=np.float32)
# Gas positions and delta
g_pos = np.array(g_data, dtype=np.float32)
g_delta = np.zeros((grid,grid,grid), dtype=np.float32)
# Stars and BH positions, no need of delta for these
S_pos = np.array(S_data, dtype=np.float32)
BH_pos= np.array(BH_data, dtype=np.float32)
BH_delta = np.zeros((grid,grid,grid), dtype=np.float32)

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

# BH
print("Constructing 3D density field of BlackHoles")
MASL.MA(BH_pos, BH_delta, BoxSize, MAS,W=BH_Mdot, verbose=verbose)
BH_delta /= np.mean(BH_delta, dtype=np.float64);
BH_delta -= 1.0
#=====================================================

#=====================================================
### compute the correlation function
print("Computing correlation function")
axis = 0
threads=14
CF     = PKL.Xi(delta, BoxSize, MAS, axis, threads)  # DM
g_CF   = PKL.Xi(g_delta, BoxSize, MAS, axis, threads) # gas
BH_G_CF= PKL.XXi(BH_delta,g_delta, BoxSize, MAS, axis, threads) # BH and Gas Cross-correlation


## get the attributes
print("Getting attributes")
# DM
r_dm      = CF.r3D      #radii in Mpc/h
xi0_dm    = CF.xi[:,0]  #correlation function (monopole)
# gas
r_gas     = g_CF.r3D      #radii in Mpc/h
xi0_gas   = g_CF.xi[:,0]  #correlation function (monopole)
# BH and gas
r_BH_gas  = BH_G_CF.r3D      #radii in Mpc/h
xi0_BH_gas= BH_G_CF.xi[:,0]  #correlation function (monopole)

## filtering data to have only non-nan values
# DM
xi0_filtered_dm = xi0_dm[~np.isnan(np.log10(xi0_dm))]
r_filterd_dm    = r_dm[~np.isnan(np.log10(xi0_dm))]
# gas
xi0_filtered_gas = xi0_gas[~np.isnan(np.log10(xi0_gas))]
r_filterd_gas    = r_gas[~np.isnan(np.log10(xi0_gas))]
# BH and gas
xi0_filtered_BH_gas = xi0_BH_gas[~np.isnan(np.log10(xi0_BH_gas))]
r_filterd_BH_gas    = r_BH_gas[~np.isnan(np.log10(xi0_BH_gas))]

#=====================================================

#=====================================================
### Computing the power spectrum
print("Computing power spectrum for: "+aux)
Pk    = PKL.Pk(delta, BoxSize,axis, MAS, threads) # DM
g_Pk  = PKL.Pk(g_delta, BoxSize,axis, MAS, threads) # gas

## get the attributes
print("Getting attributes")
# DM
k_dm      = Pk.k3D      #k in h/Mpc
Pk0_dm    = Pk.Pk[:,0]  #power spectrum (monopole)
# gas
k_gas     = g_Pk.k3D      #k in h/Mpc
Pk0_gas   = g_Pk.Pk[:,0]  #power spectrum (monopole)

## filtering data to have only non-nan values
# DM
Pk0_filtered_dm = Pk0_dm[~np.isnan(np.log10(Pk0_dm))]
k_filterd_dm    = k_dm[~np.isnan(np.log10(Pk0_dm))]
# gas
Pk0_filtered_gas = Pk0_gas[~np.isnan(np.log10(Pk0_gas))]
k_filterd_gas    = k_gas[~np.isnan(np.log10(Pk0_gas))]
#=====================================================

#=====================================================
# Saving the data to a HDF5 file
print("Saving the data to a HDF5 file")
h5f = hp5.File(save_path+'.hdf5', 'a')
# main groups: 2pCF, Pk, Cross-correlation
m1 = h5f.create_group('2pCF')
m2 = h5f.create_group('Pk')
m3 = h5f.create_group('CC')
# Subgroups for each group: DM, gas
m1.create_group('DM')
m1.create_group('gas')
m2.create_group('DM')
m2.create_group('gas')
m3.create_group('BH-Gas')
## Data Saving
# 2pCF
m1['DM'].create_dataset('xi0', data=xi0_filtered_dm)
m1['DM'].create_dataset('r', data=r_filterd_dm)
m1['gas'].create_dataset('xi0', data=xi0_filtered_gas)
m1['gas'].create_dataset('r', data=r_filterd_gas)

# Pk
m2['DM'].create_dataset('Pk0', data=Pk0_filtered_dm)
m2['DM'].create_dataset('k', data=k_filterd_dm)
m2['gas'].create_dataset('Pk0', data=Pk0_filtered_gas)
m2['gas'].create_dataset('k', data=k_filterd_gas)
# Cross-correlation
m3['BH-Gas'].create_dataset('xi0', data=xi0_filtered_BH_gas)
m3['BH-Gas'].create_dataset('r', data=r_filterd_BH_gas)
h5f.close()

#=====================================================

#==============================================================================================================
# All matter 2pCF and Pk
#==============================================================================================================
del delta, g_delta

# all matter positions
pos_matter = np.array(np.concatenate((pos, g_pos, S_pos, BH_pos)),dtype=np.float32)
# all matter masses
mass_matter = np.array(np.concatenate((dm_mass, g_mass, S_mass, BH_mass)),dtype=np.float32)

# all matter density field
print("Constructing all matter 3D density field")
delta = np.zeros((512,512,512), dtype=np.float32)
MASL.MA(pos, delta, BoxSize, MAS,W=mass_matter, verbose=verbose)
delta /= np.mean(delta, dtype=np.float64);
delta -= 1.0

#=====================================================
# all matter 2pCF and Power Spectrum
print("Computing all matter 2pCF and Power Spectrum")
CF_all_matter     = PKL.Xi(delta, BoxSize, MAS, axis, threads) 
Pk_all_matter     = PKL.Pk(delta, BoxSize,axis, MAS, threads) 

# get the attributes
print("Getting attributes")
# all matter 2pCF
all_matter_r      = CF.r3D      #radii in Mpc/h
all_matter_xi0    = CF.xi[:,0]  #correlation function (monopole)

# all matter Power Spectrum
all_matter_k      = Pk.k3D      #k in h/Mpc
all_matter_Pk0    = Pk.Pk[:,0]  #power spectrum (monopole)

#filtering data to have only non-nan values
# all matter 2pCF
all_matter_xi0_filtered = all_matter_xi0[~np.isnan(np.log10(all_matter_xi0))]
all_matter_r_filtered    = all_matter_r[~np.isnan(np.log10(all_matter_xi0))]
# all matter Power Spectrum
all_matter_Pk0_filtered = all_matter_Pk0[~np.isnan(np.log10(all_matter_Pk0))]
all_matter_k_filtered    = all_matter_k[~np.isnan(np.log10(all_matter_Pk0))]
#=====================================================

#=====================================================
# Saving the data to the previous HDF5 file
print("Saving the data to a HDF5 file")
h5f = hp5.File(save_path+'.hdf5', 'a')
# main groups: 2pCF, Pk
m1 = h5f['2pCF']
m2 = h5f['Pk']
m3 = h5f['CC']
# creating a new group for all matter
m1.create_group('all_matter')
m2.create_group('all_matter')
## Saving the data
# 2pCF
m1['all_matter'].create_dataset('xi0', data=all_matter_xi0_filtered)
m1['all_matter'].create_dataset('r', data=all_matter_r_filtered)
# Pk
m2['all_matter'].create_dataset('Pk0', data=all_matter_Pk0_filtered)
m2['all_matter'].create_dataset('k', data=all_matter_k_filtered)
print(m1.keys(), m2.keys(), m3.keys())
h5f.close()
print("Data Saved!!")
#==============================================================================================================

print("Done!!")
print("Time taken: ", round((perf_counter()-start_time)/60,3), " minutes")