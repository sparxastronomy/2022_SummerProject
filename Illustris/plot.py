import matplotlib.pyplot as plt
import numpy as np

import h5py as hp5
from tqdm import  tqdm

# modules for custom colormaps
from matplotlib import colors
from matplotlib import cm
from sctriangulate.colors import build_custom_continuous_cmap

# Pylinas core modules
import Pk_library as PKL
import MAS_library as MASL
print("All modules imported")

# reading data hdf5 snapshots
basePath = "TNG300-3/snap_099."

data = np.array([[0,0,0]], dtype=np.float32)

for i in tqdm(range(16), desc="Reading data"):
    fileName = basePath + str(i) + ".hdf5"
    f = hp5.File(fileName, 'r')
    temp = np.array(f['PartType1/Coordinates'][:], dtype=np.float32)
    data = np.append(data, temp, axis=0)
    del f, temp

# plotting
print("Plotting Matter distribution")
new_cmap = build_custom_continuous_cmap([9,9,14], [25, 65, 130], [227, 172, 70], [255,255,255])

plt.figure(figsize=(9,5.5), dpi=100)
plt.hist2d(data[:,0]/1000, data[:,1]/1000, norm =colors.LogNorm(), cmap=new_cmap, bins=2048);

plt.title('Z-axis projection of DM particles from TNG300-3 ($z=0$)\n Total DM Particles: '+str(len(data)))

plt.xlabel('x [cMpc/h]')
plt.ylabel('y [cMpc/h]')

plt.savefig('TNG300_DM.jpg', dpi=300, bbox_inches='tight')
plt.show()

# ## computing xi(r)
# print("Computing 2pCF")
# grid = 512
# BoxSize = 302.6
# verbose = True
# Np = len(data)
# MAS = 'CIC'

# print("Getting position and delta array")
# pos = np.array(data/1000, dtype=np.float32)
# delta = np.zeros((grid,grid,grid), dtype=np.float32)

# # construct 3D density field
# print("Constructing 3D density field")
# MASL.MA(pos, delta, BoxSize, MAS, verbose=verbose)
# delta /= np.mean(delta, dtype=np.float64);
# delta -= 1.0

# # compute the correlation function
# print("Computing correlation function")
# axis = 0
# threads=12
# CF     = PKL.Xi(delta, BoxSize, MAS, axis, threads)


# # get the attributes
# print("Getting attributes")
# r      = CF.r3D      #radii in Mpc/h
# xi0    = CF.xi[:,0]  #correlation function (monopole)
# xi2    = CF.xi[:,1]  #correlation function (quadrupole)
# xi4    = CF.xi[:,2]  #correlation function (hexadecapole)
# Nmodes = CF.Nmodes3D #number of modes
# print("\n\nDone!!\n")

# # log-log plot
# print("Plotting 2pCF")
# plt.figure(figsize=(9,5.5), dpi=100)
# plt.loglog(r, xi0, color='blue')
# plt.grid(alpha=0.5)

# plt.title('TNG300-3 DM Particles ($z=0$) - 2pCF')
# plt.xlabel('$r ~(Mpc~h^{-1})$',  fontsize=16)
# plt.ylabel("$\\xi(r)$", fontsize=16)

# plt.savefig('2pCF_TNG300.jpg', dpi=300, bbox_inches='tight')
# plt.show()
