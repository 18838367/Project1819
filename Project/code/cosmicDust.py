#Runs shark from python, just a test to make sure I knew how
# exists within the code directory at the highest level
#imports
import numpy as np
import subprocess
import sys
import h5py
import hdf5Common
import common
import matplotlib.pyplot as plt
import utilities_statistics as us
import analysis
#first system argument should be either an option or the path to the config file
#if just want to read data without re-running shark just put -V or -h 
#this will just give you the version of shark or the help but wont actually run
subprocess.call(["/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/build/shark", sys.argv[1]])
print(sys.argv[1])
#second and third system arguments are the timestep and subvolume respectively
galaxies=h5py.File("/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/output/mini-SURFS/my_model/"+str(sys.argv[2])+"/"+str(sys.argv[3])+"/galaxies.hdf5", 'r')

#print("KeysM: ", hdf5Common.keys(galaxies))
#an example of reading in the data
#prompt1=input('input key:   ')
#keychainM=hdf5Common.keys(galaxies)
#keychainM_1=hdf5Common.keys(galaxies[prompt1])
#print("KeysM_1: ", keychainM_1)
#these extra input parameters are for navigating the filesystem of hdf5
#prompt2=input('input key:   ')
#data=galaxies[prompt1][prompt2][()]
#print(data)

#some of this has been adapted from the plotting functions within standar_plots

#establishing data - obs
redD17d,redD17u,smdD17,err1,err2,err3,err4=np.loadtxt('shark/data/Global/Driver18_dust.dat', usecols=[1,2,3,4,5,6,7] , unpack=True)


#establishing data - model
h0=galaxies["cosmology"]['h'][()]
mHI=galaxies['global']['m_hi'][()]
mH2=galaxies['global']['m_h2'][()]
mColdMetals=galaxies['global']['mcold_metals'][()]
Zsun=0.0127
volh=galaxies["run_info"]["effective_volume"][()]
redshifts=np.loadtxt("/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/input/redshifts.txt", usecols=[1], unpack=True)

#calculating model results

mdustden=      0.006 * mColdMetals/Zsun / volh
mdustden_mol=  0.006 * mColdMetals/Zsun * mH2 / (mH2 + mHI) / volh

#calculating obs results

hobs = 0.7
xobs = (redD17d+redD17u)/2.0
yobs = smdD17 + np.log10(hobs/h0)

err = yobs*0. - 999.
err = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))

#plotting

fig=plt.figure(figsize=(5,4.5))
ax1=fig.add_subplot(111)
ax1.errorbar(us.look_back_time(xobs), yobs, yerr=[err, err])

ind = np.where(mdustden > 0)
xMod= us.look_back_time((redshifts[len(redshifts)-len(mdustden):])[ind])
yMod= np.log10(mdustden[ind]*pow(h0,2.0))

ax1.plot(xMod, yMod,'r', label ='Shark all metals')

ind = np.where(mdustden_mol > 0)
xMod_mol= us.look_back_time((redshifts[len(redshifts)-len(mdustden_mol):])[ind])
yMod_mol= np.log10(mdustden_mol[ind]*pow(h0,2.0))

ax1.plot(xMod_mol, yMod_mol,'r', linestyle = 'dashed', label ='Shark metals in molecular gas')

plt.show()

#calculating chi2 (or something vaguely similar to it)
chi2 = analysis.nonEqualChi2(xobs, xMod, yobs, yMod)
chi2_mol = analysis.nonEqualChi2(xobs, xMod_mol, yobs, yMod_mol)

print('chi2: ', chi2)
print('chi2_mol: ', chi2_mol)
