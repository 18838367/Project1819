#Runs shark from python, just a test to make sure I knew how
# exists within the code directory at the highest level
#imports

import subprocess
import sys
import h5py
import hdf5Common
#first system argument should be either an option or the path to the config file
#if just want to read data without re-running shark just put -V or -h 
#this will just give you the version of shark or the help but wont actually run
subprocess.call(["/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/build/shark", sys.argv[1]])

#second and third system arguments are the timestep and subvolume respectively
galaxies=h5py.File("/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/output/mini-SURFS/my_model/"+str(sys.argv[2])+"/"+str(sys.argv[3])+"/galaxies.hdf5", 'r')

print("KeysM: ", hdf5Common.keys(galaxies))
#an example of reading in the data
prompt1=input('input key:   ')
keychainM=hdf5Common.keys(galaxies)
keychainM_1=hdf5Common.keys(galaxies[prompt1])
print("KeysM_1: ", keychainM_1)
#these extra input parameters are for navigating the filesystem of hdf5
prompt2=input('input key:   ')
data=galaxies[prompt1][prompt2][()]
print(data)

