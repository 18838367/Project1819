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
# outdir, obsdir, h0, plotz, hist_smf, hist_smf_cen, hist_smf_sat, hist_smf_er     r, hist_smf_30kpc)
redshift_table=np.loadtxt("/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/input/redshifts.txt", usecols=[1], unpack=True)
zlist=(0, 0.5, 1, 2, 3, 4)
fields = {'galaxies': ('sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge',
                            'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                            'matom_bulge', 'mmol_bulge', 'mgas_bulge',
                            'mgas_metals_disk', 'mgas_metals_bulge',
                            'mstars_metals_disk', 'mstars_metals_bulge', 'type',
                'mvir_hosthalo', 'rstar_bulge')}
 
for index, snapshot in enumerate(redshift_table[zlist]):
    hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
    mass = prepare_data(hdf5_data, index, hist_smf, hist_smf_err, hist_smf_cen,
                         hist_smf_sat, hist_smf_30kpc, hist_HImf, hist_HImf_cen, hist_HImf_sat,
                         hist_H2mf, hist_H2mf_cen, hist_H2mf_sat, mainseq, mainseqsf,
                         sfe, mainseq_cen, mainseqsf_cen, sfe_cen, mainseq_sat,
                         mainseqsf_sat, sfe_sat, mzr, fmzr, mzr_cen, mzr_sat, plotz,
                         plotz_HImf, passive_fractions, hist_ssfr, mszr, mszr_cen,
             mszr_sat, mainseqsf_1s, mainseqHI, mainseqH2)

    h0 = hdf5_data[0]
    if index == 0:
        (sfr_disk, sfr_burst, mdisk, mbulge) = hdf5_data[2:6]
        sfr_seq = np.zeros(shape = (2, len(mdisk)))
        ind  = np.where((sfr_disk + sfr_burst > 0) & (mdisk + mbulge > 0))
        sfr_seq[0,ind] = mass[ind]
        sfr_seq[1,ind] = np.log10((sfr_disk[ind] + sfr_burst[ind]) / h0 / GyrToYr)
print('look i did it') 

