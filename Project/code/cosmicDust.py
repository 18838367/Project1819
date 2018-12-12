#Runs shark from python, just a test to make sure I knew how
# exists within the code directory at the highest level
#imports
import numpy as np
#import subprocess
#import sys
import h5py
import math
from shark.standard_plots import common
#import matplotlib.pyplot as plt
from shark.standard_plots import utilities_statistics as us
from shark.standard_plots import global_quantities
import analysis
#first system argument should be either an option or the path to the config file
#if just want to read data without re-running shark just put -V or -h 
#this will just give you the version of shark or the help but wont actually run
#subprocess.call(["/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/build/shark", sys.argv[1]])

#some of this has been adapted from the plotting functions within standar_plots

##################################
# Constants

GyrToYr = 1e9
Zsun = 0.0127
minmass = 1.0
Omegab = 0.0491
G    = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun
rho_crit = 3.0 * pow(100.0,2.0) / 8 / math.pi / G #in units of h^2*Msun/Mpc^3
sbar = rho_crit * Omegab
OmegaM = 0.3121
OmegaL = 0.6879
XH = 0.72

modeldir, outdir, redshift_table, subvols, obsdir = common.parse_args()
x=(modeldir, outdir, redshift_table, subvols, obsdir)

def cosmicDust(*args):
	modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, minmass, Omegab, G, rho_crit, sbar, OmegaM, OmegaL, XH = args
	fields = {'global': ('redshifts', 'm_hi', 'm_h2', 'mcold', 'mcold_metals',
	                     'mhot_halo', 'mejected_halo', 'mstars', 'mstars_bursts_mergers', 'mstars_bursts_diskinstabilities',
	                     'm_bh', 'sfr_quiescent', 'sfr_burst', 'm_dm', 'mcold_halo', 'number_major_mergers',
	                     'number_minor_mergers', 'number_disk_instabilities')}
	
	# Read data from each subvolume at a time and add it up
	# rather than appending it all together
	for idx, subvol in enumerate(subvols):
	    subvol_data = common.read_data(modeldir, redshift_table[0], fields, [subvol])
	    if idx == 0:
	        hdf5_data = subvol_data
	    else:
	        for subvol_datum, hdf5_datum in zip(subvol_data[3:], hdf5_data[3:]):
	            hdf5_datum += subvol_datum
	
	# Also make sure that the total volume takes into account the number of subvolumes read
	hdf5_data[1] = hdf5_data[1] * len(subvols)
	
	h0, redshifts = hdf5_data[0], hdf5_data[2]
	
	(mstar_plot, mcold_plot, mhot_plot, meje_plot,
	 mstar_dm_plot, mcold_dm_plot, mhot_dm_plot, meje_dm_plot, mbar_dm_plot,
	 sfr, sfrd, sfrb, mstarden, mstarbden_mergers, mstarbden_diskins, sfre, sfreH2, mhrat,
	 mHI_plot, mH2_plot, mH2den, mdustden, omegaHI, mdustden_mol, mcoldden, mhotden,
	 mejeden, history_interactions, mDMden) = global_quantities.prepare_data(hdf5_data, redshifts)
	
	
	#establishing data - obs
	redD17d,redD17u,smdD17,err1,err2,err3,err4=np.loadtxt('shark/data/Global/Driver18_dust.dat', usecols=[1,2,3,4,5,6,7] , unpack=True)
	
	#calculating obs results
	
	hobs = 0.7
	xobs = (redD17d+redD17u)/2.0
	yobs = smdD17 + np.log10(hobs/h0)
	
	err = yobs*0. - 999.
	err = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))
	
	#plotting
	
	#fig=plt.figure(figsize=(5,4.5))
	#ax1=fig.add_subplot(111)
	#ax1.errorbar(us.look_back_time(xobs), yobs, yerr=[err, err])
	
	ind = np.where(mdustden > 0)
	xMod= us.look_back_time((redshifts[len(redshifts)-len(mdustden):])[ind])
	yMod= np.log10(mdustden[ind]*pow(h0,2.0))
	
	#ax1.plot(xMod, yMod,'r', label ='Shark all metals')
	
	ind = np.where(mdustden_mol > 0)
	xMod_mol= us.look_back_time((redshifts[len(redshifts)-len(mdustden_mol):])[ind])
	yMod_mol= np.log10(mdustden_mol[ind]*pow(h0,2.0))
	
	#ax1.plot(xMod_mol, yMod_mol,'r', linestyle = 'dashed', label ='Shark metals in molecular gas')
	
	#plt.show()
	
	#calculating chi2 (or something vaguely similar to it)
	chi2 = analysis.nonEqualChi2(xobs, xMod, yobs, yMod)
	chi2_mol = analysis.nonEqualChi2(xobs, xMod_mol, yobs, yMod_mol)
	
	print('chi2: ', chi2)
	print('chi2_mol: ', chi2_mol)
	return chi2, chi2_mol

cosmicDust(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, minmass, Omegab, G, rho_crit, sbar, OmegaM, OmegaL, XH)
