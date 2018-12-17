#Runs shark from python, just a test to make sure I knew how
# exists within the code directory at the highest level
#imports
import numpy as np
import collections
import subprocess
import sys
import h5py
import hdf5Common
from shark.standard_plots import common
import matplotlib.pyplot as plt
from shark.standard_plots import utilities_statistics as us
import analysis
from shark.standard_plots import smf
#subprocess.call(["/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/build/shark", sys.argv[1]])

def HIMF(*args):
	observation = collections.namedtuple('observation', 'label x y yerrup yerrdn err_absolute')
	modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr = args
	
	zlist=(0, 0.5, 1, 2, 3, 4)
	
	mainseq     = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mainseq_cen = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mainseq_sat = np.zeros(shape = (len(zlist), 3, len(xmf)))
	
	mainseqsf     = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mainseqsf_cen = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mainseqsf_sat = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mainseqsf_1s  = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mainseqHI     = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mainseqH2     = np.zeros(shape = (len(zlist), 3, len(xmf)))
	
	mzr         = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mzr_cen     = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mzr_sat     = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mszr        = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mszr_cen    = np.zeros(shape = (len(zlist), 3, len(xmf)))
	mszr_sat    = np.zeros(shape = (len(zlist), 3, len(xmf)))
	
	fmzr        = np.zeros(shape = (len(zlist), 3, len(xmf)))
	
	sfe         = np.zeros(shape = (len(zlist), 3, len(xmf)))
	sfe_cen     = np.zeros(shape = (len(zlist), 3, len(xmf)))
	sfe_sat     = np.zeros(shape = (len(zlist), 3, len(xmf)))
	
	passive_fractions = np.zeros(shape = (len(zlist), 3, len(xmf2)))
	
	# Histograms
	hist_smf       = np.zeros(shape = (len(zlist), len(mbins)))
	hist_smf_30kpc = np.zeros(shape = (len(zlist), len(mbins)))
	hist_smf_err   = np.zeros(shape = (len(zlist), len(mbins)))
	hist_smf_cen   = np.zeros(shape = (len(zlist), len(mbins)))
	hist_smf_sat   = np.zeros(shape = (len(zlist), len(mbins)))
	
	plotz = np.empty(shape=(len(zlist)), dtype=np.bool_)
	hist_HImf = np.zeros(shape = (len(zlist), len(mbins)))
	hist_HImf_cen = np.zeros(shape = (len(zlist), len(mbins)))
	hist_HImf_sat = np.zeros(shape = (len(zlist), len(mbins)))
	plotz_HImf = np.empty(shape=(len(zlist)), dtype=np.bool_)
	hist_H2mf = np.zeros(shape = (len(zlist), len(mbins)))
	hist_H2mf_cen = np.zeros(shape = (len(zlist), len(mbins)))
	hist_H2mf_sat = np.zeros(shape = (len(zlist), len(mbins)))
	
	hist_ssfr = np.zeros(shape = (len(zlist), len(ssfrbins)))
	
	
	fields = {'galaxies': ('sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge',
	                       'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
	                       'matom_bulge', 'mmol_bulge', 'mgas_bulge',
	                       'mgas_metals_disk', 'mgas_metals_bulge',
	                       'mstars_metals_disk', 'mstars_metals_bulge', 'type',
	           'mvir_hosthalo', 'rstar_bulge')}
	
	for index, snapshot in enumerate(redshift_table[zlist]):
	    hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
	    mass = smf.prepare_data(hdf5_data, index, hist_smf, hist_smf_err, hist_smf_cen,
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
	
	#########################
	#take logs
	
	ind = np.where(hist_smf > 0.)
	hist_smf[ind] = np.log10(hist_smf[ind])
	ind = np.where(hist_smf_30kpc > 0.)
	hist_smf_30kpc[ind] = np.log10(hist_smf_30kpc[ind])
	ind = np.where(hist_smf_cen > 0.)
	hist_smf_cen[ind] = np.log10(hist_smf_cen[ind])
	ind = np.where(hist_smf_sat > 0.)
	hist_smf_sat[ind] = np.log10(hist_smf_sat[ind])
	ind = np.where(hist_smf_err > 0.)
	hist_smf_err[ind] = np.log10(hist_smf_err[ind])
	
	ind = np.where(hist_HImf > 0.)
	hist_HImf[ind] = np.log10(hist_HImf[ind])
	ind = np.where(hist_HImf_cen > 0.)
	hist_HImf_cen[ind] = np.log10(hist_HImf_cen[ind])
	ind = np.where(hist_HImf_sat > 0.)
	hist_HImf_sat[ind] = np.log10(hist_HImf_sat[ind])
	
	ind = np.where(hist_H2mf > 0.)
	hist_H2mf[ind] = np.log10(hist_H2mf[ind])
	ind = np.where(hist_H2mf_cen > 0.)
	hist_H2mf_cen[ind] = np.log10(hist_H2mf_cen[ind])
	ind = np.where(hist_H2mf_sat > 0.)
	hist_H2mf_sat[ind] = np.log10(hist_H2mf_sat[ind])
	
	
	
	###################################
	######################################
	#nabbed from plotting part 
#	fig=plt.figure(figsize=(5,4.5))
#	ax=fig.add_subplot(111)
	#load observations #HIPASS
	lmHI, pHI, dpHIdn, dpHIup = common.load_observation(obsdir, 'mf/GasMF/HIMF_Zwaan2005.dat', [0,1,2,3])
	
	#correct data for their choice of cosmology
	hobs = 0.75
	xobs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
	yobs = pHI + np.log10(pow(h0,3)/pow(hobs,3))
#	ax.errorbar(xobs, yobs, yerr=[dpHIdn,dpHIup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Zwaan+2005")
	xObsZwaan=xobs
	yObsZwaan=yobs
	#ALFALFA.40
	lmHI, pHI, pdnHI, pduHI = common.load_observation(obsdir, 'mf/GasMF/HIMF_Jones18.dat', [0,1,2,3])
	
	#correct data for their choice of cosmology
	dpdnHI = pHI - pdnHI
	dpduHI = pduHI - pHI
	hobs = 0.7
	xobs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
	yobs = pHI + np.log10(pow(h0,3)/pow(hobs,3))
#	ax.errorbar(xobs, yobs, yerr=[dpdnHI,dpduHI], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='x',label="Jones+2018")
	xObsJones=xobs
	yObsJones=yobs
	
	# Predicted HIMF
	y = hist_HImf[0,:]
	ind = np.where(y < 0.)
#	ax.plot(xmf[ind],y[ind],'k',  label ='all galaxies')
	xMod=xmf[ind]
	yMod=y[ind]
#	y = hist_HImf_cen[0,:]
#	ind = np.where(y < 0.)
#	ax.plot(xmf[ind],y[ind],'b', linestyle='dotted', label ='centrals')
#	y = hist_HImf_sat[0,:]
#	ind = np.where(y < 0.)
#	ax.plot(xmf[ind],y[ind],'r', linestyle='dashed', label ='satellites')
	
#	pHI_GD14 = common.load_observation(obsdir, 'shark/data/Models/SharkVariations/HIH2MassFunctions_OtherModels.dat', [0])
	
#	ind = np.where(pHI_GD14 < 0)
#	ax.plot(xmf[ind],pHI_GD14[ind],'BurlyWood', linestyle='dashdot', label = 'GD14')
#########	plt.show()	
	
	
	#############################
	#do chi2
	chi2Zwaan=analysis.nonEqualChi2(xObsZwaan, xMod, yObsZwaan, yMod, 7, 13) 
	chi2Jones=analysis.nonEqualChi2(xObsJones, xMod, yObsJones, yMod, 7, 13)
	print('chi2Zwaan :', chi2Zwaan)
	print('chi2Jones :', chi2Jones)
	print('sum :', chi2Zwaan+chi2Jones)
	return chi2Jones+chi2Zwaan
if __name__=='__main__':
	
	##################################
	# Constants
	GyrToYr = 1e9
	Zsun = 0.0127
	XH = 0.72
	MpcToKpc = 1e3
	
	##################################
	# Mass function initialization
	mlow = 5
	mupp = 14
	dm = 0.2
	mbins = np.arange(mlow,mupp,dm)
	xmf = mbins + dm/2.0
	imf   = 'cha'
	
	mlow2 = 5
	mupp2 = 14
	dm2 = 0.3
	mbins2 = np.arange(mlow2,mupp2,dm2)
	xmf2 = mbins2 + dm2/2.0
	
	ssfrlow = -6
	ssfrupp = 4
	dssfr = 0.2
	ssfrbins = np.arange(ssfrlow,ssfrupp,dssfr)
	xssfr    = ssfrbins + dssfr/2.0
	
	#########################################
	#calculation begins
	#some of these initialisations are not necessary for this exact part 
	#however they allow the use of the prepare data from the standard plots smf code
	#so they are useful in that respect 
	modeldir, outdir, redshift_table, subvols, obsdir = common.parse_args()
	HIMF(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr)
