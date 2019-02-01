#Runs shark from python, just a test to make sure I knew how
# exists within the code directory at the highest level
#imports
import numpy as np
import collections
import subprocess
import sys
import h5py
import common
import utilities_statistics as us
import smf
import matplotlib.pyplot as plt
import analysis
import collections
#subprocess.call(["/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/build/shark", sys.argv[1]])

def massFunction(*args):

	observation = collections.namedtuple('observation', 'label x y yerrup yerrdn err_absolute')
	modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr, MFOpt, zOpt=args
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
	print('stellarMF model dir : ', modeldir)	

	if subvols == 'multiple_batches':
		for index, snapshot in enumerate(redshift_table[zlist]):
		    hdf5_data = common.read_dataMB(modeldir, snapshot, fields, subvols)
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
	else:
		
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
	#load observations 
	if(MFOpt=='SMF'):
		if(zOpt==0):	
			z0obs = []
			lm, p, dpdn, dpup = np.loadtxt('/home/msammons/shark/data/mf/SMF/GAMAII_BBD_GSMFs.dat', usecols=[0,1,2,3], unpack=True)
			xobs = lm
			indx = np.where(p > 0)
			yobs = np.log10(p[indx])
			ytemp=p[indx]-dpdn[indx]
			temp=np.less(ytemp, 0)   		#
			fixed=0.0001*temp+ytemp*np.invert(temp)    #fixed a problem where there were undefined values due to log of negative values, negative values were given a minimum of 0.0001
			print('fixed', fixed)
			ydn = yobs - np.log10(fixed)	 	#
			yup = np.log10(p[indx]+dpup[indx]) - yobs
			z0obs.append((observation("Wright+2017", xobs[indx], yobs, ydn, yup, err_absolute=False), 'o'))
		
		######################################
		#nabbed from plotting part 
			y = hist_smf[0,:]
			ind = np.where(y < 0.)
			xMod=xmf[ind]
			yMod=y[ind]
			xObs=xobs[indx]
			yObs=yobs
	###############>>>>>>>
	
		elif(zOpt==1):
		
		######################################
		#nabbed from plotting part 
			y = hist_smf[2,:]
			ind = np.where(y < 0.)
			#store these ones for chi2
			xMod=xmf[ind]
			yMod=y[ind]
		
			# Moustakas (Chabrier IMF), ['Moustakas+2013, several redshifts']
			zdnM13, lmM13, pM13, dp_dn_M13, dp_up_M13 = np.loadtxt('/home/msammons/shark/data/mf/SMF/SMF_Moustakas2013.dat', usecols=[0,3,5,6,7], unpack=True)
			xobsM13 = lmM13 	
			yobsM13 = np.full(xobsM13.shape, -999.)
			lerrM13 = np.full(xobsM13.shape, -999.)
			herrM13 = np.full(xobsM13.shape, 999.)
			indx = np.where( pM13 < 1)
			yobsM13[indx] = (pM13[indx])
			indx = np.where( dp_dn_M13 > 0)
			lerrM13[indx]  = dp_dn_M13[indx] 
			indx = np.where( dp_up_M13 > 0)
			herrM13[indx]  = dp_up_M13[indx]
		
		
			# Muzzin (Kroupa IMF), ['Moustakas+2013, several redshifts']
			zdnMu13,zupMu13,lmMu13,pMu13,dp_dn_Mu13,dp_up_Mu13 = np.loadtxt('/home/msammons/shark/data/mf/SMF/SMF_Muzzin2013.dat', usecols=[0,1,2,4,5,5], unpack=True)
			# -0.09 corresponds to the IMF correction
			xobsMu13 = lmMu13 - 0.09
			yobsMu13 = np.full(xobsMu13.shape, -999.)
			lerrMu13 = np.full(xobsMu13.shape, -999.)
			herrMu13 = np.full(xobsMu13.shape, 999.)
			indx = np.where( pMu13 < 1)
			yobsMu13[indx] = (pMu13[indx])
			indx = np.where( dp_dn_Mu13 > 0)
			lerrMu13[indx]  = dp_dn_Mu13[indx] 
			indx = np.where( dp_up_Mu13 > 0)
			herrMu13[indx]  = dp_up_Mu13[indx]
			
			
			
			# Wright et al. (2018, several reshifts). Assumes Chabrier IMF.
			zD17, lmD17, pD17, dp_dn_D17, dp_up_D17 = np.loadtxt('/home/msammons/shark/data/mf/SMF/Wright18_CombinedSMF.dat', usecols=[0,1,2,3,4], unpack=True)
			hobs = 0.7
			pD17 = pD17 - 3.0*np.log10(hobs) 
			lmD17= lmD17 - np.log10(hobs)
		
			
			# z1 obs
		#	in_redshift = np.where(zdnM13 == 0.8)
		#	in_redshift = np.where(zdnMu13 == 1)
			in_redshift = np.where(zD17 == 1)
			xObs=lmD17[in_redshift]
			yObs=pD17[in_redshift]
			ydn=dp_dn_D17[in_redshift]
			yup=dp_up_D17[in_redshift]
		else:
			print('zOpt was not a supported option')


##############>>>>>>>>
	elif(MFOpt=='HIMF'):
		if(zOpt==0):
			###################################
			######################################
			#nabbed from plotting part 
			#load observations #HIPASS
			lmHI, pHI, dpHIdn, dpHIup = np.loadtxt('/home/msammons/shark/data/mf/GasMF/HIMF_Zwaan2005.dat', usecols=[0,1,2,3], unpack=True)
			
			#correct data for their choice of cosmology
			hobs = 0.75
			xobs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
			yobs = pHI + np.log10(pow(h0,3)/pow(hobs,3))
			xObsZwaan=xobs
			yObsZwaan=yobs
			ydnZwaan=dpHIdn
			yupZwaan=dpHIup
			#ALFALFA.40
			lmHI, pHI, pdnHI, pduHI = np.loadtxt('/home/msammons/shark/data/mf/GasMF/HIMF_Jones18.dat', usecols=[0,1,2,3], unpack=True)
			
			#correct data for their choice of cosmology
			dpdnHI = pHI - pdnHI
			dpupHI = pduHI - pHI
			hobs = 0.7
			xobs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
			yobs = pHI + np.log10(pow(h0,3)/pow(hobs,3))
			xObsJones=xobs
			yObsJones=yobs
			ydnJones=dpdnHI
			yupJones=dpupHI
			
			# Predicted HIMF
			y = hist_HImf[0,:]
			ind = np.where(y < 0.)
			xMod=xmf[ind]
			yMod=y[ind]
			xObs=xObsJones
			yObs=yObsJones
			ydn=ydnJones
			yup=yupJones
		else:
			print('zOpt was not a supported option')
	else:
		print('MFOpt was not a supported option')


	print('xMod & yMod :', xMod, yMod)
	print('xObs & yObs :', xObs, yObs)
	#print('studentT :', studentT)
	return xObs, xMod, yObs, yMod, ydn, yup


if __name__=='__main__':
	
	observation = collections.namedtuple('observation', 'label x y yerrup yerrdn err_absolute')
	
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

	stellarMF(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr)
