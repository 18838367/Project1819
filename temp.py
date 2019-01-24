import collections
import functools
import math

import numpy as np

import common
import utilities_statistics as us


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

    (h0, volh, sfr_disk, sfr_burst, mdisk, mbulge, rstar_disk, mBH, mHI, mH2, 
     mgas_disk, mHI_bulge, mH2_bulge, mgas_bulge, mgas_metals_disk, mgas_metals_bulge, 
     mstars_metals_disk, mstars_metals_bulge, typeg, mvir_hosthalo, rstar_bulge) = hdf5_data

h0=galaxies['cosmology']['h'][()]
volh=galaxies['run_infoo']['']
#establishing data - model
h0=galaxies["cosmology"]['h'][()]
mHI=galaxies['global']['m_hi'][()]
mH2=galaxies['global']['m_h2'][()]
mColdMetals=galaxies['global']['mcold_metals'][()]
volh=galaxies["run_info"]["effective_volume"][()]
sfr_disk=galaxies['galaxies']['sfr_disk'][()]
sfr_burst=galaxies['galaxies']['sfr_burst'][()]
mdisk=galaxies['']
redshifts=np.loadtxt("/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/input/redshifts.txt", usecols=[1], unpack=True)
    mgas = mgas_disk+mgas_bulge
    mgas_metals = mgas_metals_disk+mgas_metals_bulge

    mass          = np.zeros(shape = len(mdisk))
    mass_30kpc    = np.zeros(shape = len(mdisk))
    massd_30kpc   = np.zeros(shape = len(mdisk))
    massb_30kpc   = np.zeros(shape = len(mdisk))
    mass_atom     = np.zeros(shape = len(mdisk))
    mass_atom_cen = np.zeros(shape = len(mdisk))
    mass_atom_sat = np.zeros(shape = len(mdisk))

    mass_mol = np.zeros(shape = len(mdisk))
    mass_mol_cen = np.zeros(shape = len(mdisk))
    mass_mol_sat = np.zeros(shape = len(mdisk))

    ind = np.where((mdisk+mbulge) > 0.0)
    mass[ind] = np.log10(mdisk[ind] + mbulge[ind]) - np.log10(float(h0))
    print('number of galaxies with mstars>0 and max mass: %d, %d' % (len(mass[ind]), max(mass[ind])))

    H, _ = np.histogram(mass,bins=np.append(mbins,mupp))
    hist_smf[index,:] = hist_smf[index,:] + H
    ran_err = np.random.normal(0.0, 0.25, len(mass))
    mass_err = mass + ran_err
    H, _ = np.histogram(mass_err,bins=np.append(mbins,mupp))
    hist_smf_err[index,:] = hist_smf_err[index,:] + H

    #Calculate the stellar mass contained in 30pkpc, assuming an exponential profile for the disk and a Plummer profile for the bulge.
    ind = np.where((mdisk > 0.0)  & (rstar_disk > 0))
    massd_30kpc[ind] = mdisk[ind] * (1.0 - (1.0 + 30.0/(rstar_disk[ind]/1.67/h0 * MpcToKpc)) * np.exp(-30.0/(rstar_disk[ind]/1.67/h0 * MpcToKpc)))
    ind = np.where((mbulge > 0.0)  & (rstar_bulge > 0))
    massb_30kpc[ind] = mbulge[ind] * pow(30.0, 3.0) / pow((pow(30.0, 2.0) + pow(rstar_bulge[ind]/1.3/h0 * MpcToKpc, 2.0)), 3.0/2.0)

    ind = np.where((massd_30kpc + massb_30kpc) > 0)
    mass_30kpc[ind] = np.log10(massd_30kpc[ind] + massb_30kpc[ind]) - np.log10(float(h0))
    H, _ = np.histogram(mass_30kpc,bins=np.append(mbins,mupp))
    hist_smf_30kpc[index,:] = hist_smf_30kpc[index,:] + H


    ind = np.where(typeg == 0)
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_cen[index,:] = hist_smf_cen[index,:] + H
    ind = np.where(typeg > 0)
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_sat[index,:] = hist_smf_sat[index,:] + H



    if volh > 0:
        vol = volh/pow(h0,3.)  # In Mpc^3
        hist_smf[index,:]  = hist_smf[index,:]/vol/dm
        hist_smf_30kpc[index,:]= hist_smf_30kpc[index,:]/vol/dm
        hist_smf_err[index,:]  = hist_smf_err[index,:]/vol/dm
        hist_smf_cen[index,:]  = hist_smf_cen[index,:]/vol/dm
        hist_smf_sat[index,:]  = hist_smf_sat[index,:]/vol/dm

        hist_HImf[index,:] = hist_HImf[index,:]/vol/dm
        hist_HImf_cen[index,:] = hist_HImf_cen[index,:]/vol/dm
        hist_HImf_sat[index,:] = hist_HImf_sat[index,:]/vol/dm

        hist_H2mf[index,:] = hist_H2mf[index,:]/vol/dm
        hist_H2mf_cen[index,:] = hist_H2mf_cen[index,:]/vol/dm
        hist_H2mf_sat[index,:] = hist_H2mf_sat[index,:]/vol/dm

        plotz[index]     = True
        plotz_HImf[index]= True
    else:
        plotz[index]     = False
        plotz_HImf[index]= False

