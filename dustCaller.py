import numpy as np
from PSOcode.pyswarm.pyswarm import pso
import subprocess
import sys  
from shark.standard_plots import common
import cosmicDust
import math
import time
import multiprocessing 

###########SWARM PARAMS#############
ss = 4      ## swarmsize           #
mi = 4      ## maximum iterations  #
prc = 4		## number of processes #
####################################
#############GLOBAL PARAMS########
count=0

def callCosmic(x, *args):
	global count
	count = count + 1
	modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, minmass, Omegab, G, rho_crit, sbar, OmegaM, OmegaL, XH=args
	fgas_dissipation, cgal=x
	subprocess.call(['shark/build/shark', 'sample.cfg', '-o galaxy_mergers.fgas_dissipation='+str(fgas_dissipation), '-o galaxy_mergers.cgal='+str(cgal)])
	chi2, chi2_mol=cosmicDust.cosmicDust(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, minmass, Omegab, G, rho_crit, sbar, OmegaM, OmegaL, XH)
	return chi2


if __name__ == '__main__':
	modeldir, outdir, redshift_table, subvols, obsdir = common.parse_args()
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
	
	args=(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, minmass, Omegab, G, rho_crit, sbar, OmegaM, OmegaL, XH)
	ub=[1.5, 0.5]
	lb=[0, 0.45]
	tStart=time.time()
	xopt, fopt=pso.pso(callCosmic, lb, ub, args=args, swarmsize=ss, maxiter=mi, processes=prc) 
	tEnd=time.time()
	print('xopt: ', xopt, 'fopt: ', fopt)
	f=open('PSOoutput/results/Dustgm57_'+str(ss)+'_'+str(mi)+'.log', 'a+')
	f.write('Time for PSO ='+str(tEnd-tStart)+'\n')
	f.write('number of processes = '+str(prc)+'\n')
	f.write('Count ='+str(count)+'\n')
	f.write('Searched parameters = fgas_dissipation, cgal \n')
	f.write('xopt ='+str(xopt)+'\n')
	f.write('fopt ='+ str(fopt)+'\n')
	f.write('swarmsize ='+str(ss)+'\n')
	f.write('maxiter ='+str(mi)+'\n')
	f.write('lb ='+str(lb)+'\n')
	f.write('ub ='+str(ub)+'\n')
	f.write('-----------------------------------------------------\n')
	f.close()
