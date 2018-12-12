import numpy as np
from pyswarm import pso
import subprocess
import sys  
from shark.standard_plots import common
import cosmicDust
import math

###########SWARM PARAMS###########
ss = 5      ## swarmsize         #
mi = 1      ## maximum iterations#
##################################



def callCosmic(x, *args):
	print(x)
	modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, minmass, Omegab, G, rho_crit, sbar, OmegaM, OmegaL, XH=args
	fgas_dissipation=x
	subprocess.call(['shark/build/shark', 'sample.cfg', '-o fgas_dissipation='+str(fgas_dissipation)])
	chi2, chi2_mol=cosmicDust.cosmicDust(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, minmass, Omegab, G, rho_crit, sbar, OmegaM, OmegaL, XH)
	return chi2


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
ub=[1.5]
lb=[0]                        
xopt, fopt=pso(callCosmic, lb, ub, args=args, swarmsize=ss, maxiter=mi) 
print('xopt: ', xopt, 'fopt: ', fopt)
f=open('PSOoutput/results/fgas_dissipation_'+str(ss)+'_'+str(mi)+'.log', 'a+')
f.write('xopt ='+str(xopt)+'\n')
f.write('fopt ='+ str(fopt)+'\n')
f.write('swarmsize ='+str(ss)+'\n')
f.write('maxiter ='+str(mi)+'\n')
f.write('lb ='+str(lb)+'\n')
f.write('ub ='+str(ub)+'\n')
f.close()
