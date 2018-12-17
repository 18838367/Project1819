import numpy as np
from PSOcode.pyswarm.pyswarm import pso
import subprocess
import sys  
from shark.standard_plots import common
import stellarMF
import math
import time
import multiprocessing
###########SWARM PARAMS#############
ss = 28      ## swarmsize           #
mi = 10      ## maximum iterations  #
prc = 4     ## number of processes # 
####################################
#############GLOBAL PARAMS########
count=0

def callSMF(x, *args):
	global count
	count = count + 1
	print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX', count, 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
	modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr = args
	fgas_dissipation, cgal=x
	modeldir='/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/output/output'+str(str(count)+'_'+str(multiprocessing.current_process()._identity))+'/mini-SURFS/my_model'
	print('--------------------------------------', x, '-----------------------------------------')
	print('#######################################', multiprocessing.current_process()._identity, '#########################')
	subprocess.call(['shark/build/shark', 'sample.cfg', '-o galaxy_mergers.fgas_dissipation='+str(fgas_dissipation), '-o galaxy_mergers.cgal='+str(cgal), '-o execution.output_directory=/Users/mawsonsammons/Documents/ICRARInternship/Project/code/shark/output/output'+str(count)+'_'+str(multiprocessing.current_process()._identity)])
	chi2 = stellarMF.stellarMF(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2										, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr)
	return chi2


if __name__ == '__main__':
	modeldir, outdir, redshift_table, subvols, obsdir = common.parse_args()
	print(modeldir)
	print(outdir)
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

	args=(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp,
		 	dssfr, ssfrbins, xssfr)
	ub=[1.5, 0.5]
	lb=[0, 0.45]
	tStart=time.time()
	xopt, fopt=pso.pso(callSMF, lb, ub, args=args, swarmsize=ss, maxiter=mi, processes=prc) 
	tEnd=time.time()
	print('xopt: ', xopt, 'fopt: ', fopt)
	f=open('PSOoutput/results/SMFgm57_'+str(ss)+'_'+str(mi)+'.log', 'a+')
	f.write('Time for PSO ='+str(tEnd-tStart)+'\n')
	f.write('Number of processes = '+str(prc)+'\n')
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
