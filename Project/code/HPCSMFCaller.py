import numpy as np
import psoHPC
import subprocess
import sys  
import common
import stellarMF
import math
import time
import HPCCommon
import os
import multiprocessing
###########SWARM PARAMS#############
ss = 24      ## swarmsize           #
mi = 10      ## maximum iterations  #
prc = 1     ## number of processes # 
####################################
#############GLOBAL PARAMS########
count=0

def HPCCallSMF(x, *args):
	global count 
	count = count + 1
	chi2=np.zeros(ss)
	modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr = args
	names=np.genfromtxt('/home/msammons/aux/fields.txt', dtype='str')
	if names.size==1 :
		temp=names
		names=np.zeros(1, dtype='object_')
		names[0]=temp
	looper=0
	f2=open('/home/msammons/aux/particlePositions.ssv', 'w+')
	for i in range(len(x[:,0])):
		for j in range(len(names)):
			f2.write(str(str(names[j])+'='+str(x[i, j]))+'\n')
			looper=looper+1
	f2.close()
	subprocess.call(['./specialistSharkSubmit', '-a', 'Pawsey0119', '-S', '../build/shark', '-w', '4:00', '-m', '1000M', '-c', '1', '-n', 'PSOSMF'+str(count), '-O', '/mnt/su3ctm/mawson/sharkOut/PSOoutput/PSOSMF'+str(count), '-N', str(len(names)), '-V', '1', '../sample.cfg'])
	time.sleep(10)	
	#above will submit the shark instances for each PSO particle
	#now need to ping queue until the jobs are done before returning the values 
	while(HPCCommon.count_jobs('msammons')>0):
		time.sleep(10)
	
	#####
	#one ping, one ping only
	######
	runNum=np.genfromtxt('/home/msammons/aux/SOD.txt', dtype='str')
	for i in range(ss):
		modeldir=str(runNum)+'/'+str(i)+'/mini-SURFS/my_model'
		chi2[i]=stellarMF.stellarMF(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr)
	return chi2


if __name__ == '__main__':
	modeldir, outdir, redshift_table, subvols, obsdir = common.parse_args()
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
	bounds=np.genfromtxt('/home/msammons/aux/bounds.txt')
	if(len(bounds.shape)==1):
		ub=np.zeros(1)
		lb=np.zeros(1)
		ub[0]=bounds[0]
		lb[0]=bounds[1]
	else:
		ub=bounds[0,:]
		lb=bounds[1,:]
	tStart=time.time()
	os.chdir('shark/hpc')
	xopt, fopt=psoHPC.pso(HPCCallSMF, lb, ub, args=args, swarmsize=ss, maxiter=mi, processes=prc) 
	tEnd=time.time()
	print('xopt: ', xopt, 'fopt: ', fopt)
	f=open('/mnt/su3ctm/mawson/sharkOut/PSOoutput/results/SMFgm5_'+str(ss)+'_'+str(mi)+'.log', 'a+')
	f.write('Time for PSO ='+str(tEnd-tStart)+'\n')
	f.write('Number of processes = '+str(prc)+'\n')
	f.write('Count ='+str(count)+'\n')
	f.write('Searched parameters = fgas_dissipation \n')
	f.write('xopt ='+str(xopt)+'\n')
	f.write('fopt ='+ str(fopt)+'\n')
	f.write('swarmsize ='+str(ss)+'\n')
	f.write('maxiter ='+str(mi)+'\n')
	f.write('lb ='+str(lb)+'\n')
	f.write('ub ='+str(ub)+'\n')
	f.write('-----------------------------------------------------\n')
	f.close()
