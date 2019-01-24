import numpy as np
import psoHPC
import subprocess
import sys  
import common
import HIMF
import stellarMF
import stellarMF1
import math
import time
import HPCCommon
import os
import multiprocessing
###########SWARM PARAMS#############
ss = 16     ## swarmsize           #
mi = 10      ## maximum iterations  #
prc = 1     ## number of processes # 
####################################
#############GLOBAL PARAMS########
count=0

def HPCCallSMF(x, *args):
	global count 
	count = count + 1
	sT=np.zeros([ss, 3])
	modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr = args
	names=np.genfromtxt('/home/msammons/aux/fields.txt', dtype='str')
	if names.size==2 :
		temp=names
		names=np.zeros([1,1], dtype='object_')
		names[0,0]=temp
	f2=open('/home/msammons/aux/particlePositions.ssv', 'w+')
	for i in range(len(x[:,0])):
		for j in range(len(names)):
			if(names[j,1]=='1'):
				val=10**x[i,j]
				f2.write(' -o "'+str(names[j, 0])+'='+str(val)+'"')
			else:
				f2.write(' -o "'+str(names[j, 0])+'='+str(x[i, j])+'"')
		f2.write('\n')
	f2.close()
	subprocess.call(['./shark-submit', '-a', 'Pawsey0119', '-S', '../build/shark', '-w', '4:00', '-m', '1500M', '-c', '2', '-n', 'PSOSMF'+str(count), '-O', '/mnt/su3ctm/mawson/sharkOut/PSOoutput/PSOSMF'+str(count), '-E', str(ss), '-V', '0 1', '../sample.cfg'])
	time.sleep(10)	
	subvols="multiple_batches"
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
		sT[i, 0]=HIMF.HIMF(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr)
		sT[i, 1]=stellarMF.stellarMF(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr)
		sT[i, 2]=stellarMF1.stellarMF1(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr)
			
	ssT=np.sum(sT,1)
	np.save('/home/msammons/aux/tracks/track'+str(count), ssT)
	np.save('/home/msammons/aux/tracks/trackP'+str(count), x)
	return ssT


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
	f=open('/mnt/su3ctm/mawson/sharkOut/PSOoutput/results/SMFri123_'+str(ss)+'_'+str(mi)+'.log', 'a+')
	f.write('Time for PSO ='+str(tEnd-tStart)+'\n')
	f.write('Number of processes = '+str(prc)+'\n')
	f.write('Count ='+str(count)+'\n')
	f.write('Searched parameters = tau_reinc mhalo_norm halo_mass_power \n')
	f.write('xopt ='+str(xopt)+'\n')
	f.write('fopt ='+ str(fopt)+'\n')
	f.write('swarmsize ='+str(ss)+'\n')
	f.write('maxiter ='+str(mi)+'\n')
	f.write('lb ='+str(lb)+'\n')
	f.write('ub ='+str(ub)+'\n')
	f.write('-----------------------------------------------------\n')
	f.close()
