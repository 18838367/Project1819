#bread and butter
import os
import sys  
import subprocess
import multiprocessing
import math
import time
import numpy as np
#edited but based on public code
import pso
#developed internally
import common
import analysis
import constraints
import HPCCommon

###########SWARM PARAMS#############
ss = 16     ## swarmsize           #
mi = 10    ## maximum iterations  #
prc = 0     ## number of processes, set to 0 if using shark submit #
####################################
#############GLOBAL PARAMS########
count=0
def HPCCallPSO(x, *args):
	"""
	- Handler function for running PSO on Shark on a SLURM based cluster. 
	- Swarm size and number of iterations need to be set within the script for now
	- Function needs the relative path to a Shark config file under the -c option
	- For now the subprocess call within must be altered if you are changing shark submit options
	- To find appropriate memory allocations peruse the initial output lines of each particle
	"""
	# this global count just tracks the number of iterations so they can be saved to different files
	global count 
	count = count + 1
	subvols=''
	sT=np.zeros([ss, 3])
	modeldir, outdir, redshift_table, subvolsList, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr, statTest, MFOpt, zOpt, domainUP, domainLW = args
	space=np.genfromtxt('/home/msammons/Project1819/aux/searchSpace.txt', dtype='str')
	
	#this simply corrects an annoying python habit to put things into 1 D arrays if there is only 1 row in file
	if space.size==4 :
		temp=space
		space=np.zeros([1,4], dtype='object_')
		space[0,:]=temp
	
	#This file is read by the shark run/submit? script to determine which values Shark will be run for
	f2=open('/home/msammons/Project1819/aux/particlePositions.ssv', 'w+')
	for i in range(len(x[:,0])):
		for j in range(len(space)):
			if(space[j,1]=='1'):
				val=10**x[i,j]
				f2.write(' -o "'+str(space[j, 0])+'='+str(val)+'"')
			else:
				f2.write(' -o "'+str(space[j, 0])+'='+str(x[i, j])+'"')
			
		f2.write('\n')
	f2.close()
	#calls the submit process and hands in the correct arguments, due to the use of the common.parse_args they can't simply be moved to command line arguments
	#perhaps they could be put in like a config file and read in. Anyway that may be a problem to solve another day
	for i in range(len(subvolsList)):
		subvols=str(subvols)+' '+str(subvolsList[i])
			

	subprocess.call(['./shark-submit', '-a', 'Pawsey0119', '-S', '../build/shark', '-w', '30:00', '-m', '1500M', '-c', '1', '-n', 'PSOSMF'+str(count), '-O', '/mnt/su3ctm/mawson/sharkOut/PSOoutput/PSOSMF'+str(count), '-E', '/home/msammons/Project1819/aux/particlePositions.ssv', '-V', subvols, '../sample.cfg'])
	time.sleep(10)	
	if(len(subvolsList)!=1):
		subvols="multiple_batches"
	else:
		subvols=str(subvolsList[0])
	#above will submit the shark instances for each PSO particle
	#now need to ping queue until the jobs are done before returning the values 
	while(HPCCommon.count_jobs('msammons')>0):
		time.sleep(10)
	

	#below performs the statistical computation of the "goodness" of fit

	for i in range(ss):
		#this directory will eventually need to be made generic, even more config file options is probably the best way forward
		modeldir='/mnt/su3ctm/mawson/sharkOut/PSOoutput/PSOSMF'+str(count)+'/'+str(i)+'/mini-SURFS/my_model'
		for j in range(len(MFOpt)):
			xObs, xMod, yObs, yMod, ydn, yup=constraints.massFunction(modeldir, outdir, redshift_table, subvols, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp, dssfr, ssfrbins, xssfr, MFOpt[j], zOpt[j])	
			sT[i, j] = statTest(xObs, xMod, yObs, yMod, ydn, yup, domainLW[j], domainUP[j]) 
			
	ssT=np.sum(sT,1)
	np.save('/home/msammons/Project1819/aux/tracks/track'+str(count), ssT)
	np.save('/home/msammons/Project1819/aux/tracks/trackP'+str(count), x)
	return ssT


if __name__ == '__main__':

	modeldir, outdir, redshift_table, subvols, obsdir = common.parse_args()
	#some parameters that the use needs to set
	#################################
	MFOpt=['HIMF', 'SMF', 'SMF']
	zOpt=[0,0,1]
	domainUP=[13,13,13]
	domainLW=[7,8,8]
	statTest='StudentT'
	#statTest='Chi2'
	subvolsList=[0]
	#################################

	if statTest == 'StudentT':
		statTest = analysis.nonEqualStudentT
	elif statTest == 'Chi2':
		statTest = analysis.nonEqualChi2
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

	args=(modeldir, outdir, redshift_table, subvolsList, obsdir, GyrToYr, Zsun, XH, MpcToKpc, mlow, mupp, dm, mbins, xmf, imf, mlow2, mupp2, dm2, mbins2, xmf2, ssfrlow, ssfrupp,
		 	dssfr, ssfrbins, xssfr, statTest, MFOpt, zOpt, domainUP, domainLW)
	space=np.genfromtxt('/home/msammons/Project1819/aux/searchSpace.txt')
	if(space.size==4):
		ub=np.zeros(1)
		lb=np.zeros(1)
		ub[0]=space[3]
		lb[0]=space[2]
	else:
		ub=space[:,3]
		lb=space[:,2]
	tStart=time.time()
	os.chdir('../shark/hpc')
	xopt, fopt=pso.pso(HPCCallPSO, lb, ub, args=args, swarmsize=ss, maxiter=mi, processes=prc) 
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
