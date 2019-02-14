"""
This script produces the diagnostic plots used to help visualise the PSO process
One is a 3D representation of the swarm movement in the parameter space the other
is the Log Liklihood over iteration for each particle
"""


import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D

ss=int(sys.argv[1])
it=int(sys.argv[2])

sT=np.zeros([ss, it])
for i in range(it):
    sT[:,i]=np.load('/home/msammons/Project1819/aux/tracks/track'+str(i+1)+'.npy')
pos=np.zeros([ss,3,it])
for i in range(it):
    pos[:,:,i]=np.load('/home/msammons/Project1819/aux/tracks/trackP'+str(i+1)+'.npy')


##Performance diagnostic plot
fig = plt.figure()
ax = fig.add_subplot(111)
plt.ylabel('Log Likelihood')
plt.xlabel('Number of Iterations')
for i in range(ss):
    ax.plot(np.arange(it), sT[i,:], label='Particle '+str(i))
fig.savefig('performance.pdf')


##3D map diagnostic plot
#configure colours to show sT
cm = plt.get_cmap('plasma')
cNorm = matplotlib.colors.Normalize(vmin=np.amin(sT), vmax=np.amax(sT))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
fig = plt.figure()
ax = Axes3D(fig)
for i in range(ss):
    ax.scatter(pos[i,0,:], pos[i,1,:], pos[i,2,:], c=scalarMap.to_rgba(sT[i,:]), s=(np.flip(np.arange(20), 0)/10)**9)
    #ax.plot(pos[i, 0, :], pos[i,1,:], pos[i, 2, :], linewidth=1, color='grey')
scalarMap.set_array(sT)
cbar=fig.colorbar(scalarMap)
cbar.set_label('Log Likelihood', fontsize=16)
ax.set_xlabel('\n$\log_{10}(\\tau_{reinc})$', fontsize=16)
ax.set_ylabel('\n$\log_{10}(M_{norm})$', fontsize=16)
ax.set_zlabel('\n$\\gamma$', fontsize=16)
ax.set_zticks([0, -1, -2])
ax.set_xticks([0, 0.4, 0.8, 1.2])
ax.set_yticks([9, 10, 11])
fig.savefig('3DPSOC.pdf')

