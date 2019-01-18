import numpy as np
import matplotlib.pyplot as plt
import sys

num_particles=sys.argv[0]
dimension=sys.argv[1]
iterations=sys.argv[2]
positions=np.zeros([num_particles, dimensions, iterations])
for i in range(iterations):
	positions[:, :, i]=np.load('/home/msammons/aux/tracks/track'+str(i+1))
np.save('positions', positions)	
