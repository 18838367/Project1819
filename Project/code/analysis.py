#this analysis file contains the python methods used to quantify the "goodness" of fit
#between model and observed parameters, to run just call _____ and give it the model 
# data and obs data, (x and y values of course)

import numpy as np 

#linear interpolator
def linIntp(x1, x2, y1, y2, xR):
#does a linear interpolation of yR at xR from the line constructed between points
#x1,y1 and x2,y2
	#initialising what might be lists to actual arrays
	x1=np.asarray(x1)
	x2=np.asarray(x2)
	y1=np.asarray(y1)
	y2=np.asarray(y2)
	xR=np.asarray(xR)
	#finding exceptions
	temp=np.equal(x1,x2)
	
	#doing interpolation 
	grad=(y2-y1)/(x2-x1)
	yR=(xR-x1)*grad+y1
	#correcting exceptions
	yR[temp]=y1[temp]
	return yR


def nearestNeighbours(xObs, xMod):
#finds the two nearest xMod values for each xObs value, used for finding 
#closest model x values for linear interpolating a model y value at each obs x value
	xObs=np.asarray(xObs)
	xMod=np.asarray(xMod)
	kept=np.copy(xMod)
	LObs=len(xObs)
	LMod=len(xMod)
	xObs=np.expand_dims(xObs, axis=1)
	xMod=np.expand_dims(xMod, axis=1)
	xObs=np.repeat(xObs, LMod, axis=1)
	xMod=np.repeat(xMod, LObs, axis=1)
	xMod=xMod.T
	diffs=xObs-xMod
	#interesting point: the smallest point (the one you are looking for) will be
	#the point just before the first negative value in a row
	#this could be used in an alternative method much to your advantage
	temp=np.greater(diffs,0)
	altered=temp*diffs + np.invert(temp)*(10**30)
	mins=altered.min(1)	
	mins=np.expand_dims(mins, axis=1)
	mins=np.repeat(mins, LMod, axis=1)
	placed=np.equal(mins, diffs)*np.repeat(np.expand_dims(np.arange(0,LMod), axis=1), LObs, axis=1).T
	placed1=np.sum(placed, axis=1)
	#closest1=kept[placed1]
	placed2=np.add(placed1,1)
	#below deals with the fringe case; when there is no model x value greater than
	#a specific observation x value 
	temp=np.where(placed2 > (len(kept)-1))
	placed2[temp]=placed2[temp]-1
	#closest2=kept[placed]
	#print("-----------------")
	#print(closest1, closest2)
	return placed1, placed2

def myChi2(Obs, Mod):
	chi2=(Mod-Obs)**2/Obs
	return chi2

def nonEqualChi2(xObs, xMod, yObs, yMod):
#invokes above methods to get a chi2 value from two sets of data with different x
#values
# CAUTION if models x values do not sufficiently resolve curve then this will 
#be dead wrong
	x1, x2 = nearestNeighbours(xObs, xMod)
	yR = linIntp(xMod[x1], xMod[x2], yMod[x1], yMod[x2], xObs)
	chi2 = myChi2(yObs, yR)
	return np.sum(chi2)
