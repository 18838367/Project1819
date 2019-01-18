#this analysis file contains the python methods used to quantify the "goodness" of fit
#between model and observed parameters, to run just call _____ and give it the model 
# data and obs data, (x and y values of course)

import numpy as np 
import scipy.special
import math

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
	closest1=kept[placed1]
	placed2=np.add(placed1,1)
	#below deals with the fringe case; when there is no model x value greater than
	#a specific observation x value 
	temp=np.where(placed2 > (len(kept)-1))
	placed2[temp]=placed2[temp]-1
	closest2=kept[placed]
	#print("-----------------")
	#print(closest1, closest2)
	return placed1, placed2, closest1

def myChi2(Obs, Mod, err):
	chi2=(Mod-Obs)**2/(err**2)
	return chi2

def myStudentT(Obs, Mod, err):
	var=sum((Obs-Mod)**2)/len(Obs)
	nu=(2*var)/(var-1)
	print('nu :', nu)
	x=(Mod-Obs)**2/err
	print('x :', x)
	sT=(scipy.special.gamma((nu+1)/2.0))/((nu*math.pi)**(0.5)*scipy.special.gamma(nu/2.0))*(1+x/nu)**(-1*(nu+1)/2.0)
	print('st :', sT)
	return sT	

def nonEqualChi2(xObs, xMod, yObs, yMod, ydn, yup, above=666, below=666):
#invokes above methods to get a chi2 value from two sets of data with different x
#values
# CAUTION if models x values do not sufficiently resolve curve then this will 
#be dead wrong
	if above==666 and below==666:
		above=np.amin(xObs)-1.0
		below=np.amax(xObs)+1.0
	x1, x2, closest1= nearestNeighbours(xObs, xMod)
	yR = linIntp(xMod[x1], xMod[x2], yMod[x1], yMod[x2], xObs)
	tempC=np.greater(yR,yObs)
	err=yup*tempC+ydn*np.invert(tempC)
	print('yup & ydn :', yup, ydn)
	print('tempC :', tempC)
	print('err :', err)
	chi2 = myChi2(yObs, yR, err)
	tempA=np.greater(xObs,above)
	tempB=np.less(xObs,below)
	chi2=chi2*tempA*tempB
	sumChi2=np.sum(chi2)
	return sumChi2

def nonEqualT(xObs, xMod, yObs, yMod, ydn, yup, above=666, below=666):
#invokes above methods to get a chi2 value from two sets of data with different x
#values
# CAUTION if models x values do not sufficiently resolve curve then this will 
#be dead wrong
        if above==666 and below==666:
                above=np.amin(xObs)-1.0
                below=np.amax(xObs)+1.0
        x1, x2, closest1= nearestNeighbours(xObs, xMod)
        yR = linIntp(xMod[x1], xMod[x2], yMod[x1], yMod[x2], xObs)
        tempC=np.greater(yR,yObs)
        err=yup*tempC+ydn*np.invert(tempC)
        print('yup & ydn :', yup, ydn)
        print('tempC :', tempC)
        print('err :', err)
        sT = myStudentT(yObs, yR, err)
        tempA=np.greater(xObs,above)
        tempB=np.less(xObs,below)
        sT=sT*tempA*tempB+np.invert(tempA*tempB)
	print('sT after: ', sT)
        LL=np.sum(-1*np.log(sT))
	print('LL :', LL)
        return LL
