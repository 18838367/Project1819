#this analysis file contains the python methods used to quantify the "goodness" of fit
#between model and observed parameters, to run just call _____ and give it the model 
# data and obs data, (x and y values of course)

import numpy as np 

#linear interpolator
def linIntp(x1, x2, y1, y2, xR):
	x1=np.asarray(x1)
	x2=np.asarray(x2)
	y1=np.asarray(y1)
	y2=np.asarray(y2)
	xR=np.asarray(xR)
	grad=(y2-y1)/(x2-x1)
	yR=(xR-x1)*grad+y1
	return yR

