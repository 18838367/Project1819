import numpy as np 
from pyswarm import pso
import sys

def myFunc(x, *args):
	a, b, c=x
	beta, eta, rho=args
	y=a**b-3*b*a+c
	return y 

def condition(x, *args):
	a, b, c=x
	beta, eta, rho=args
	y=c**beta
	return y


beta=float(sys.argv[1])
eta=float(sys.argv[2])
rho=float(sys.argv[3])
args= (beta, eta, rho)

ub = [100000, 100000, 10000]
lb = [0, 0, 0]
xopt, fopt = pso(myFunc, lb, ub, f_ieqcons=condition, args=args)
print(xopt, fopt)
