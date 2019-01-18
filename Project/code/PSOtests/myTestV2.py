import numpy as np 
from pyswarm import pso


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


beta=0.7
eta=-1
rho=0.19
args= (beta, eta, rho)

ub = [100000, 100000, 10000]
lb = [0, 0, 0]
xopt, fopt = pso(myFunc, lb, ub, f_ieqcons=condition, args=args)
print(xopt, fopt)
