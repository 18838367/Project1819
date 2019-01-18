import numpy as np 
from pyswarm import pso


def quadratic(x):
	y=x**2
	return y 
ub = [10]
lb = [0]
xopt, fopt = pso(quadratic, lb, ub)
print(xopt, fopt)
