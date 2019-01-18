import subprocess
def caller(x, *args)
	beta, eta, rho = x
	xopt, fopt=subprocess.call(['python', '/Users/mawsonsammons/Documents/ICRARInternship/Project/code/PSOtests/myTestV2.py', beta, eta, rho])
	return xopt, fopt
