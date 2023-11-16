from scipy import optimize
import numpy as np
import sys

#
# Reference URL: https://www.monte-carlo-note.com/2017/04/How-to-solve-equation-with-python.html
#

args = sys.argv
dirname = args[1]
DGs = []
fname_refs = dirname +  '/state02.dat'
fname_soln = dirname +  '/state01.dat'

w_refs = np.loadtxt(fname_refs)[:,0]
w_soln = -1.0E0*np.loadtxt(fname_soln)[:,1]

def f_refs(w,x):
    return 1.0E0/(1.0E0+np.exp(beta*(w+x)))

def f_soln(w,x):
    return 1.0E0/(1.0E0+np.exp(-beta*(w+x)))

def sum_f_refs(w,x):
    return np.sum(f_refs(w,x))

def sum_f_soln(w,x):
    return np.sum(f_soln(w,x))

def diff(x):
    return sum_f_refs(w_refs,x) - sum_f_soln(w_soln,x)


T = 300.0 #K
kBT = 2.479 * T/298.0 #kJ/mol
beta = 1.0E0/kBT

# BAR
D = optimize.fsolve(diff,0)

NA = len(w_refs)
NB = len(w_soln) # assume the same number
print(NA,NB)
print(D)
print(-kBT * np.log(NA/NB))
DeltaG_bar = -kBT * np.log(NB/NA) - D

outf = dirname +  '/DeltaG_BAR.dat'
with open(outf, 'wt') as f:
    f.write("kJ/mol\tkcal/mol\n")
    l = '{0:6.3f}\t{1:6.3f}\n'.format(DeltaG_bar[0],DeltaG_bar[0]*0.2390060)
    f.write(l)
