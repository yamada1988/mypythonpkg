from pytrr import GroTrrReader
import sys
import math
import numpy as np

#
# Calculate rho(k,q,t) and <rho(k,q> from Gromacs trrfile.
# This script requires pytrr. (https://github.com/andersle/pytrr)
# install pytrr using "pip insall pytrr".
#

args = sys.argv
fname = args[1]
T = fname.split('_')[1].split('.')[0]

pi = math.pi

# Mass Parameters
mO = 16.00 * 1.661e-27 #(kg)
mH =  1.08 * 1.661e-27 #(kg)

# Phisical Parameters
N = 3000
tN = 100
L = 4.37864 #(nm) 
#T = 200 #(K)
M = mO + mH + mH #(u)
kB = 1.3801e-23 #(J K^-1)

# Wave number Parameters
k0 = 2.0e0*pi / L # (nm^-1)
dk = 2.0e0*pi / L # (nm^-1)
#vth = math.sqrt(3.0e0*kB*T/M) # (nm/fs)^-1 = (km/s)^-1
vth = 2.0e0*pi
alpha = 2.0
dq = alpha * 2.0e0*pi / vth
q0 = 0.010
kN = 1
qN = 10

K = [1/math.sqrt(3) *np.array([k0+i*dk, k0+i*dk, k0+i*dk]) for i in range(kN)]
Q = [1/math.sqrt(3) * np.array([q0+i*dq, q0+i*dq, q0+i*dq]) for i in range(qN)]

# rho_i(k,q,t,i)
rho_i = [[[[0.0e0 for i in range(N)] for t in range(tN+1)] for q in range(qN)] for k in range(kN)]

# rho(k,q,t)
rho = [[[0.0e0 for t in range(tN+1)] for q in range(qN)] for k in range(kN)]

Rs = []
R = [0 for i in range(N)]
Vs = []
V = [0 for i in range(N)]

icount = -1

with GroTrrReader(fname) as trrfile:
    for frame in trrfile:
        icount += 1
#        print(icount)
        if icount >= tN:
            break
#        print(frame['step'])
        frame_data = trrfile.get_data()
        for i in range(N):
#            print('molnum:', i+1)
            iO  = 3*i 
            iH1 = 3*i + 1
            iH2 = 3*i + 2

# Calculate single H2O molecule's center of mass (R)
            xO  = np.array(frame_data['x'][iO])
            xH1 = np.array(frame_data['x'][iH1])
            xH2 = np.array(frame_data['x'][iH2])
            R[i] = (mO * xO + mH * xH1 + mH * xH2)/ M 
#            print('coordinate:',R[i])

# Calculate single H2O molecule's velocity (V)
            vO  = np.array(frame_data['v'][iO])
            vH1 = np.array(frame_data['v'][iH1])
            vH2 = np.array(frame_data['v'][iH2])
            V[i] = (mO * vO + mH * vH1 + mH * vH2)/ M
#            print('velocity:',V[i])

# Calculate rho(k,q,t)
            for ik,k in enumerate(K):
                for iq,q in enumerate(Q):
                    theta = np.dot(k,R[i]) + np.dot(q,V[i])
                    rho_i[ik][iq][icount][i] = complex(math.cos(theta), -math.sin(theta))
        for ik,k in enumerate(K):
            for iq,q in enumerate(Q):
                rho[ik][iq][icount] = 0.0e0
                for i in range(N): 
                    rho[ik][iq][icount] += rho_i[ik][iq][icount][i]
                rho[ik][iq][icount] /= N
#                print(k,q,rho[ik][iq][icount])

# Calculate <rho(k,q)>
print("##### Ensemble Averaged density ({0:3d}K) ####".format(int(T)))
rho_ens = [[0.0e0 for iq in range(qN)] for ik in range(kN)]
for ik,k in enumerate(K):
    for iq,q in enumerate(Q):
        rho_ens[ik][iq] = 0.0e0
        for i in range(tN):
            rho_ens[ik][iq] += rho[ik][iq][i]
        rho_ens[ik][iq] /= tN
        print(k,q,rho_ens[ik][iq])
