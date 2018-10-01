from pytrr import GroTrrReader
import sys
import math
import numpy as np
import time

#
# Time measuring version
#
# time.sh: 
##!/bin/bash

#for t in 100 200 300 400 500 600 700 800 900 1000 1500 2000
#do

# t=`printf %5d $t`
# echo "===num_frames:${t}==="
# python script/calc_rhokq.py MD/sample0001_180.trr 1 1 $t
# echo " "
#done
#
#bash time.sh > timelog.dat
#

#
# Calculate rho(k,q,t), <rho(k,q)>, and <drho(k,q,t)drho(-k,-q,0)> from Gromacs trrfile.
# This script requires pytrr(https://github.com/andersle/pytrr).
# Install pytrr using "pip insall pytrr".
#

# Sample:
# python script/calc_rhokq.py MD/md_300.trr 1 2
#

args = sys.argv
fname = args[1]
T = float(fname.split('_')[1].split('.')[0])
nk = int(args[2])
nq = int(args[3])

pi = math.pi

# Mass Parameters
mO = 16.00 * 1.661e-27 #(kg)
mH =  1.08 * 1.661e-27 #(kg)

# Phisical Parameters
N =  3000
tN =  int(args[4])
dt = 0.010e0 #(ps)
L = 4.37864 #(nm) 
M = mO + mH + mH #(kg)
kB = 1.3801e-23 #(J K^-1)

# Wave number Parameters
k0 = 2.0e0*pi / L # (nm^-1)
dk = 1.0*2.0e0*pi / L # (nm^-1)
vth = math.sqrt(3.0e0*kB*T/M) # (nm/fs)^-1 = (km/s)^-1
#vth = 10.0e0*pi
alpha = 40.0
dq = alpha * 2.0e0*pi / vth
q0 = 0.0010
kN =  1
qN =  1

K = [1/math.sqrt(3) *np.array([k0+nk*dk, k0+nk*dk, k0+nk*dk]) for i in range(kN)]
Q = [1/math.sqrt(3) * np.array([q0+nq*dq, q0+nq*dq, q0+nq*dq]) for i in range(qN)]
for k in K:
    for ki in k:
        ki = round(ki, 6)
for q in Q:
    for qi in q:
        qi = round(qi, 6)

# rho_i(k,q,t,i)
rho_i = [[[[0.0e0 for i in range(N)] for t in range(tN+1)] for q in range(qN)] for k in range(kN)]

# rho(k,q,t)
rho = [[[0.0e0 for t in range(tN+1)] for q in range(qN)] for k in range(kN)]

R = [[0.0e0 for it in range(tN)] for i in range(N)]
V = [[0.0e0 for it in range(tN)] for i in range(N)]

t0 = time.time()
icount = -1
with GroTrrReader(fname) as trrfile:
    for frame in trrfile:
        icount += 1
        if icount >= tN:
            break
        frame_data = trrfile.get_data()
        for i in range(N):
            iO  = 3*i 
            iH1 = 3*i + 1
            iH2 = 3*i + 2

# Calculate single H2O molecule's center of mass (R)
            xO  = np.array(frame_data['x'][iO])
            xH1 = np.array(frame_data['x'][iH1])
            xH2 = np.array(frame_data['x'][iH2])
            R[i][icount] = (mO * xO + mH * xH1 + mH * xH2)/ M 

# Calculate single H2O molecule's velocity (V)
            vO  = np.array(frame_data['v'][iO])
            vH1 = np.array(frame_data['v'][iH1])
            vH2 = np.array(frame_data['v'][iH2])
            V[i][icount] = (mO * vO + mH * vH1 + mH * vH2)/ M

print("Read trrfile time:", time.time()-t0)
t1 = time.time()

# Calculate rho(k,q,t)
for it in range(tN):
    for iq,q in enumerate(Q):
        for ik,k in enumerate(K):
            rho[ik][iq][it] = 0.0e0
            for i in range(N):
                theta = np.dot(k,R[i][it]) + np.dot(q,V[i][it])
                rho_i[ik][iq][it][i] = complex(math.cos(theta), -math.sin(theta))
            rho[ik][iq][it] = np.sum(rho_i[ik][iq][it], axis=0)


# Calculate <rho(k,q)>
rho_ens = [[0.0e0 for iq in range(qN)] for ik in range(kN)]
for iq,q in enumerate(Q):
    for ik,k in enumerate(K):
        rho_ens[ik][iq] = np.sum(rho[ik][iq], axis=0)
        rho_ens[ik][iq] /= tN

print("Calculate <rho(k,q)> time:", time.time()-t1)


t2 = time.time()
# Calculate C(k,q,t) = <drho(k,q,t)drho(-k,q,0)>
thalf = int(tN/5)
C_kqt = np.array([[[complex(0.0e0, 0.0e0) for it in range(thalf+1)] for iq in range(qN)] for ik in range(kN)])
icount = [0 for i in range(thalf+1)]
for it0 in range(tN):
    for t_interval in range(tN-it0+1):
        if t_interval > thalf:
            continue
        icount[t_interval] += 1
        for ik in range(kN):
            for iq in range(qN):
                ct = rho[ik][iq][it0+t_interval]-rho_ens[ik][iq]
                c0 = rho[ik][iq][it0]-rho_ens[ik][iq] 
                C_kqt[ik][iq][t_interval] += ct * c0.conjugate()


for ik,k in enumerate(K):
    for iq,q in enumerate(Q):
        for t_interval in range(thalf+1):
            C_kqt[ik][iq][t_interval] /= icount[t_interval]
        for t_interval in range(1, thalf+1):
            C_kqt[ik][iq][t_interval] /= C_kqt[ik][iq][0]
        C_kqt[ik][iq][0] /= C_kqt[ik][iq][0]

print("Calculate C(k,q,t) time:", time.time()-t2)

t3 = time.time()
for ik,k in enumerate(K):
    for iq,q in enumerate(Q):
        ofname = 'DAT/Cnk{0:02d}nq{1:02d}_{2:d}.dat'.format(nk,nq,int(T))
        with open(ofname, 'wt') as f:
            f.write('# k=[{0[0]},{0[1]},{0[2]}]\n# q=[{1[0]},{1[1]},{1[2]}]\n# t(ps) C(k,q,t)\n'.format(k,q))

        for t_interval in range(thalf+1):
            t = '{0:5.2f}'.format(t_interval * dt)
            with open(ofname, 'a+') as f:
                f.write('{0:5.2f}\t{1:9.7f}\n'.format(float(t), C_kqt[ik][iq][t_interval]))

print("write file time:", time.time()-t3)
print("total time:", time.time()-t0)
