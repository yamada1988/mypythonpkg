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
#result:
#===num_frames:  100===
#Read trrfile time: 4.0275678634643555
#Calculate <rho(k,q)> time: 0.7862162590026855
#Calculate C(k,q,t) time: 0.0072154998779296875
#write file time: 0.006284236907958984
#total time: 4.827349662780762
 
#===num_frames:  200===
#Read trrfile time: 8.450029134750366
#Calculate <rho(k,q)> time: 1.577906608581543
#Calculate C(k,q,t) time: 0.027846574783325195
#write file time: 0.012517690658569336
#total time: 10.068350076675415
 
#===num_frames:  300===
#Read trrfile time: 12.38418173789978
#Calculate <rho(k,q)> time: 2.380192279815674
#Calculate C(k,q,t) time: 0.06324076652526855
#write file time: 0.016291379928588867
#total time: 14.843965291976929
 
#===num_frames:  400===
#Read trrfile time: 16.628856658935547
#Calculate <rho(k,q)> time: 3.169315814971924
#Calculate C(k,q,t) time: 0.10999298095703125
#write file time: 0.03341960906982422
#total time: 19.930981636047363
 
#===num_frames:  500===
#Read trrfile time: 20.03941035270691
#Calculate <rho(k,q)> time: 3.935575008392334
#Calculate C(k,q,t) time: 0.1746840476989746
#write file time: 0.026569604873657227
#total time: 24.1763014793396
 
#===num_frames:  600===
#Read trrfile time: 24.21495819091797
#Calculate <rho(k,q)> time: 4.74502158164978
#Calculate C(k,q,t) time: 0.2578587532043457
#write file time: 0.03341960906982422
#total time: 29.251312255859375
 
#===num_frames:  700===
#Read trrfile time: 28.1339271068573
#Calculate <rho(k,q)> time: 5.788179159164429
#Calculate C(k,q,t) time: 0.34491825103759766
#write file time: 0.03717231750488281
#total time: 34.30425810813904
 
#===num_frames:  800===
#Read trrfile time: 33.28687334060669
#Calculate <rho(k,q)> time: 6.744075775146484
#Calculate C(k,q,t) time: 0.4389834403991699
#write file time: 0.04378199577331543
#total time: 40.51376724243164
 
#===num_frames:  900===
#Read trrfile time: 36.22404408454895
#Calculate <rho(k,q)> time: 7.4806249141693115
#Calculate C(k,q,t) time: 0.5475361347198486
#write file time: 0.04688286781311035
#total time: 44.299171447753906
 
#===num_frames: 1000===
#Read trrfile time: 40.170594453811646
#Calculate <rho(k,q)> time: 8.021784543991089
#Calculate C(k,q,t) time: 0.6901025772094727
#write file time: 0.05596137046813965
#total time: 48.9384970664978
 
#===num_frames: 1500===
#Read trrfile time: 63.82360339164734
#Calculate <rho(k,q)> time: 12.624643325805664
#Calculate C(k,q,t) time: 1.56199049949646
#write file time: 0.08169889450073242
#total time: 78.09200096130371
 
#===num_frames: 2000===
#Read trrfile time: 81.69437336921692
#Calculate <rho(k,q)> time: 16.01133918762207
#Calculate C(k,q,t) time: 2.703766345977783
#write file time: 0.10345339775085449
#total time: 100.51300549507141



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
