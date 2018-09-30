from pytrr import GroTrrReader
import sys
import math
import numpy as np

#
# Calculate rho(k,q,t), <rho(k,q)>, and <drho(k,q,t)drho(-k,-q,0)> from Gromacs trrfile.
# This script requires pytrr(https://github.com/andersle/pytrr).
# Install pytrr using "pip insall pytrr".
#

args = sys.argv
fname = args[1]
T = fname.split('_')[1].split('.')[0]

pi = math.pi

# Mass Parameters
mO = 16.00 * 1.661e-27 #(kg)
mH =  1.08 * 1.661e-27 #(kg)

# Phisical Parameters
N =  3000
tN = 1000
dt = 0.010e0 #(ps)
L = 4.37864 #(nm) 
#T = 200 #(K)
M = mO + mH + mH #(u)
kB = 1.3801e-23 #(J K^-1)

# Wave number Parameters
k0 = 2.0e0*pi / L # (nm^-1)
dk = 5.0*2.0e0*pi / L # (nm^-1)
#vth = math.sqrt(3.0e0*kB*T/M) # (nm/fs)^-1 = (km/s)^-1
vth = 10.0e0*pi
alpha = 5.0
dq = alpha * 2.0e0*pi / vth
q0 = 0.010
kN =  3
qN =  3

K = [1/math.sqrt(3) *np.array([k0+i*dk, k0+i*dk, k0+i*dk]) for i in range(kN)]
Q = [1/math.sqrt(3) * np.array([q0+i*dq, q0+i*dq, q0+i*dq]) for i in range(qN)]

# rho_i(k,q,t,i)
rho_i = [[[[0.0e0 for i in range(N)] for t in range(tN+1)] for q in range(qN)] for k in range(kN)]

# rho(k,q,t)
rho = [[[0.0e0 for t in range(tN+1)] for q in range(qN)] for k in range(kN)]

R = [[0.0e0 for it in range(tN)] for i in range(N)]
V = [[0.0e0 for it in range(tN)] for i in range(N)]

icount = -1

print('Read trrfile start.')
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

print('Read trrfile end.')
print('rho(k,q,t) calculation start.')
# Calculate rho(k,q,t)
for it in range(tN):
    for iq,q in enumerate(Q):
        for ik,k in enumerate(K):
            rho[ik][iq][it] = 0.0e0
            for i in range(N):
                theta = np.dot(k,R[i][it]) + np.dot(q,V[i][it])
                rho_i[ik][iq][it][i] = complex(math.cos(theta), -math.sin(theta))
            rho[ik][iq][it] = np.sum(rho_i[ik][iq][it], axis=0)/N
#                print(k,q,rho[ik][iq][icount])

# Calculate <rho(k,q)>
print("##### Ensemble Averaged density ({0:3d}K) ####".format(int(T)))
print('# k (wave number)                   q (vel fluc)                        rho(k,q)')
rho_ens = [[0.0e0 for iq in range(qN)] for ik in range(kN)]
for iq,q in enumerate(Q):
    for ik,k in enumerate(K):
        rho_ens[ik][iq] = np.sum(rho[ik][iq], axis=0)
        rho_ens[ik][iq] /= tN
        print(k,q,rho_ens[ik][iq])

# Calculate C(k,q,t) = <drho(k,q,t)drho(-k,q,0)>
thalf = int(tN/5)
C_kqt = np.array([[[0.0e0 for it in range(thalf+1)] for iq in range(qN)] for ik in range(kN)])
icount = [0 for i in range(thalf+1)]
for it0 in range(tN):
    for t_interval in range(tN-it0+1):
        if t_interval > thalf:
            continue
        icount[t_interval] += 1
        for ik in range(kN):
            for iq in range(qN):
                c = rho[ik][iq][it0]-rho_ens[ik][iq]
                C_kqt[ik][iq][t_interval] += (rho[ik][iq][it0+t_interval]-rho_ens[ik][iq]) * c.conjugate()
#                 c_t[it0] = (rho[ik][iq][it0] - rho_ens[ik][iq])
#                 c_0[it0] = (rho[ik][iq][it0] - rho_ens[ik][iq])


for ik,k in enumerate(K):
    for iq,q in enumerate(Q):
        for t_interval in range(thalf+1):
            C_kqt[ik][iq][t_interval] /= icount[t_interval]
        for t_interval in range(1, thalf+1):
            C_kqt[ik][iq][t_interval] /= C_kqt[ik][iq][0]
        C_kqt[ik][iq][0] /= C_kqt[ik][iq][0]

print("##### C(k,q,t) ({0:3d}K) ####".format(int(T)))
print('# k (wave number)                   q (vel fluc)                       t(ps)  C(k,q,t)')
for ik,k in enumerate(K):
    for iq,q in enumerate(Q):
        if ik ==kN-1 and iq == qN-1:
            with open('Ckqt_{0:3d}.dat'.format(int(T)), 'wt') as f:
                f.write('#k=[{0[0]},{0[1]},{0[2]}] q=[{1[0]},{1[1]},{1[2]}]\n'.format(k,q))

        for t_interval in range(thalf+1):
            t = '{0:5.2f}'.format(t_interval * dt)
            print(k,q,t,C_kqt[ik][iq][t_interval])
            if ik == kN-1 and iq == qN-1:
                with open('Ckqt_{0:3d}.dat'.format(int(T)), 'a+') as f:
                    f.write('{0:5.2f}\t{1:9.7f}\n'.format(float(t), C_kqt[ik][iq][t_interval]))
