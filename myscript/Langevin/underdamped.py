import numpy as np
import sys
import random

args = sys.argv
index = args[1]

N = 10
L = 50.0E0
T = 1.0E0
U0 = 13.0*T
W0 = 13.0*T
Y0 = 13.0*T
V0 = 5.0*T
k0 = 0.0*T
gamma = 1.0E0
m = 1.0E0
dt = 1.0e-5
f =  5.0E0
k1 = 3.0E0
k2 = 7.0E0
k3 = 5.0E0
sigma = L/100.0E0
drp = 1.0e-12
drv = 1.0e-15
TN =  250000000
T0 =       5000
lmbda = np.array([0.50]*N)
#lmbda = np.array([1.0,1.0,1.0,1.0,1.0])*L
r0 = 0.5*L
rs = np.array([0.0,0.0,0.0,0.0,0.0])*L
#TN =   1000000
#T0 =       100

for prpty in ['position', 'momentum', 'heat_current', 'escape_rate', 'energy_current', 'force_pot', 'force_int']:
    outf = 'DAT/status-stratonovich-'+prpty+index+'.data'
    line = '#time\t\t'
    with open(outf, 'wt') as of:
        for i in range(N):
            line += "p{0:02d}\t".format(i)
        line += "\n"
        of.write(line)

def U(l0,r):
    return l0*U0*np.cos(2.0*k1*np.pi*r/L)+(1.0-l0)*W0*np.cos(2.0*k2*np.pi*r/L)+Y0*np.cos(2.0*k3*np.pi*r/L)

def Fp(l0,rx):
    Urpdr = U(l0,rx+drp)
    Urmdr = U(l0,rx-drp)
    return -1.0E0*(Urpdr-Urmdr)/(2.0E0*drp)

def V(rv):
    return V0*4.0*((sigma/rv)**12.0E0+1.0E0*(sigma/rv)**6.0E0)

def Fv(r):
    Vrpdr = V(r+drv)
    Vrmdr = V(r-drv)
    return -1.0E0*(Vrpdr-Vrmdr)/(2.0E0*drv)

def set_initpos(N0,L0,dx0):
    Nbin0 = int(L0/dx0)
    pool = [i for i in range(Nbin0)]
    ix0 = random.sample(pool, N0)
    x0 = np.array(ix0)*dx0
    return x0

#initial condition
dx = 2*sigma
x0 = set_initpos(N, L, dx)
p0 = np.array([1.0E0*np.random.normal() for i in range(N)])
xi0 = np.array([np.sqrt(2.0*gamma*T)*np.random.normal() for i in range(N)])
#x0 = np.array([2.50,2.0])

fij = np.zeros((N,N))
x = np.zeros(N)
p = np.zeros(N)
fint = np.zeros(N)
fp = np.zeros(N)

for it in range(1,TN):
    t = it*dt
    fint0 = fint
    fint = np.zeros(N)
    Wfv0 = np.zeros(N)
    J0 = np.zeros(N)
    K0 = np.zeros(N)
    f0 = np.zeros(N)

    xi = np.array([np.sqrt(2.0*gamma*T)*np.random.normal() for i in range(N)])
    for i in range(N):
        for j in range(i+1,N):
            xij = x0[j] - x0[i]
            xij -= np.round(xij/L)*L
            rij = np.sqrt(np.sum(xij**2.0E0))
            
            fij[i,j] = Fv(rij)*(-1.0E0)*xij/rij
            fij[j,i] = -fij[i,j]
            fint[i] += fij[i,j]
            fint[j] += fij[j,i]
        fp[i] = Fp(lmbda[i],x0[i])
        f0[i] = -gamma*(p0[i]/m) + f + fp[i] + fint[i] 
        p[i] = p0[i] + f0[i]*dt + xi[i]*dt**0.50E0
        x[i] = x0[i] + (p0[i]/m)*dt  

        if x[i]>L:
            x[i]-=L
        elif x[i]<0:
             x[i]+=L

        if it == 1:
            continue

    t0 = t-dt
    Wfv0 = f*(p0/m)
    J0 = gamma*(p0/m)*(p0/m) - (gamma*T/m) - xi0*(p0/m)
    je0 = Wfv0 - J0
    K0 = np.sum(p0**2.0E0)/(2.0E0*m)
    Vint0 = 0.0E0
    for i in range(N):
        for j in range(i+1,N):
            xij0 = x0[i]-x0[j]
            xij0 -= np.round(xij0/L)*L
            rij0 = np.sqrt(np.sum(xij0**2.0E0))
            Vint0 += V(rij0)
    Ux0 = 0.0E0
    for i in range(N):
        Ux0 += U(lmbda[i],x0[i])
    Hv0 = Ux0 + Vint0
    E0 = K0 + Hv0
    er0 = (((p-p0)/dt)**2.0E0)/(2.0E0*gamma)+(gamma/m)*(p0**2.0E0)/(2.0E0*m)+(f0**2.0E0)/(2.0E0*gamma)-T*(gamma/m)

    if it > T0 and it % 5000 == 0:
        for prpty in ['position', 'momentum', 'heat_current', 'escape_rate', 'energy_current', 'force_pot', 'force_int']:
            outf = 'DAT/status-stratonovich-'+prpty+index+'.data'
            if prpty == 'position':
                data = x0
            elif prpty == 'momentum':
                data = p0
            elif prpty == 'heat_current':
                data = J0
            elif prpty == 'escape_rate':
                data = er0
            elif prpty == 'energy_current':
                data = je0
            elif prpty == 'force_pot':
                data = fp0
            elif prpty == 'force_int':
                data = fint0
            with open(outf, 'a+') as of:
                status = '{0:.6e}'.format(t0)
                status += '\t{0:5.3f}\t{1:5.3f}\t{2:5.3f}\t{3:5.3f}\t{4:5.3f}'.format(data[0],data[1],data[2],data[3],data[4])
                status += '\t{0:5.3f}\t{1:5.3f}\t{2:5.3f}\t{3:5.3f}\t{4:5.3f}'.format(data[5],data[6],data[7],data[8],data[9])
                #status += '\t{0:5.3f}\t{1:5.3f}\t{2:5.3f}\t{3:5.3f}\t{4:5.3f}'.format(data[10],data[11],data[12],data[13],data[14])
                #status += '\t{0:5.3f}\t{1:5.3f}\t{2:5.3f}\t{3:5.3f}\t{4:5.3f}'.format(data[15],data[16],data[17],data[18],data[19])
                #print(status) 
                status += '\n'
                of.write(status)


    if Vint0 > 1000000:
        sys.exit()
    x0 = x
    p0 = p
    xi0 = xi   
    fp0 = fp
