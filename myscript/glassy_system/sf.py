import numpy as np


fname = 'rdf.xvg'

rdfs = np.loadtxt(fname)

r = rdfs[:,0]
rdf = rdfs[:,1]

N = 10000
L = 6.773 #nm
rho = N/(L**3.0E0)
Nr = len(r)
dr = r[1]-r[0]
dk = 0.10
Nk = 5000

k = np.arange(dk, dk*(Nk+1),dk)

kr = k[:,np.newaxis]*r[np.newaxis,:]
sinkr_kr = np.sin(kr)/(kr)
print(sinkr_kr.shape)
Sk = np.zeros(Nk)
Sk += 1.0E0

rdf[int(0.6*Nr):] = 1.0E0
for ir in range(1,Nr):
    f = 0.50E0*(sinkr_kr[:,ir-1]+sinkr_kr[:,ir])
    g = 0.50E0*(rdf[ir-1]+rdf[ir])
    r2 = ((r[ir-1]+r[ir])/2.0)**2.0E0
    #print(ir,g)
    Sk += rho*4*np.pi*r2*dr*f*(g-1.0E0)


of = 'Sk.dat'
with open(of, 'wt') as f:
    for ik in range(Nk):
        l = '{0:10.6f}\t{1:8.6f}\n'.format(k[ik], Sk[ik])
        #print(k[ik], Sk[ik])    
        f.write(l)
