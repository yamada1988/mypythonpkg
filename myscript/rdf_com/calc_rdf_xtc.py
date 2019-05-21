import mdtraj as md
import numpy as np
import sys
args = sys.argv
 
MassTable = {'C':12.0,'O':24.0,'Os':24.0,'VS':1.008,'H':1.008}
target = 100
dr = 0.020 #nm
Lstd = 12.0
gr = np.array([0.0e0 for i in range(int(Lstd/2.0e0*1/dr))])


def calc_com(dataf):
    R = np.array([0.0e0, 0.0e0, 0.0e0])
    m0 = 0.0e0
    for index, d in dataf.iterrows():
        #print(d)
        e = d['element']
        m = MassTable[e]
        r = np.array([d['x'], d['y'], d['z']])
        #print(e,m)
        R += m*r
        m0 += m
    R /= m0
    print(R)
    return R

def write_f(fname, imol, index, Rvec):
    with open(fname, 'a+') as f:
        f.write('{0:5d}COM     c3{1:5d}{2:8.4f}{3:8.4f}{4:8.4f}\n'.format(imol, index, Rvec[0], Rvec[1], Rvec[2]))

Nmol = 100
Nchain = 1502
dl = 14

fname = args[1]
topname = args[2]
outf = topname.split('system')[0] + 'com' + topname.split('system')[1]
num = topname.split('system')[1].split('.gro')[0]
num = int(num)

with open(outf, 'wt') as f:
    f.write('Center of Mass for each monomers\n')
    f.write('{0:7d}\n'.format(Nmol * 100))

ts = md.load(fname,top=topname)
it = 0
ng = 50
for t in ts:
    it += 1
    print(it)
    Rs = [0.0e0 for i in range(10000)]
    top = t.topology
    df, b = top.to_dataframe()
    pos = t.xyz
    box = t.unitcell_lengths[0,0]

    df['x'] = pos[0,:,0]
    df['y'] = pos[0,:,1]
    df['z'] = pos[0,:,2]

    for k in range(0,Nmol):
        ind = 100 * k + 1
        nc000 = Nchain * k
        nc00 = nc000 + 1
        nc0 = nc00 + dl + 1
        l0 = "a {0:d}-{1:d}\n".format(nc00,nc0)
        n1 = nc0
        #print(l0)
        df_ = df[nc00-1:nc0]
        #print(df_)
        Rs[ind-1] = calc_com(df_)
        #write_f(outf,k+1,ind,R)
    
        for i in range(1,98+1):
            n0 = n1 + 1
            n1 = n0 + dl
            l="a {0:d}-{1:d}\n".format(n0,n1)
            #print(l)
            df_ = df[n0-1:n1]
            #print(df_)
            ind += 1
            Rs[ind-1] = calc_com(df_)
            #write_f(outf,k+1,ind,R)
    
        ncN = nc000 + Nchain - dl - 2
        ncNN = nc00 + Nchain - 1
        lN0="a {0:d}-{1:d}\n".format(ncN,ncNN)
        #print(lN0)
        df_ = df[ncN:ncNN]
        #print(df_)
        ind += 1
        Rs[ind-1] = calc_com(df_)
        #write_f(outf,k+1,ind,R)

    L = t.unitcell_lengths[0,0]
    print(Rs)
    print(L)
    rho = 10000/(L**3)
    Lhalf = L/2.0e0
    Nmax = int(Lhalf/dr)
    r0 = Rs[target-1]
    r0 = np.array(r0)
    rslv = Rs[target:]
    npos = 10000 - 100
    for rs in rslv:
        rs = np.array(rs)
        r = rs - r0
        #print(r)
        #print(np.round(r/L)) 
        r -= np.round(r/L)*L
        #print(r)
        rnorm = np.linalg.norm(r)
        index = int(rnorm/dr)
        #print(rnorm, index)
        if index <= Nmax:
            gr[index] += 1
    print(i)

for ir in range(Nmax):
    dv = 4.0e0/3.0e0 * np.pi * (ir**3-(ir-1)**3)*dr**3
    gr[ir] /= rho * dv * it

with open('rdf{0:04d}.xvg'.format(num), 'wt') as f:
    f.write('# r \t gr\n')
with open('rdf{0:04d}.xvg'.format(num), 'a+') as f:
    for ir in range(Nmax):
        f.write('{0:6.3f}\t{1:6.4f}\n'.format(ir*dr, gr[ir]))
