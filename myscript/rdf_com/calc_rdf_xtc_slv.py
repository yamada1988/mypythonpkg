import mdtraj as md
import numpy as np
import sys
args = sys.argv
 
MassTable = {'C':12.0,'O':16.0,'Os':24.0,'VS':1.008,'H':1.008,'F':18.99}
target = 80
NP1  =  33
Nmol = 200
NP2  = Nmol - NP1
N1   = 100
N2   = 100
Ncom = NP1*N1+NP2*N2
Ncomslt = NP1*N1
Ncomslv = NP2*N2
Nchain1 = 602
Nchain2 = 602
dl1 = 5
dl2 = 5
Nskip = 1
Ntot = 100
Nit = Ntot/Nskip

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
        #print(m)
        R += m*r
        m0 += m
    R /= m0
    #print(R)
    return R

def write_f(fname, imol, index, Rvec):
    with open(fname, 'a+') as f:
        f.write('{0:5d}COM     c3{1:5d}{2:8.4f}{3:8.4f}{4:8.4f}\n'.format(imol, index, Rvec[0], Rvec[1], Rvec[2]))


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
    if it%Nskip != 0:
        continue
    if it == Ntot+1:
        break
    #if it > 2:
    #    sys.exit('test done.')

    Rs = [0.0e0 for i in range(Ncom)]
    top = t.topology
    df, b = top.to_dataframe()
    pos = t.xyz
    box = t.unitcell_lengths[0,0]

    df['x'] = pos[0,:,0]
    df['y'] = pos[0,:,1]
    df['z'] = pos[0,:,2]

    # calc rcom for P1
    for k in range(0, NP1):
        #print(k)
        ind = N1 * k + 1
        nc000 = Nchain1 * k
        nc00 = nc000 + 1
        nc0 = nc00 + dl1 + 1
        l0 = "a {0:d}-{1:d}\n".format(nc00,nc0)
        n1 = nc0
        #print(l0)
        df_ = df[nc00-1:nc0]
        #print(df_)
        Rs[ind-1] = calc_com(df_)
        #write_f(outf,k+1,ind,R)
    
        for i in range(1, N1-1):
            n0 = n1 + 1
            n1 = n0 + dl1
            l="a {0:d}-{1:d}\n".format(n0,n1)
            #print(l)
            df_ = df[n0-1:n1]
            #print(df_)
            ind += 1
            Rs[ind-1] = calc_com(df_)
            #write_f(outf,k+1,ind,R)
    
        ncN = nc000 + Nchain1 - dl1 - 2
        ncNN = nc00 + Nchain1 - 1
        lN0="a {0:d}-{1:d}\n".format(ncN,ncNN)
        #print(lN0)
        df_ = df[ncN:ncNN]
        #print(df_)
        ind += 1
        Rs[ind-1] = calc_com(df_)
        #write_f(outf,k+1,ind,R)
    ind_tmp = ind
    tmp = ncNN

    # calc rcom for P2
    for k in range(NP1, Nmol):
        #print(k)
        ind = ind_tmp + N2*(k-NP1) + 1
        nc000 = tmp + Nchain2 * (k-NP1)
        nc00 = nc000 + 1
        nc0 = nc00 + dl2 + 1
        l0 = "a {0:d}-{1:d}\n".format(nc00,nc0)
        n1 = nc0
        #print(l0)
        df_ = df[nc00-1:nc0]
        #print(df_)
        Rs[ind-1] = calc_com(df_)
        #write_f(outf,k+1,ind,R)
    
        for i in range(1, N2-1):
            n0 = n1 + 1
            n1 = n0 + dl2
            l="a {0:d}-{1:d}\n".format(n0,n1)
            #print(l)
            df_ = df[n0-1:n1]
            #print(df_)
            ind += 1
            Rs[ind-1] = calc_com(df_)
            #write_f(outf,k+1,ind,R)
    
        ncN = nc000 + Nchain2 - dl2 - 2
        ncNN = nc00 + Nchain2 - 1
        lN0="a {0:d}-{1:d}\n".format(ncN,ncNN)
        #print(lN0)
        df_ = df[ncN:ncNN]
        #print(df_)
        ind += 1
        Rs[ind-1] = calc_com(df_)
        #write_f(outf,k+1,ind,R)


    L = t.unitcell_lengths[0,0]
    #print(Rs)
    #print(L)
    rho = (Ncomslv)/(L**3)
    Lhalf = L/2.0e0
    Nmax = int(Lhalf/dr)
    r0 = Rs[target-1]
    r0 = np.array(r0)
    rslv = Rs[Ncomslt:]
    #npos = Ncomslt - N1
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
    #print(i)

for ir in range(Nmax):
    dv = 4.0e0/3.0e0 * np.pi * (ir**3-(ir-1)**3)*dr**3
    gr[ir] /= rho * dv * Nit

with open('DAT/rdf{0:04d}-PEPVDF.xvg'.format(num), 'wt') as f:
    f.write('# r \t gr\n')
with open('DAT/rdf{0:04d}-PEPVDF.xvg'.format(num), 'a+') as f:
    for ir in range(Nmax):
        f.write('{0:6.3f}\t{1:6.4f}\n'.format(ir*dr, gr[ir]))
