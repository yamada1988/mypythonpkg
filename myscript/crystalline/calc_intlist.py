import numpy as np
import pickle
from joblib import Parallel, delayed
from scipy import interpolate
import os
import sys

def load_dumps(f):
    obj = []
    while 1:
        try:
            obj.append(pickle.load(f))
        except:
            break
    return obj


def calc_gr(r0, dr, dists, N, rho):
    dv = 4.0*np.pi*r0**(2.0)*dr
    if r0 > 0.40:
        return float(np.count_nonzero((r0<=d) & (d<r0+dr)))/(rho*dv*N/2.0e0)
    else:
        return 0.0


def spline_3d(x, y):
    print(len(x), len(y))
    spline_func = interpolate.CubicSpline(x, y)
    print(spline_func.c)
    x_res = np.linspace(0, x[-1], 500)
    y_res = spline_func(x_res)

    return x_res, y_res


def write_rdf(rmin, x, y, outdir, fname='rdf.xvg'):
    with open(outdir+fname, 'wt') as f:
        f.write('# d = {0:5.4f} nm\n'.format(float(rmin)))
        f.write('# r \t gr\n')
    with open(outdir+fname, 'a+') as f:
        for ir in range(len(x)):
            f.write('{0:6.3f}\t{1:6.4f}\n'.format(x[ir], y[ir]))


def write_intlist(it_, rm, outdir, obj):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    if it_ == 0:
        with open(outdir+'/intlist.pickle', mode='wb') as f:
            pickle.dump(rm, f)
        with open(outdir+'/intlist.pickle', mode='ab') as f:
            pickle.dump(obj, f)
    else:
        with open(outdir+'/intlist.pickle', mode='ab') as f:
            pickle.dump(obj, f)


args = sys.argv
mode = args[1].split('-')[1]

pdir = 'npt_r_gpu_pickles/'
posb = pdir + 'box.pickle'
posc = pdir + 'com_pos.pickle'
dr = 0.010 #nm
L = 20.0 #nm
Lh = L/2.0 #nm
d0 = 1.0 #nm
r_ = np.array([(i+0.50)*dr for i in range(int(Lh/dr))])
ind0 = int(d0/dr)
with open(posb, 'rb') as f:
    boxes = load_dumps(f)

with open(posc, 'rb') as f:
    com_pos = load_dumps(f)

Nchain = len(com_pos[0])
Ncom = len(com_pos[0][0])
N = Nchain * Ncom
Frames = len(boxes)

ds = ['' for it in range(Frames)]
grs = ['' for it in range(Frames)]
for it in range(Frames):
    #print(it)
    box = boxes[it]
    r_ = r_[:int(0.50*box/dr)]
    rho = N/(box**3)
    c_pos = np.array(com_pos[it])
    d_pos = c_pos[:,:,np.newaxis, np.newaxis] - c_pos[np.newaxis, np.newaxis, :, :]
    d_pos -= np.trunc(d_pos/(box/2.0))
    d = np.sqrt(np.sum(d_pos**2.0, axis=4))

    if mode in ['gs', 'sg']:
        # Serial calculation
        gr = [float(np.count_nonzero((r<=d) & (d<r+dr)))/(rho*4.0/3.0*np.pi*r**3.0*dr*N/2.0e0) for r in r_]
    elif mode in ['gp', 'pg']:
        # Parallel calulation
        result = Parallel(n_jobs=-1, verbose=10, backend="threading")([delayed(calc_gr)(r, dr, d, N, rho) for r in r_] )
        a = np.array(result)
        grs[it] = np.array(result)
        #for ir, r in enumerate(r_):
        #    print(r, grs[it][ir])
    ds[it] = d

if'g' in mode:
    gr = np.array([0.0 for r in r_])
    for ir, r in enumerate(r_):
       for it in range(Frames):
           gr[ir] += grs[it][ir]
    gr /= Frames
    
    Nmax =  int((np.min(boxes)/2.0)/dr)
    r_sp, gr_sp = spline_3d(r_[:Nmax], gr[:Nmax])
    #for ir in range(len(r_sp)):
    #    print(r_sp[ir], gr_sp[ir])
    
    gmax = np.amax(gr_sp)
    maxid = np.where(gr_sp >= gmax)[0][0]
    from scipy import signal
    #print(gr_sp[maxid:maxid+100])
    minId = signal.argrelmin(gr_sp[maxid:maxid+100])
    minId += maxid
    #print(r_sp[minId])
    #print(gr_sp[minId])
    gmin = gr_sp[minId][0][0]
    minid = 0
    for ig, g_ in enumerate(gr_sp[minId][0]):
        if g_ < gmin:
            gmin = g_
            minid = ig
    
    #print(gmin)
    # not clustered
    d_ = r_sp[minId][0][minid]
    # clustered
    d_ = 1.0

    #print('gsp-max:', gmax)
    #print('maxid-sp:', maxid)
    #print('gr(spline)]', gr_sp[maxid:-50])
    #print('minid-sp:', minid)
    #print(r_sp[maxid], gr_sp[maxid], d_, gmin)
    minid = int(d_/dr)
    #print(minid, r_[minid])

    write_rdf(d_, r_[:Nmax], gr[:Nmax], pdir, fname='rdf.xvg')
    write_rdf(d_, r_sp, gr_sp, pdir, fname='rdf-spline3d.xvg')

#clustered
d_ = 1.0

# truncate intlist[r>d]
for it,d in enumerate(ds):
    #print('it:', it)
    inds = np.where( (0.20<d) & (d<d_))
    intlist = zip(*inds)
    write_intlist(it, d_, pdir, intlist)
