import mdtraj as md
import numpy as np
import pickle
import os
import sys

def load_MDinfo(fname):
    try:
        os.path.islink('HISTORY')
    except:
        sys.exit('Symbolic link HISTORY does not exist.')
    sysxtc = os.readlink('HISTORY')

    with open(fname, 'rt') as f:
        ls = f.readline().split()
        Nspec = ls[0]
        Frames = ls[1]
        Nchains = f.readline().split()
        Nmons = f.readline().split()
    return sysxtc, Nspec, Frames, Nchains, Nmons


def judge_wrapped(pos):
    wrap_mode = 0
    for p in pos:
        p = np.array(p)
        for i in range(len(p)):
            TorF = p[i]>box
            if True in TorF:
                wrap_mode = 1
            TorF = p[i]<0.0
            if True in TorF:
                wrap_mode = 1
    return wrap_mode 


def write_comdirc(it, box, posc, posd):
    outdir = sysxtc.split('.xtc')[0].split('/')[-1] + '_pickles'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    if it == 1:
        with open(outdir+'/com_pos.pickle', mode='wb') as f:
            pickle.dump(posc, f)
        with open(outdir+'/dirc_pos.pickle', mode='wb') as f:
            pickle.dump(posd, f)
        with open(outdir+'/box.pickle', mode='wb') as f:
            pickle.dump(box, f)
    else:
        with open(outdir+'/com_pos.pickle', mode='ab') as f:
            pickle.dump(posc, f)
        with open(outdir+'/dirc_pos.pickle', mode='ab') as f:
            pickle.dump(posd, f)
        with open(outdir+'/box.pickle', mode='ab') as f:
            pickle.dump(box, f)



def write_comdirc_grofmt(it, box, posc, posd):
    outdir = sysxtc.split('.xtc')[0].split('/')[-1] + '_grofmt'
    outgro = outdir + sysxtc.split('.xtc')[0] + '_director_{0:05d}.gro',format(it)
    resnum = 1
    print(len(p), len(p[0]))
    numatm = len(p)*len(p[0])
    with open(outgro, 'wt') as f:
        l = 'director_gro\n{0:8d}\n'.format(numatm)
        f.write(l)
        for i,pi in enumerate(posc):
            for ii,pii in enumerate(pi):
                try:
                    dx = posd[i][ii,0]
                    dy = posd[i][ii,1]
                    dz = posd[i][ii,2]
                except:
                    dx,dy,dz = 0.0, 0.0, 0.0
                l = '{0:5d}PE_20   c3{1:5d}'.format(resnum, ((resnum-1)*len(p[0])+ii)% 99999 + 1) \
                    + '{0:8.3f}{1:8.3f}{2:8.3f}{3:8.3f}{4:8.3f}{5:8.3f}\n'.format(pii[0], pii[1], pii[2], dx, dy, dz)
                f.write(l)
            resnum += 1
        l = '{0:6.4f}\t{0:6.4f}\t{0:6.4f}\n'.format(box)
        f.write(l)


try:
    sysxtc, Nspec, Frames, Nchains, Nmons = load_MDinfo('./MDinfo')
except OSError:
    sys.exit('MDinfo does not exist!')

if not os.path.islink('HISTORY'):
    sys.exit(sysxtc + ' does not exist!')


Frames = int(Frames)
Nchain = int(Nchains[0])
Nmon = int(Nmons[0])
sysgro = './system_300.gro'


it = 1
for t in md.iterload(sysxtc, top=sysgro):
    if it > Frames:
        break
    print('{0:>5d}/{1:>5d})'.format(it, Frames))
    pos = t.xyz
    top = t.topology
    df, b = top.to_dataframe()
    pos = t.xyz
    box = t.unitcell_lengths[0,0]

    df['x'] = pos[0,:,0]
    df['y'] = pos[0,:,1]
    df['z'] = pos[0,:,2]

    c3s = ['' for i in range(Nchain)]
    indexes = ['' for i in range(Nchain)]
    c3s_pos = ['' for i in range(Nchain)]
    for j in range(Nchain):
        c3sj = df[df['name'] == 'c3' ]
        c3s[j] = c3sj[c3sj['resSeq'] == j+1]
        indexes[j] = np.array(c3s[j].index)
        c3s_pos[j] = pos[0][indexes[j]]

    # Check PBC treatment 
    if it == 1:
        wrap_mode = judge_wrapped(c3s_pos)
        if wrap_mode == 0:
            print('wrapped.')
        elif wrap_mode == 1:
            print('unwrapped.')


    xyz = ['x', 'y', 'z']
    com_pos = ['' for i in range(Nchain)]
    for i in range(Nchain):
        com_pos[i] = np.array([np.array([0.0e0 for k in xyz]) for j in range(2*Nmon-1)])
        if wrap_mode == 0:
            for l in range(len(c3s_pos[i])-1):
                print(i, l)
                print(c3s_pos[i][l+1])
                #print(c3s_pos[i][l])
                TorF1 = c3s_pos[i][l+1]-c3s_pos[i][l] > box/2.0
                scale = np.array([-1.0*int(x) for x in TorF1])
                b =  c3s_pos[i][l+1] + box*scale
                print(TorF1, scale)
                TorF2 = c3s_pos[i][l+1]-c3s_pos[i][l] < -1.0*box/2.0
                scale = np.array([int(x) for x in TorF2])
                b += + box*scale
                print(TorF2, scale)
                com_pos[i][l] = 0.50e0*(b+c3s_pos[i][l])
                print(b)
                print(c3s_pos[i][l])
                print(com_pos[i][l])

        elif wrap_mode == 1:
            for l in range(len(c3s_pos[i])-1):
                com_pos[i][l] = 0.50e0*(c3s_pos[i][l+1]+c3s_pos[i][l])
        Ncom = len(com_pos[0])

    dirc_pos = ['' for i in range(Nchain)]
    for i in range(Nchain):
        dirc_pos[i] = np.array([np.array([0.0e0 for k in xyz]) for j in range(Ncom-1)])
        if wrap_mode == 0:
            for l in range(Ncom-1):
                #print(i, l)
                #print(com_pos[i][l+1])
                TorF1 = com_pos[i][l+1]-com_pos[i][l] > box/2.0
                scale = np.array([-1.0*int(x) for x in TorF1])
                #print(TorF1, scale)
                b = com_pos[i][l+1] + box * scale
                TorF2 = com_pos[i][l+1]-com_pos[i][l] < -1.0*box/2.0
                scale = np.array([int(x) for x in TorF2])
                #print(TorF2, scale)
                b += + box * scale
                dirc_pos[i][l] = (b - com_pos[i][l])
                #print(b)
                #print(com_pos[i][l])
                #print(dirc_pos[i][l])
        elif wrap_mode == 1:
            for l in range(Ncom-1):    
                dirc_pos[i][l] = (com_pos[i][l+1]-com_pos[i][l])

    if wrap_mode == 1:
        for i in range(Nchain):
            # PBC treatment
            for k in range(len(xyz)):
                TorF = com_pos[i][:,k]>box
                tof0 = TorF[0]
                co = 0
                print(i,TorF)
                for tof_ in TorF:
                     if tof_ == True:
                         com_pos[i][co,k] -= box
                     if tof_ != tof0:
                         tof0 = tof_
                     co += 1

                TorF = com_pos[i][:,k]<0.0
                tof0 = TorF[0]
                co = 0
                for tof_ in TorF:
                    if tof_ == True:
                        com_pos[i][co,k] += box
                    if tof_ != tof0:
                        tof0 = tof_
                    co += 1
 
    write_comdirc(it, box, com_pos, dirc_pos)
    it += 1
