import mdtraj as md
import numpy as np

def write_compos(pos):
    outgro = sysgro.split('.gro')[0] + '_com.gro'
    resnum = 1
    print(len(pos), len(pos[0]))
    numatm = len(pos)*len(pos[0])
    with open(outgro, 'wt') as f:
        l = 'comgro\n{0:8d}\n'.format(numatm)
        f.write(l)
        for pi in pos:
            for ii,pii in enumerate(pi):
                l = '{0:5d}PE_20   c3{1:5d}'.format(resnum, ((resnum-1)*len(pos[0])+ii)% 99999 + 1) + '{0:8.3f}{1:8.3f}{2:8.3f}\n'.format(pii[0], pii[1], pii[2])
                f.write(l)
            resnum += 1
        l = '{0:6.4f}\t{0:6.4f}\t{0:6.4f}\n'.format(box)
        f.write(l)


def write_pointsdirc(p, a):
    outgro = sysgro.split('.gro')[0] + '_director.gro'
    resnum = 1
    print(len(p), len(p[0]))
    numatm = len(p)*len(p[0])
    with open(outgro, 'wt') as f:
        l = 'director_gro\n{0:8d}\n'.format(numatm)
        f.write(l)
        for i,pi in enumerate(p):
            for ii,pii in enumerate(pi):
                l = '{0:5d}PE_20   c3{1:5d}'.format(resnum, ((resnum-1)*len(p[0])+ii)% 99999 + 1) \
                    + '{0:8.3f}{1:8.3f}{2:8.3f}{3:8.3f}{4:8.3f}{5:8.3f}\n'.format(pii[0], pii[1], pii[2], a[i][ii,0], a[i][ii,1],a[i][ii,2])
                f.write(l)
            resnum += 1
        l = '{0:6.4f}\t{0:6.4f}\t{0:6.4f}\n'.format(box)
        f.write(l)

sysgro = './md_run.gro'
Nchain = 72
Nmol = 20

t = md.load(sysgro)
pos = t.xyz
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]
print(df)

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]

c3s = ['' for i in range(Nchain)]
indexes = ['' for i in range(Nchain)]
c3s_pos = ['' for i in range(Nchain)]
for j in range(Nchain):
    c3sj = df[df['name'] == 'c3' ]
    c3s[j] = c3sj[c3sj['resSeq'] == j+2]
    indexes[j] = np.array(c3s[j].index)
    #print(indexes[j], pos[0][indexes[j]])
    c3s_pos[j] = pos[0][indexes[j]]
    #print(j,c3s_pos[j])

com_pos = ['' for i in range(Nchain)]
for i in range(Nchain):
    com_pos[i] = np.array([np.array([0.0e0 for k in ['x', 'y', 'z']]) for j in range(2*Nmol-1)])
    for l in range(len(c3s_pos[i])-1):
        #print(l, c3s_pos[i][l],c3s_pos[i][l+1])
        com_pos[i][l] = 0.50e0*(c3s_pos[i][l+1]+c3s_pos[i][l])
    #print(i,com_pos[i])

write_compos(com_pos)


Ncom = len(com_pos[0])
xyz = ['x', 'y', 'z']

align_pos = ['' for i in range(Nchain)]
for i in range(Nchain):
    align_pos[i] = np.array([[0.0e0 for k in xyz] for j in range(Ncom-1)])
    for l in range(Ncom-1):
        align_pos[i][l] = com_pos[i][l+1]-com_pos[i][l]

    # PBC treatment
    for k in range(len(xyz)):
        TorF = com_pos[i][:,k]>box
        tof0 = TorF[0]
        co = 0
        for tof_ in TorF:
            if tof_ == True:
                com_pos[i][co,k] -= box
            if tof_ != tof0:
                align_pos[i][co-1] = np.zeros(3)
                tof0 = tof_
            co += 1
 
        TorF = com_pos[i][:,k]<0.0
        tof0 = TorF[0]
        co = 0
        for tof_ in TorF:
            if tof_ == True:
                com_pos[i][co,k] += box
            if tof_ != tof0:
                align_pos[i][co-1] = np.zeros(3)
                tof0 = tof_
            co += 1
 

points_pos = ['' for i in range(Nchain)]
for i in range(Nchain):
    points_pos[i] = np.array([np.array([0.0e0 for k in ['x', 'y', 'z']]) for j in range(Ncom-1)])
    for l in range(Ncom-1):
        points_pos[i][l] = 0.50e0*(com_pos[i][l+1]+com_pos[i][l])
        #print(l, points_pos[i][l], align_pos[i][l])

write_pointsdirc(points_pos, align_pos)
