import mdtraj as md
import numpy as np

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return 180.0e0/np.pi * np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def write_compos(pos):
    outgro = sysgro.split('.gro')[0] + '_com.gro'
    resnum = 1
    numatm = len(pos)*len(pos[0])
    with open(outgro, 'wt') as f:
        l = 'comgro\n{0:8d}\n'.format(numatm)
        f.write(l)
        for pi in pos:
            for ii,pii in enumerate(pi):
                l = '{0:5d}PVA10   c3{1:5d}'.format(resnum, ((resnum-1)*len(pos[0])+ii)% 99999 + 1) + '{0:8.3f}{1:8.3f}{2:8.3f}\n'.format(pii[0], pii[1], pii[2])
                f.write(l)
            resnum += 1
        l = '{0:6.4f}\t{0:6.4f}\t{0:6.4f}\n'.format(box)
        f.write(l)


sysgro = './sys/solution.gro'
Nchain = 100
Nmol = 100

t = md.load(sysgro)
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
    #print(indexes[j], pos[0][indexes[j]])
    c3s_pos[j] = pos[0][indexes[j]]
    #print(j,c3s_pos[j])

com_pos = ['' for i in range(Nchain)]
align_pos = ['' for i in range(Nchain)]
for i in range(Nchain):
    align_pos[i] = [[0.0e0 for k in ['x', 'y', 'z']] for j in range(200)]
    com_pos[i] = [[0.0e0 for k in ['x', 'y', 'z']] for j in range(200)]
    for l in range(len(c3s_pos[i])-1):
        #print(l, c3s_pos[i][l],c3s_pos[i][l+1])
        com_pos[i][l] = 0.50e0*(c3s_pos[i][l+1]+c3s_pos[i][l])
    for l in range(len(c3s_pos[i])-2):
        align_pos[i][l] = 0.50e0*(com_pos[i][l+1]-com_pos[i][l])

write_compos(com_pos)
 
thetas = ['' for i in range(Nchain)]
for i in range(Nchain):
    thetas[i] = np.array([0.0 for l in range(len(align_pos[i])-2)])
    for l in range(len(align_pos[i])-2):
        thetas[i][l] = angle_between(align_pos[i][l], align_pos[i][l+1])
        
    print(i, thetas[i]>40.0)


import sys
sys.exit()




nu[17:04:19]:~/calspa/OpenMM/POLYMER_MIXTURE/PVA/PVA100/2019-06-04/Iso-W/phi_0.85_grid>
nu[17:05:08]:~/calspa/OpenMM/POLYMER_MIXTURE/PVA/PVA100/2019-06-04/Iso-W/phi_0.85_grid>cat calc_clystalline.py 
import mdtraj as md
import numpy as np

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return 180.0e0/np.pi * np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

sysgro = './sys/solution_com.gro'
Nchain = 100
Nmol = 200

t = md.load(sysgro)
pos = t.xyz
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]

com = ['' for i in range(Nchain)]
indexes = ['' for i in range(Nchain)]
com_pos = ['' for i in range(Nchain)]
for j in range(Nchain):
    com_ = df[df['name'] == 'c3' ]
    com[j] = com_[com_['resSeq'] == j+1]
    indexes[j] = np.array(com[j].index)
    com_pos[j] = pos[0][indexes[j]]

align_pos = ['' for i in range(Nchain)]
for i in range(Nchain):
    align_pos[i] = [[0.0e0 for k in ['x', 'y', 'z']] for j in range(200)]
    for l in range(len(com_pos[i])-2):
        align_pos[i][l] = 0.50e0*(com_pos[i][l+1]-com_pos[i][l])

 
thetas = ['' for i in range(Nchain)]
for i in range(Nchain):
    thetas[i] = np.array([0.0 for l in range(len(align_pos[i])-2)])
    for l in range(len(align_pos[i])-2):
        thetas[i][l] = angle_between(align_pos[i][l], align_pos[i][l+1])
        
    #print(i, thetas[i]>40.0)

t0 = 30.0
Rods = [[[]] for i in range(Nchain)]
for i in range(Nchain):
    k = 0
    ik = 0
    v1 = align_pos[i][k]
    while True:
        if Rods[i][ik] == []:
            Rods[i][ik].append(k)
        v2 = align_pos[i][k+1]
        theta = angle_between(v1, v2)
        #print(i, ik, v1, v2, theta)
        if theta <= t0:
            Rods[i][ik].append(k+1)
            k += 1
        else:
            v1 = align_pos[i][k+1]
            ik += 1
            k += 1
            Rods[i].append([])
        if k+1 == len(align_pos[i]):
            del Rods[i][-1]
            break

for i in range(Nchain):
    for ik in range(len(Rods[i])):
        print(i, ik, Rods[i][ik])



import sys
sys.exit()
