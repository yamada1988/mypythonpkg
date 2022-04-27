import mdtraj as md
from pandas import DataFrame
import os

def printHbond(list, atom_pairs, frame, t, dt, index, selfflag=False):
    print('print HObonds...')
    if not selfflag:
        print('avoid self-chain...')
        for f in range(frame):
            time = (t-1) * frame * dt + f * dt
            for i in range(len(list[f])):
                if list[f][i] < 0.30:
                    print('pairs'+index+':{0:d}\t{1:d}'.format(atom_pairs[i][0], atom_pairs[i][1]))
                    print('time:{0:8.4f}ps \t r:{1:7.4f}nm'.format(time, list[f][i]))
    else:
        print('print self-chain...')
        for f in range(frame):
            time = (t-1) * frame * dt + f * dt
            for i in range(len(list[f])):
                if list[f][i] < 0.30 and abs(atom_pairs[i][0]- atom_pairs[i][1]) > n_monomer * 2:
                    print('pairs'+index+':{0:d}\t{1:d}'.format(atom_pairs[i][0], atom_pairs[i][1]))
                    print('time:{0:8.4f}ps \t r:{1:7.4f}nm'.format(time, list[f][i]))

def make_nopbc(list, Lbox):
    for i in [0, 1, 2]:
        if list[i] > Lbox:
            list[i] -= Lbox
        elif list[i] < 0:
            list[i] += Lbox
    return list

    
def calcHbondinfo(traj, list, atom_pairs, frame, t, dt, index, selfflag=False):
    print('check Hbonds of '+index+'...')
    if not selfflag:
        bondinfo = []
        for f in range(frame):
            time = (t-1) * frame * dt + f * dt
            Lbox = traj.unitcell_vectors[f][0][0]
            for i in range(len(list[f])):
                if list[f][i] < 0.30 :
                    atom_pair1 = atom_pairs[i][0]
                    atom_pair2 = atom_pairs[i][1]
                    r = '{0:5.3f}'.format(list[f][i])
                    r1 = make_nopbc(traj.xyz[f][atom_pair1], Lbox)
                    r2 = make_nopbc(traj.xyz[f][atom_pair2], Lbox)
                    list_ = [time, atom_pair1, atom_pair2, r, r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], Lbox]
                    line = '\t'.join(map(str, list_))
                    bondinfo.append(line)         
    else:
        bondinfo = []
        for f in range(frame):
            time = (t-1) * frame * dt + f * dt
            Lbox = traj.unitcell_vectors[f][0][0]
            for i in range(len(list[f])):
                if list[f][i] < 0.30 and abs(atom_pairs[i][0]- atom_pairs[i][1]) > n_monomer * 2:
                    atom_pair1 = atom_pairs[i][0]
                    atom_pair2 = atom_pairs[i][1]
                    r = '{0:6.4f}'.format(list[f][i])
                    r1 = make_nopbc(traj.xyz[f][atom_pair1], Lbox)
                    r2 = make_nopbc(traj.xyz[f][atom_pair2], Lbox)
                    list_ = [time, atom_pair1, atom_pair2, r, r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], Lbox]
                    line = '\t'.join(map(str, list_))
                    bondinfo.append(line) 
    return bondinfo


def writeHbond(list, fname):
    with open(fname, 'a+') as f:
        for i,line in enumerate(list):
            f.write(line+'\n')


dt     = 0.10
nstart = 1
nend   = 16
ncore  = 752
n_monomer = 15
tot_frame = 5000
nsep   = 500
i1 = 1
i2 = 1
tname  = '../pdbdir/em.pdb'

print('load topology: '+tname+' ...')
pdb = md.load_pdb(tname)
topology = pdb.topology
td = topology.to_dataframe()

print('select atom-pairs for hydrogen bonding...')
atom_indicesH = [a.index for a in topology.atoms if a.element.symbol == 'H' ]
atom_indicesO = [a.index for a in topology.atoms if a.element.symbol == 'O' ]

ghost = [i for i in range(nstart, nend+1)]
core  = [i for i in range(nend+1, ncore)]

atom_indicesH_g = [ index for index in atom_indicesH if index in ghost ]
atom_indicesO_g = [ index for index in atom_indicesO if index in ghost ]
for k in range(i1, i2+1):
    fname  = 'MD/md{0:03d}_1.xtc'.format(k)
    honame = 'DAT/hobondnetworkinfo_{0:03d}.dat'.format(k)
    ohname = 'DAT/ohbondnetworkinfo_{0:03d}.dat'.format(k)

    str_ = '# data loaded in '+fname+'\n# time\tpair1\tpair2\tr (nm)\tx1\ty1\tz1\tx2\ty2\tz2\tLbox'
    for fname in [honame, ohname]:
        with open(fname, 'wt') as f:
            f.write(str_+'\n')


atom_indicesH_c = [ index for index in atom_indicesH if not index in ghost and index in core]
atom_indicesO_c = [ index for index in atom_indicesO if not index in ghost and index in core]

atom_indicesH_v = [ index for index in atom_indicesH if not index in ghost and not index in core]
atom_indicesO_v = [ index for index in atom_indicesO if not index in ghost and not index in core]

atom_indicesH = atom_indicesH_g + atom_indicesH_c + atom_indicesH_v
atom_indicesO = atom_indicesO_g + atom_indicesO_c + atom_indicesO_v

atom_pairsHO = topology.select_pairs(atom_indicesH, atom_indicesO) 
atom_pairsOH = topology.select_pairs(atom_indicesO, atom_indicesH)


for k in range(i1, i2+1):
    fname  = 'MD/md{0:03d}_1.xtc'.format(k)
    honame = 'DAT/hobondnetworkinfo_{0:03d}.dat'.format(k)
    ohname = 'DAT/ohbondnetworkinfo_{0:03d}.dat'.format(k)

    if not os.path.exists(fname):
        print(fname+' does not exist! skip this xtc file...')
        continue

    print('load trajectory: ' + fname+' ...')
    for t,chunk in enumerate(md.iterload(fname, top=tname, chunk=tot_frame/nsep)):
        t += 1
        print('chunk:', t)
        print('calculate rHO and rOH...')
        rHO = md.compute_distances(chunk, atom_pairsHO)
        rOH = md.compute_distances(chunk, atom_pairsOH)

        print('check hbond information...')
        frame = chunk.n_frames 
        print('frame:', frame)
#    printHbond(rHO, atom_pairsHO, frame, t, dt, 'HO')
#    printHbond(rOH, atom_pairsOH, frame, t, dt, 'OH')

#    printHbond(rHO_self, atom_pairsHO_self, frame, t, dt, 'HO_self', selfflag=True)
#    printHbond(rOH_self, atom_pairsOH_self, frame, t, dt, 'OH_self', selfflag=True)

        hbondHO = calcHbondinfo(chunk, rHO, atom_pairsHO, frame, t, dt, 'HO', selfflag = True)
        writeHbond(hbondHO, honame) 

        hbondOH = calcHbondinfo(chunk, rOH, atom_pairsOH, frame, t, dt, 'OH', selfflag = True)
        writeHbond(hbondOH, ohname) 

