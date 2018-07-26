import mdtraj as md
from pandas import DataFrame

def printHbond(list, atom_pairs, frame, t, dt, index, selfflag=False):
    print('print HObonds...')
    if not selfflag:
        print('avoid self-chain...')
        for f in range(frame):
            time = t * frame * dt + f * dt
            for i in range(len(list[f])):
                if list[f][i] < 0.30:
                    print('pairs'+index+':{0:d}\t{1:d}'.format(atom_pairs[i][0], atom_pairs[i][1]))
                    print('time:{0:8.5f}ps \t r:{1:7.5f}nm'.format(time, list[f][i]))
    else:
        print('print self-chain...')
        for f in range(frame):
            time = t * frame * dt + f * dt
            for i in range(len(list[f])):
                if list[f][i] < 0.30 and abs(atom_pairs[i][0]- atom_pairs[i][1]) > n_monomer * 2:
                    print('pairs'+index+':{0:d}\t{1:d}'.format(atom_pairs[i][0], atom_pairs[i][1]))
                    print('time:{0:8.5f}ps \t r:{1:7.5f}nm'.format(time, list[f][i]))

def calcHbondinfo(list, atom_pairs, frame, t, dt, index, selfflag=False):
    print('check Hbonds of '+index+'...')
    if not selfflag:
        bondinfo = []
        for f in range(frame):
            time = t * frame * dt + f * dt
            for i in range(len(list[f])):
                if list[f][i] < 0.30 :
                    atom_pair1 = atom_pairs[i][0]
                    atom_pair2 = atom_pairs[i][1]
                    r = list[f][i]
                    line = '\t'.join(map(str, [time, atom_pair1, atom_pair2, r]))
                    bondinfo.append(line)         
    else:
        bondinfo = []
        for f in range(frame):
            time = t * frame * dt + f * dt
            for i in range(len(list[f])):
                if list[f][i] < 0.30 and abs(atom_pairs[i][0]- atom_pairs[i][1]) > n_monomer * 2:
                    atom_pair1 = atom_pairs[i][0]
                    atom_pair2 = atom_pairs[i][1]
                    r = list[f][i]
                    line = '\t'.join(map(str, [time, atom_pair1, atom_pair2, r]))
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
tot_frame = 500
nsep   = 10
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
    honame = 'DAT/hobondinfo_{0:03d}.dat'.format(k)
    ohname = 'DAT/ohbondinfo_{0:03d}.dat'.format(k)
    hoselfname = 'DAT/hobondselfinfo_{0:03d}.dat'.format(k)
    ohselfname = 'DAT/ohbondselfinfo_{0:03d}.dat'.format(k)

    str_ = '# data loaded in '+fname+'\n# time\tpair1\tpair2\tr (nm)'
    with open(honame, 'wt') as f:
        f.write(str_+'\n')
    with open(ohname, 'wt') as f:
        f.write(str_+'\n')
    with open(hoselfname, 'wt') as f:
        f.write(str_+'\n')
    with open(ohselfname, 'wt') as f:
        f.write(str_+'\n')


atom_indicesH_c = [ index for index in atom_indicesH if not index in ghost and index in core]
atom_indicesO_c = [ index for index in atom_indicesO if not index in ghost and index in core]

atom_indicesH_v = [ index for index in atom_indicesH if not index in ghost and not index in core]
atom_indicesO_v = [ index for index in atom_indicesO if not index in ghost and not index in core]

atom_pairsHO = topology.select_pairs(atom_indicesH_g, atom_indicesO_v) 
atom_pairsOH = topology.select_pairs(atom_indicesO_g, atom_indicesH_v)

atom_pairsHO_self = topology.select_pairs(atom_indicesH_g, atom_indicesO_c) 
atom_pairsOH_self = topology.select_pairs(atom_indicesO_g, atom_indicesH_c)

for k in range(i1, i2+1):
    fname  = 'MD/md{0:03d}_1.xtc'.format(k)
    honame = 'DAT/hobondinfo_{0:03d}.dat'.format(k)
    ohname = 'DAT/ohbondinfo_{0:03d}.dat'.format(k)
    hoselfname = 'DAT/hobondselfinfo_{0:03d}.dat'.format(k)
    ohselfname = 'DAT/ohbondselfinfo_{0:03d}.dat'.format(k)

    print('load trajectory: ' + fname+' ...')
    for t,chunk in enumerate(md.iterload(fname, top=tname, chunk=tot_frame/nsep)):
        t += 1
        print('chunk:', t)
        print('calculate rHO and rOH...')
        rHO = md.compute_distances(chunk, atom_pairsHO)
        rOH = md.compute_distances(chunk, atom_pairsOH)

        rHO_self = md.compute_distances(chunk, atom_pairsHO_self)
        rOH_self = md.compute_distances(chunk, atom_pairsOH_self)

        print('check hbond information...')
        frame = chunk.n_frames 
        print('frame:', frame)
#    printHbond(rHO, atom_pairsHO, frame, t, dt, 'HO')
#    printHbond(rOH, atom_pairsOH, frame, t, dt, 'OH')

#    printHbond(rHO_self, atom_pairsHO_self, frame, t, dt, 'HO_self', selfflag=True)
#    printHbond(rOH_self, atom_pairsOH_self, frame, t, dt, 'OH_self', selfflag=True)

        hbondHO = calcHbondinfo(rHO, atom_pairsHO, frame, t, dt, 'HO')
        writeHbond(hbondHO, honame) 

        hbondOH = calcHbondinfo(rOH, atom_pairsOH, frame, t, dt, 'OH')
        writeHbond(hbondOH, ohname) 

        hbondHOself = calcHbondinfo(rHO_self, atom_pairsHO_self, frame, t, dt, 'HO', selfflag=True)
        writeHbond(hbondHOself, hoselfname)

        hbondOHself = calcHbondinfo(rOH_self, atom_pairsOH_self, frame, t, dt, 'OH', selfflag=True)
        writeHbond(hbondOHself, ohselfname)

