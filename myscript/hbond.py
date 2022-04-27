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

def calcHbond(list, atom_indices, atom_pairs, frame, t, dt, index, selfflag=False):
    print('check Hbonds...')
    if not selfflag:
        nbond = [[0 for raw in range(len(atom_indices)+1)] for f in range(frame)]
        for f in range(frame):
            time = t * frame * dt + f * dt
            nbond[f][0] = time
            for i in range(len(list[f])):
                for j,atm in enumerate(atom_indices):
                    if list[f][i] < 0.30 and atom_pairs[i][0] == atm:
                        nbond[f][j+1] += 1
    else:
        nbond = [[0 for raw in range(len(atom_indices)+1)] for f in range(frame)]
        for f in range(frame):
            time = t * frame * dt + f * dt
            nbond[f][0] = time
            for i in range(len(list[f])):
                for j,atm in enumerate(atom_indices):
                    if list[f][i] < 0.30 and abs(atom_pairs[i][0]- atom_pairs[i][1]) > n_monomer * 2 and atom_pairs[i][0] == atm:
                        nbond[f][j+1] += 1
    return nbond


def writeHbond(list, fname):
    with open(fname, 'a+') as f:
        for i,line in enumerate(list):
            str_ = '\t'.join(map(str, list[i]))
            f.write(str_+'\n')


def calc_ave(fname):
    with open(fname, 'rt') as f:
        lines = [map(float, line.split()) for line in f if not line.startswith('#')]
        n_data = len(lines)
        aves = []
        for r in range(1, len(lines[0])):
            raw = [raw[r] for raw in lines]
            print('raw')
            print(raw)
            sum_ = sum(raw)
            ave_ = sum_ / n_data
            aves.append(ave_)
        print(aves)
    with open(fname, 'a+') as f:
        str_ = '\t' + '\t'.join(map(str, aves))
        f.write(str_ + '\n')
         

dt     = 0.10
nstart = 1
nend   = 16
ncore  = 752
n_monomer = 15
tot_frame = 40
nsep   = 10
tname  = '../pdbdir/em.pdb'

fname  = 'MD/md001_sep.xtc'
honame = 'DAT/hobond_sep.dat'
ohname = 'DAT/ohbond_sep.dat'
hoselfname = 'DAT/hobondself_sep.dat'
ohselfname = 'DAT/ohbondself_sep.dat'

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
str_Hg = '# time\t'+'\t'.join(map(str, atom_indicesH_g))
with open(honame, 'wt') as f:
    f.write(str_Hg+'\n')
str_Og = '# time\t'+'\t'.join(map(str, atom_indicesO_g))
with open(ohname, 'wt') as f:
    f.write(str_Og+'\n')
with open(hoselfname, 'wt') as f:
    f.write(str_Hg+'\n')
str_Og = '# time\t'+'\t'.join(map(str, atom_indicesO_g))
with open(ohselfname, 'wt') as f:
    f.write(str_Og+'\n')


atom_indicesH_c = [ index for index in atom_indicesH if not index in ghost and index in core]
atom_indicesO_c = [ index for index in atom_indicesO if not index in ghost and index in core]

atom_indicesH_v = [ index for index in atom_indicesH if not index in ghost and not index in core]
atom_indicesO_v = [ index for index in atom_indicesO if not index in ghost and not index in core]

atom_pairsHO = topology.select_pairs(atom_indicesH_g, atom_indicesO_v) 
atom_pairsOH = topology.select_pairs(atom_indicesO_g, atom_indicesH_v)

atom_pairsHO_self = topology.select_pairs(atom_indicesH_g, atom_indicesO_c) 
atom_pairsOH_self = topology.select_pairs(atom_indicesO_g, atom_indicesH_c)

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

    hbondHO = calcHbond(rHO, atom_indicesH_g, atom_pairsHO, frame, t, dt, 'HO')
    writeHbond(hbondHO, honame) 
    hbondOH = calcHbond(rOH, atom_indicesO_g, atom_pairsOH, frame, t, dt, 'OH')
    writeHbond(hbondOH, ohname) 

    hbondHOself = calcHbond(rHO_self, atom_indicesH_g, atom_pairsHO_self, frame, t, dt, 'HO', selfflag=True)
    writeHbond(hbondHOself, hoselfname)

    hbondOHself = calcHbond(rOH_self, atom_indicesO_g, atom_pairsOH_self, frame, t, dt, 'OH', selfflag=True)
    writeHbond(hbondOHself, ohselfname)

for fname in [honame, ohname, hoselfname, ohselfname]:
    calc_ave(fname)
