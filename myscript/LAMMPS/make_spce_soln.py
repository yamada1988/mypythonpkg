import numpy as np
import sys

#
# make ethanol + water .data file for LAMMPS
#

atom_dic = {'C1':-0.0971, 'C2':0.13140,'O':-0.6028,'HC':0.04437, 'HT':0.01870, 'HO':0.39800,'HW':0.424, 'OW':-0.848}

args = sys.argv
enum = int(args[1])
wnum = int(args[2])
mkboxname = args[3]

atomnum_slt = 9
skipnum = atomnum_slt * enum
bondnum_slt = 8 * enum
anglenum_slt = 13 * enum
dihedralnum_slt = 12 * enum
impropnum_slt = 0 * enum
skipindex = 9 * enum

atomtype_slt = 5
bondtype_slt = 5
angletype_slt = 7
dihedraltype_slt = 4
improptype_slt = 0

ofname = 'test.data'

#header:nums

with open(ofname, 'wt') as of:
    of.write('\n')
    atomnum_water = wnum * 3
    bondnum_water = wnum * 2
    anglenum_water = wnum
    dihedralnum_water = 0
    impropnum_water = 0

    atomnum = atomnum_slt + atomnum_water
    bondnum = bondnum_slt + bondnum_water
    anglenum = anglenum_slt + anglenum_water
    dihedralnum = dihedralnum_slt + dihedralnum_water
    impropnum = impropnum_slt + impropnum_water

    l = '\t{0:5d} atoms\n\t{1:5d} bonds\n\t{2:5d} angles\n\t{3:5d} dihedrals\n\t{4:5d} impropers\n'.format(atomnum, bondnum, anglenum, dihedralnum, impropnum)
    of.write(l)

#header:types
    of.write('\n')
    atomtype_water = 2
    bondtype_water = 2
    angletype_water = 1
    dihedraltype_water = 0
    improptype_water = 0

    atomtype = atomtype_slt + atomtype_water
    bondtype = bondtype_slt + bondtype_water
    angletype = angletype_slt + angletype_water
    dihedraltype = dihedraltype_slt + dihedraltype_water
    improptype = improptype_slt + improptype_water

    l = '\t{0:5d} atom types\n\t{1:5d} bond types\n\t{2:5d} angle types\n\t{3:5d} dihedral types\n\t{4:5d} improper types\n'.format(atomtype, bondtype, angletype, dihedraltype, improptype)
    of.write(l)

#box:
    with open(mkboxname,'rt') as mkf:
        lines = [line.split() for line in mkf]
    boxline = lines[-1]

    of.write('\n')
    xlo = 0.00000
    xhi = float(boxline[0])*10.0
    ylo = 0.00000
    yhi = float(boxline[1])*10.0
    zlo = 0.0000
    zhi = float(boxline[2])*10.0
    l = '{0:11.8f}\t{1:11.8f} xlo xhi\n'.format(xlo, xhi)
    of.write(l)
    l = '{0:11.8f}\t{1:11.8f} ylo yhi\n'.format(ylo, yhi)
    of.write(l)
    l = '{0:11.8f}\t{1:11.8f} zlo zhi\n'.format(zlo, zhi)
    of.write(l)

#masses:ethanol
#1 12.010000
#2 1.008000
#3 16.000000
#4 1.008000
#5 1.008000
    of.write('\n Masses\n\n')
    l1 = ' 1\t12.01000\n'
    l2 = l1 + ' 2\t1.00800\n'
    l3 = l2 + ' 3\t16.0000\n'
    l4 = l3 + ' 4 1.008000\n'
    l5 = l4 + ' 5 1.008000\n'
    of.write(l5)

#masses:water
    owindex = atomtype_slt + 1
    hwindex = owindex + 1
    owmass = 15.9994
    hwmass = 1.008
    l = '{0:2d}\t{1:10.9f} # OW\n'.format(owindex, owmass)
    of.write(l)
    l = '{0:2d}\t{1:10.9f} # HW\n'.format(hwindex, hwmass)
    of.write(l)

#pairtype:ethanol
#1      0.109400                3.399670
#2      0.000000                0.000000
#3      0.210400                3.066473
#4      0.015700                2.471353
#5      0.015700                2.649533      
    of.write('\n Pair Coeffs # lj/cut/coul/long/omp\n\n')
    l1 = ' 1\t0.109400\t3.399670\n'
    l2 = l1 + ' 2\t0.000000\t0.000000\n'
    l3 = l2 + ' 3\t0.210400\t3.066473\n'
    l4 = l3 + ' 4\t0.015700\t2.471353\n'
    l5 = l4 + ' 5\t0.015700\t2.649533\n'
    of.write(l5)

#pairtype:water
    l = '\t{0:2d}\t0.155504063097514\t3.1655700000\t #OW,spce\n'.format(owindex)
    of.write(l)
    l = '\t{0:2d}\t0.000000000000000\t0.0000000000\t #HW,spce\n'.format(hwindex)
    of.write(l)

#bondtype:ethanol
#1 303.100000 1.535000
#2 314.100000 1.426000
#3 337.300000 1.092000
#4 335.900000 1.093000
#5 369.600000 0.974000
    of.write('\n Bond Coeffs # spce\n\n')
    l1 = ' 1\t303.100000\t1.535000\n'
    l2 = l1 + ' 2\t314.100000\t1.426000\n'
    l3 = l2 + ' 3\t337.300000\t1.092000\n'
    l4 = l3 + ' 4\t335.900000\t1.093000\n'
    l5 = l4 + ' 5\t369.600000\t0.974000\n'
    of.write(l5)

#bondtype:water
    oh1index = bondtype_slt + 1
    oh2index = oh1index + 1
    l = '\t{0:2d}\t0.0000000\t1.0000000000\t #OW,spce\n'.format(oh1index)
    of.write(l)
    l = '\t{0:2d}\t0.0000000\t1.0000000000\t #HW,spce\n'.format(oh2index)
    of.write(l)

#angletype:ethanol
#1 46.400000 110.050000
#2 39.400000 108.350000
#3 67.700000 109.430000
#4 46.400000 110.070000
#5 51.000000 109.880000
#6 39.200000 109.550000
#7 47.100000 108.160000
    of.write('\n Angle Coeffs # spce\n\n')
    l1 = ' 1\t46.400000\t110.050000\n'
    l2 = l1 + ' 2\t39.400000\t108.350000\n'
    l3 = l2 + ' 3\t67.700000\t109.430000\n'
    l4 = l3 + ' 4\t46.400000\t110.070000\n'
    l5 = l4 + ' 5\t51.000000\t109.880000\n'
    l6 = l5 + ' 6\t39.200000\t109.550000\n'
    l7 = l6 + ' 7\t47.100000\t108.160000\n'
    of.write(l7)

#angletype:water
    hohindex = angletype_slt + 1
    l = '\t{0:2d}\t0.0000000\t109.47000000\t #OW,spce\n'.format(hohindex)
    of.write(l)

#dihedraltype:ethanol
#1  0.155556 -0.466667 0.000000 0.622222 0.000000
#2  0.250000 0.250000 0.000000 0.000000 0.000000
#3  0.166667 -0.500000 0.000000 0.666667 0.000000
#4  0.410000 -0.230000 0.000000 0.640000 0.000000
    of.write('\n Dihedral Coeffs # solute\n\n')
    l1 = ' 1\t0.155556\t-0.466667\t0.000000\t0.622222\t0.000000\n'
    l2 = l1 + ' 2\t0.250000\t0.250000\t0.000000\t0.000000\t0.000000\n'
    l3 = l2 + ' 3\t0.166667\t-0.500000\t0.000000\t0.666667\t0.000000\n'
    l4 = l3 + ' 4\t0.410000\t-0.230000\t0.000000\t0.640000\t0.000000\n'
    of.write(l4)

#atoms:
    of.write('\n Atoms \n\n')
    with open(mkboxname,'rt') as mkf:
        lines = [line for line in mkf]
        atomlines = lines[2:-1]
        for al in atomlines:
            #print(al)
            molindex = int(al[0:5])
            moltype = str(al[5:10])
            atomtype = str(al[10:15].split()[0])
            index = int(al[15:20])
            #print(molindex, moltype, atomtype, index)
            if atomtype == 'C1':
                atomtype_lammps = 1
            elif atomtype == 'C2':
                atomtype_lammps = 1
            elif atomtype == 'O':
                 atomtype_lammps = 3
            elif atomtype == 'HC':
                 atomtype_lammps = 5
            elif atomtype == 'HT':
                 atomtype_lammps = 4
            elif atomtype == 'HO':
                 atomtype_lammps = 2
            elif atomtype == 'OW':
                 atomtype_lammps = 6
            elif atomtype == 'HW':
                 atomtype_lammps = 7
            else:
                atomtype_lammps = 0
            charge = atom_dic[atomtype]
            x = float(al[20:28])*10.0 #nm => A
            y = float(al[28:36])*10.0 #nm => A
            z = float(al[36:44])*10.0 #nm => A
            l = '{0:4d}\t{1:4d}\t{2:2d}\t{3:5.4f}\t{4:8.5f}\t{5:8.5f}\t{6:8.5f}\n'.format(index,molindex,atomtype_lammps,charge,x,y,z)
            of.write(l)

#Bonds:ethanol
    of.write('\n Bonds\n\n')
    eindex = 1
    for n_e in range(1,enum+1):
        c1index = 9*(n_e-1)+1
        c2index = c1index + 1
        oh1index = c2index + 1
        hc1index = oh1index + 1
        hc2index = hc1index + 1
        h11index = hc2index + 1
        ho1index = h11index + 1
        hc3index = ho1index + 1
        h12index = hc3index + 1
#12334534
        l1 = '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\n'.format(eindex,1,c1index,c2index)
        eindex += 1
        l2 = l1 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\n'.format(eindex,2,c2index,oh1index)
        eindex += 1
        l3 = l2 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\n'.format(eindex,3,hc1index,c1index)
        eindex += 1
        l4 = l3 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\n'.format(eindex,3,c1index,hc2index)
        eindex += 1 
        l5 = l4 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\n'.format(eindex,4,oh1index,c2index)
        eindex += 1
        l6 = l5 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\n'.format(eindex,5,oh1index,ho1index)
        eindex += 1
        l7 = l6 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\n'.format(eindex,3,hc3index,c1index)
        eindex += 1
        l8 = l7 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\n'.format(eindex,4,h12index,c2index)
        eindex += 1
        of.write(l8)

#Bonds:water
    bindex = bondnum_slt + 1
    for n_w in range(1,wnum+1):
        oindex = 3*(n_w-1)+1+skipindex
        for itype in [1,2]:
            hindex = oindex + itype
            l = '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\n'.format(bindex,itype+bondtype_slt,oindex,hindex)
            of.write(l)
            bindex += 1

#Angles:ethanol
#1 1 2 1 4
#2 1 2 1 5
#3 1 2 1 8
#4 2 4 1 5
#5 2 4 1 8
#6 2 5 1 8
#7 3 1 2 3
#8 4 1 2 6
#9 4 1 2 9
#10 5 3 2 6
#11 5 3 2 9
#12 6 6 2 9
#13 7 2 3 7
    of.write('\n Angles\n\n')
    eindex = 1
    for n_e in range(1,enum+1):
        c1index = 9*(n_e-1)+1
        c2index = c1index + 1
        oh1index = c2index + 1
        hc1index = oh1index + 1
        hc2index = hc1index + 1
        h11index = hc2index + 1
        ho1index = h11index + 1
        hc3index = ho1index + 1
        h12index = hc3index + 1
        l1 = '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,1,c2index,c1index,hc1index)
        eindex += 1
        l2 = l1 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,1,c2index,c1index,hc2index)
        eindex += 1
        l3 = l2 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,1,c2index,c1index,hc3index)
        eindex += 1
        l4 = l3 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,2,hc1index,c1index,hc2index)
        eindex += 1 
        l5 = l4 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,2,hc1index,c1index,hc3index)
        eindex += 1
        l6 = l5 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,2,hc2index,c1index,hc3index)
        eindex += 1
        l7 = l6 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,3,c1index,c2index,oh1index)
        eindex += 1
        l8 = l7 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,4,c1index,c2index,h11index)
        eindex += 1
        l9 = l8 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,4,c1index,c2index,h12index)
        eindex += 1
        l10 = l9 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,5,oh1index,c2index,h11index)
        eindex += 1
        l11 = l10 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,5,oh1index,c2index,h12index)
        eindex += 1
        l12 = l11 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,6,h11index,c2index,h12index)
        eindex += 1
        l13 = l12 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,7,c2index,oh1index,ho1index)
        eindex += 1
        of.write(l13)

#Angles:water
    aindex = anglenum_slt + 1
    for n_w in range(1,wnum+1):
        oindex = 3*(n_w-1)+1+skipindex
        itype = 1
        hindex_1 = oindex + 1
        hindex_2 = oindex + 2
        l = '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(aindex,itype+angletype_slt,hindex_1,oindex,hindex_2)
        of.write(l)
        aindex += 1


#Dihedral:ethanol
#1 2 4 1 2 3
#2 1 4 1 2 6
#3 1 4 1 2 9
#4 2 5 1 2 3
#5 1 5 1 2 6
#6 1 5 1 2 9
#7 2 8 1 2 3
#8 1 8 1 2 6
#9 1 8 1 2 9
#10 4 1 2 3 7
#11 3 6 2 3 7
#12 3 9 2 3 7
    of.write('\n Dihedrals\n\n')
    eindex = 1
    for n_e in range(1,enum+1):
        c1index = 9*(n_e-1)+1   #1
        c2index = c1index + 1   #2
        oh1index = c2index + 1  #3
        hc1index = oh1index + 1 #4
        hc2index = hc1index + 1 #5
        h11index = hc2index + 1 #6
        ho1index = h11index + 1 #7
        hc3index = ho1index + 1 #8
        h12index = hc3index + 1 #9
        l1 = '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\t{5:4d}\n'.format(eindex,       2,hc1index,c1index,c2index,oh1index)
        eindex += 1
        l2 = l1 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,  1,hc1index,c1index,c2index,h11index)
        eindex += 1
        l3 = l2 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,  1,hc1index,c1index,c2index,h12index)
        eindex += 1
        l4 = l3 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,  2,hc2index,c1index,c2index,oh1index)
        eindex += 1 
        l5 = l4 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,  1,hc2index,c1index,c2index,h11index)
        eindex += 1
        l6 = l5 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,  1,hc3index,c1index,c2index,h11index)
        eindex += 1
        l7 = l6 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,  2,hc3index,c1index,c2index,oh1index)
        eindex += 1
        l8 = l7 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,  1,hc3index,c1index,c2index,h11index)
        eindex += 1
        l9 = l8 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,  1,hc3index,c1index,c2index,h12index)
        eindex += 1
        l10 = l9 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex, 4,c1index,c2index,oh1index,ho1index)
        eindex += 1
        l11 = l10 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,3,h11index,c2index,oh1index,ho1index)
        eindex += 1
        l12 = l11 + '\t{0:4d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:4d}\n'.format(eindex,3,h12index,c2index,oh1index,ho1index)
        eindex += 1
        of.write(l12)

    of.write('\n')
