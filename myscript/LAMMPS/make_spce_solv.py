import numpy as np
import sys

#
# make ethanol + water .data file for LAMMPS
#

atom_dic = {'C1':-0.0971, 'C2':0.13140,'O':-0.6028,'HC':0.04437, 'HT':0.01870, 'HO':0.39800,'HW':0.424, 'OW':-0.848}

args = sys.argv
wnum = int(args[1])
mkboxname = args[2]

enum = 0
atomnum_slt = 9 * enum
skipnum = atomnum_slt 
bondnum_slt = 8 * enum
anglenum_slt = 13 * enum
dihedralnum_slt = 12 * enum
impropnum_slt = 0 * enum
skipindex = 9 * enum

atomtype_slt = 0
bondtype_slt = 0
angletype_slt = 0
dihedraltype_slt = 0
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
    of.write('\n Masses\n\n')
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
    of.write('\n Pair Coeffs # lj/cut/coul/long/omp\n\n')
#pairtype:water
    l = '\t{0:2d}\t0.155504063097514\t3.1655700000\t #OW,spce\n'.format(owindex)
    of.write(l)
    l = '\t{0:2d}\t0.000000000000000\t0.0000000000\t #HW,spce\n'.format(hwindex)
    of.write(l)

#bondtype:ethanol
    of.write('\n Bond Coeffs # spce\n\n')
#bondtype:water
    oh1index = bondtype_slt + 1
    oh2index = oh1index + 1
    l = '\t{0:2d}\t0.0000000\t1.0000000000\t #OW,spce\n'.format(oh1index)
    of.write(l)
    l = '\t{0:2d}\t0.0000000\t1.0000000000\t #HW,spce\n'.format(oh2index)
    of.write(l)

#angletype:ethanol
    of.write('\n Angle Coeffs # spce\n\n')
#angletype:water
    hohindex = angletype_slt + 1
    l = '\t{0:2d}\t0.0000000\t109.47000000\t #OW,spce\n'.format(hohindex)
    of.write(l)

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
                 atomtype_lammps = 1
            elif atomtype == 'HW':
                 atomtype_lammps = 2
            else:
                atomtype_lammps = 0
            charge = atom_dic[atomtype]
            x = float(al[20:28])*10.0 #nm => A
            y = float(al[28:36])*10.0 #nm => A
            z = float(al[36:44])*10.0 #nm => A
            l = '{0:4d}\t{1:4d}\t{2:2d}\t{3:8.6f}\t{4:8.5f}\t{5:8.5f}\t{6:8.5f}\n'.format(index,molindex,atomtype_lammps,charge,x,y,z)
            of.write(l)

#Bonds:ethanol
    of.write('\n Bonds\n\n')
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
    of.write('\n Angles\n\n')
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

