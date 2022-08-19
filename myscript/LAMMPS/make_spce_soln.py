import numpy as np
import sys

#
# make ethanol + water .data file for LAMMPS
#

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
    of.write('\n')
    xlo = 0.00000
    xhi = 96.00000
    ylo = 0.00000
    yhi = 96.0000
    zlo = 0.0000
    zhi = 96.0000
    l = '{0:11.8f}\t{1:11.8f} xlo xhi\n'.format(xlo, xhi)
    of.write(l)
    l = '{0:11.8f}\t{1:11.8f} ylo yhi\n'.format(ylo, yhi)
    of.write(l)
    l = '{0:11.8f}\t{1:11.8f} zlo zhi\n'.format(zlo, zhi)
    of.write(l)

#masses:water
    of.write('\n Masses\n\n')
    owindex = atomtype_slt + 1
    hwindex = owindex + 1
    owmass = 15.9994
    hwmass = 1.008
    l = '{0:2d}\t{1:10.9f} # OW\n'.format(owindex, owmass)
    of.write(l)
    l = '{0:2d}\t{1:10.9f} # HW\n'.format(hwindex, hwmass)
    of.write(l)

#pairtype:water
    of.write('\n Pair Coeffs # lj/cut/coul/long/omp\n\n')
    l = '\t{0:2d}\t0.155504063097514\t3.1655700000\t #OW,spce\n'.format(owindex)
    of.write(l)
    l = '\t{0:2d}\t0.000000000000000\t0.0000000000\t #HW,spce\n'.format(hwindex)
    of.write(l)

#bondtype:water
    oh1index = 1
    oh2index = 2
    of.write('\n Bond Coeffs # spce\n\n')
    l = '\t{0:2d}\t0.0000000\t1.0000000000\t #OW,spce\n'.format(oh1index)
    of.write(l)
    l = '\t{0:2d}\t0.0000000\t1.0000000000\t #HW,spce\n'.format(oh2index)
    of.write(l)

#angletype:water
    hohindex = 1
    of.write('\n Angle Coeffs # spce\n\n')
    l = '\t{0:2d}\t0.0000000\t109.47000000\t #OW,spce\n'.format(hohindex)
    of.write(l)

#atoms:
    of.write('\n Atoms \n\n')
    with open(mkboxname,'rt') as mkf:
        lines = [line.split() for line in mkf]
        atomlines = lines[2:-1]
        for al in atomlines:
            #print(al)
            index = int(al[2])
            try:
                resnum = int(al[0].split('SOL')[0])
            except:
                resnum = int(al[0].split('MOL')[0])
            if 'OW' in al[1]:
                atomtype = owindex
                charge = -0.848 #spce
            elif 'HW' in al[1]:
                atomtype = hwindex
                charge = 0.424 #spce
            else:
                atomtype = 0
                charge = 0.0
            x = float(al[3])*10.0 #nm => A
            y = float(al[4])*10.0 #nm => A
            z = float(al[5])*10.0 #nm => A
            l = '{0:4d}\t{1:4d}\t{2:2d}\t{3:5.4f}\t{4:8.5f}\t{5:8.5f}\t{6:8.5f}\n'.format(index,resnum,atomtype,charge,x,y,z)
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
