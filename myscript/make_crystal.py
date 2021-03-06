import sys
import numpy as np

inpf = 'PE024.gro'
with open(inpf, 'rt') as f:
    lines = [line.split() for line in f]
Nline = len(lines)
natom = int(lines[1][0])
Lbox = list(map(float, lines[-1]))
zlength = float(lines[-2][5]) - float(lines[2][5])

rotf = 'PE024_rot.gro'
with open(rotf, 'rt') as f:
    rotlines = [line.split() for line in f]


newlines = []
line0 = '    1PVA24   '
for line in lines[2:-1]:
    if line[1] == 'C':
        line[1] = 'c3'
    elif line[1] == 'H':
        line[1] = 'hc'
    #newline = line0 + line[1] 
    newline = line[1]
    newlines.append(newline)

print(newlines)

rvec = ['' for i in range(Nline-2)]
for i in range(2, Nline-1):
    print(i, lines[i])
    #print(np.array(list(map(float, [lines[i][3], lines[i][4], lines[i][5]]))))
    rvec[i-2] = np.array(list(map(float, [lines[i][3], lines[i][4], lines[i][5]])))
    #print(rvec[i])
#print(rvec)

for i in range(Nline-3):
   print(i, rvec[i])

hvec = np.array([0.370, 0.2465, 0.0])
rotrvec = ['' for i in range(Nline-2)]
for i in range(2, Nline-1):
    print(i, rotlines[i])
    rotrvec[i-2] = np.array(list(map(float, [rotlines[i][3], rotlines[i][4], rotlines[i][5]]))) + hvec


translate_vec = np.array([0.740, 0.493, zlength + 0.30])
nx, ny, nz = [6,6,1]
repnum = nx*ny*nz
centnum = (nx-1)*(ny-1)*nz
totalnum = repnum + centnum
newvec = np.zeros((totalnum, natom, 3))
for iz in range(nz):
    for iy in range(ny):
        for ix in range(nx):
            i = ny*nx*iz+ nx*iy + ix
            tvec = np.array([ix * translate_vec[0], iy * translate_vec[1], iz*translate_vec[2]])
            for j in range(natom):
            #print(rvec[j], tvec)
            #print(rvec[j]+tvec)
                newvec[i][j] = np.array(rvec[j]) + tvec
            #print(i, newvec[i])

icount = repnum - 1
for iz in range(nz):
    for iy in range(ny-1):
        for ix in range(nx-1):
            i = ny*nx*iz+ nx*iy + ix
            icount += 1
            print(i, icount)
            tvec = np.array([ix * translate_vec[0], iy * translate_vec[1], iz*translate_vec[2]])
            for j in range(natom):
                newvec[icount][j] = rotrvec[j] + tvec


outf = 'crystal_' + inpf
with open(outf, 'wt') as f:
    f.write('crystal structure generated by python\n')
    f.write('{0:7d}\n'.format(totalnum*natom))
    for nres in range(totalnum):
        for j in range(natom):
            natm = nres*natom + j+1
            l = '{0:5d}PVA24   '.format(nres+1) + newlines[j] + '{0:5d}{1:8.3f}{2:8.3f}{3:8.3f}\n'.format(natm,newvec[nres][j][0], newvec[nres][j][1], newvec[nres][j][2])
            #print(l)
            f.write(l)
 
    f.write('{0:5.3f} {1:5.3f} {2:5.3f}'.format((nx)*translate_vec[0]+0.20, (ny)*translate_vec[1]+0.20, (nz)*Lbox[2]+0.50))

cc_xy = np.array([rvec[4][0]-rvec[0][0], rvec[4][1]-rvec[0][1]])
ncc_xy = cc_xy / np.sqrt(cc_xy[0]**2+cc_xy[1]**2)
theta = np.arccos(ncc_xy[1])
print(ncc_xy, theta, np.rad2deg(theta))


#for i in range(Nline): 
#    print(lines[i])
#    print(newlines[i])
