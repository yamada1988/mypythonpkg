import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import time


class Cell:
    def __init__(self, index, size, nlims, Lbox, x0=0.0, y0=0.0, z0=0.0, pbc='xyz'):
        self.id = index[2] * nlims[1]*nlims[0] + index[1]*nlims[0] + index[0] 
        self.ix, self.iy, self.iz = index[0], index[1], index[2]
        self.index = [self.ix, self.iy, self.iz]
        self.size = size
        self.volume = size[0]*size[1]*size[2]
        self.pos = np.array([x0+size[0]*index[0], y0+size[1]*index[1], z0+size[2]*index[2]])
        self.nxlim = nlims[0]
        self.nylim = nlims[1]
        self.nzlim = nlims[2]
        self.pbc = pbc
        if self.pbc == 'xyz':
            pbcx = [-1,0,1]
            pbcy = [-1,0,1]
            pbcz = [-1,0,1]
        elif pbc == 'xy':
            pbcx = [-1,0,1]
            pbcy = [-1,0,1]
            pbcz = [0]
        self.neighbor_index3d = [[(index[0]+i)%nlims[0], (index[1]+j)%nlims[1], (index[2]+k)%nlims[2]] for k in pbcz for j in pbcy for i in pbcx]
        self.neighbor_index = [n3d[2]*self.nylim*self.nxlim + n3d[1]*self.nxlim + n3d[0] for n3d in self.neighbor_index3d]
        self.intlist = []
        self.box = Lbox

    def gen_subcell(self, ncell=10):
        self.poreindex = []
        self.ds = float(self.size[0]) / ncell
        x0, y0, z0 = self.pos[0], self.pos[1], self.pos[2]
        self.subpos = np.array([[[[x0+ix*ds, y0+iy*ds, z0+iz*ds] for iz in range(ncell)] for iy in range(ncell)] for ix in range(ncell)], dtype=np.float32)


    def load_intlist(self, ilist):
        nl = [inb for inb in self.neighbor_index3d]
        for n in nl:
            #print(n)
            #print([i for i in ilist[n[0], n[1], n[2]]])
            self.intlist += [i for i in ilist[n[0], n[1], n[2]]]

    def calc_pore(self, pos, vdwradii):
        print(self.id)
        atm_pos = pos[self.intlist]
        d_pos = self.subpos[:,:,:,np.newaxis] - atm_pos[np.newaxis,:]
        d_pos -= self.box*np.trunc(d_pos/(self.box*0.50e0))
        dc = np.sqrt(np.sum(d_pos**2, axis=4))
        vdw_test = self.ds * 0.50
        drc = vdw_test + vdwradii
        ind_dc = np.where(dc < drc)[0:3]
        self.poreindex.append(list(set(list(zip(*ind_dc)))))
        


# unit:nm
vdw_table = {'c3': 0.170, 'hc':0.120, 'oh':0.152}
vdw_test = 0.050

fname = 'npt03_0001.gro'
t = md.load(fname)
top = t.topology
c3_inds = top.select("name c3")
hc_inds = top.select("name hc")
pos = t.xyz[0]
pos_c3 =  np.array(pos[c3_inds],dtype=np.float32)
pos_hc =  np.array(pos[hc_inds],dtype=np.float32)
box = t.unitcell_lengths[0]
Lx = box[0]
Ly = box[1]
zh = np.max(pos[:,2])
zl = np.min(pos[:,2])
Lz = zh-zl
ds = vdw_test*2.0 #nm
dl = 1.0 #nm
boxh = 0.50*box
Nx = int(Lx/1.0) + 1
Ny = int(Ly/1.0) + 1
Nz = int(Lz/1.0) + 1

Cells = np.array([Cell([i,j,k], [dl,dl,dl], [Nx,Ny,Nz], box, z0=zl, pbc='xy') for k in range(Nz) for j in range(Ny) for i in range(Nx)])
#for i,c in enumerate(Cells):
#    print(i, c)
#    print('volume:', c.volume)
#    print('index:', c.index)
#    print('pos:', c.pos)
#    print('neighbor index3d:', c.neighbor_index3d)
#    print('neighbor index:', c.neighbor_index)

def make_intlist(pos_, dl, nx, ny, nz, origin):
    b = np.trunc((pos_-origin)/dl)
    intlist = np.array([[[np.where((b[:,0] == i) & (b[:,1] == j) & (b[:,2] == k))[0] for k in range(nz)] for j in range(ny)] for i in range(nx)])
    return intlist

print('make_intlist:')
orig = np.array([0.0, 0.0, zl])
intlist = make_intlist(pos_c3, dl, Nx, Ny, Nz, orig)
#for i,c in enumerate(Cells):
#    if i == 0:
#        print(c.id, c.index)
#        print(intlist[c.ix, c.iy, c.iz])
 

# load intlist to each cell class
[c.load_intlist(intlist) for c in Cells]
#print(Cells[0].pos, Cells[0].intlist)
#print(pos_c3[Cells[0].intlist])

[c.gen_subcell() for c in Cells]

pores = [c.calc_pore(pos_c3, vdw_table['c3']) for c in Cells]
npore = 0
for i,c in enumerate(Cells):
    print(i, c.poreindex, len(c.poreindex))
    npore += len(c.poreindex)
print(npore)
sys.exit()
