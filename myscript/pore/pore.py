import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import time

class System:
    def __init__(self, dl, Nx, Ny, Nz, box, xl, yl, zl, pbc):
        self.cells = np.array([Cell([i,j,k], [dl,dl,dl], [Nx,Ny,Nz], box, x0=xl, y0=yl, z0=zl, pbc=pbc) for k in range(Nz) for j in range(Ny) for i in range(Nx)])
        self.dl = dl
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.origin = np.array([xl, yl, zl])


    def make_intlist(self, posinfo):
        self.posinfo = posinfo
        self.intlist_dict = {}
        for pname, pos_ in self.posinfo.items():
            b = np.trunc((pos_ - self.origin)/self.dl)
            self.intlist_dict[pname] = np.array([[[np.where((b[:,0] == i) & (b[:,1] == j) & (b[:,2] == k))[0] 
                                               for k in range(self.Nz)] for j in range(self.Ny)] for i in range(self.Nx)])
    
    def load_intlist(self):
        [c.load_intlist(self.intlist_dict) for c in self.cells]

    def calc_occupiedself):
        [c.gen_subcell() for c in self.cells]
        self.subposes = np.array([np.array(c.subpos) for c in self.cells])
        [c.calc_pore(self.posinfo) for c in self.cells]
        self.occupied_indexes = [zip(*c.occupied_index) for c in self.cells]


    def get_cellinfo(self):
        for i,c in enumerate(self.cells):
            print(i, c)
            print('volume:', c.volume)
            print('index:', c.index)
            print('pos:', c.pos)
            print('neighbor index3d:', c.neighbor_index3d)
            print('neighbor index:', c.neighbor_index)


class Cell:
    vdw_table = {'c3': 0.170, 'hc':0.120, 'oh':0.152}
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
        self.pore_index = np.array([[i,j,k] for k in range(self.nzlim) for j in range(self.nylim) for i in range(self.nxlim)])
        self.ds = float(self.size[0]) / ncell
        x0, y0, z0 = self.pos[0], self.pos[1], self.pos[2]
        self.subpos = np.array([[[[x0+ix*ds, y0+iy*ds, z0+iz*ds] for iz in range(ncell)] for iy in range(ncell)] for ix in range(ncell)], dtype=np.float32)


    def load_intlist(self, ilist_dict):
        nl = [inb for inb in self.neighbor_index3d]
        #print(nl)
        self.intlist_dict = {}
        for pname, ilist in ilist_dict.items():
            self.intlist_dict[pname] = []
            for n in nl:
                #print(n)
                #print([i for i in ilist_dict[pname][n[0], n[1], n[2]]])
                self.intlist_dict[pname] += [i for i in ilist[n[0], n[1], n[2]]]
        #print(self.id, self.intlist_dict)

    def get_uniques(self, seqs):
        seen = []
        for seq in seqs:
            for x in seq:
                if x not in seen:
                    seen.append(x)
        return seen

    def calc_occupied(self, posinfo):
        self.occupied_index = []
        print(self.id)
        #print(self.intlist_dict)
        ind_d = {}
        for pname, pos in posinfo.items():
            atm_pos = pos[self.intlist_dict[pname]]
            d_pos = self.subpos[:,:,:,np.newaxis] - atm_pos[np.newaxis,:]
            d_pos -= self.box*np.trunc(d_pos/(self.box*0.50e0))
            d = np.sqrt(np.sum(d_pos**2, axis=4))
            vdw_test = self.ds * 0.50
            datm = vdw_test + vdw_table[pname]
            ind_d[pname] = np.where(d < datm)[0:3]
            ind_d[pname] = list(set(list(zip(*ind_d[pname]))))
            #print(pname, ind_d[pname])
            self.occupied_index.append(ind_d[pname])
        #print(self.occupied_index)
        self.occupied_index = self.get_uniques(self.occupied_index)
        #print(len(self.occupied_index))


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
posinfo = {'c3': pos_c3, 'hc': pos_hc}

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

Nx = Ny = Nz = 4
Cells = np.array([Cell([i,j,k], [dl,dl,dl], [Nx,Ny,Nz], box, z0=zl, pbc='xy') for k in range(Nz) for j in range(Ny) for i in range(Nx)])
system = System(dl, Nx, Ny, Nz, box, 0.0, 0.0, zl, 'xy')
system.make_intlist(posinfo)
#print(system.intlist_dict['c3'])
system.load_intlist()
system.calc_pores()
for i,pi in enumerate(system.occupied_indexes):
    print(i, len(pi))
for i,ss in enumerate(system.subposes):
    print(ss[system.occupied_indexes[i]])

sys.exit()
