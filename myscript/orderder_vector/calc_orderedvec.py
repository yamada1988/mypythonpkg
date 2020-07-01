import numpy as np
import mdtraj as md
from scipy.spatial import distance_matrix
import sys

groname = 'PVA-20000.gro'

gro = md.load(groname)
L = gro.unitcell_vectors[0,0,0]
top = gro.top
selection = top.select('name c3')
groc = gro.atom_slice(selection)

Nchain = groc.n_residues
Natoms = groc.n_atoms
nc = int(Natoms/Nchain)

# rcs(t=0)
rcs = groc.xyz[0]
# rc[ichain, jatom, axis]
rc = rcs.reshape(Nchain,nc,3)
print(rc[0])
# rc_2k-1 (k=1,2,\cdots,n)
rc_odd = rc[:,::2,:]
print(rc_odd[0,:3])

# ordered vector vc_i,jatom = rc_ichain,jatom+2 - rc_ichain,jatom
vc_odd_0 = np.diff(rc_odd, axis=1)
# recalculate the ordered vector to consider periodic bounrady condition
vc_odd = vc_odd_0 - L * np.round(vc_odd_0/L)

natm_odd = len(rc_odd[0])
# reshape rc_odd array to calculate distance matrix (dim=2, rc_odd_re(X,Y))
rc_odd_re = rc_odd.reshape(Nchain*natm_odd,3)
# distance matrix for k-l pair, k = ik chain at jk monomer, l = il chain at jl monomer
drs_0 = np.sqrt(np.sum(((rc_odd_re[:,np.newaxis] - rc_odd_re[np.newaxis]) - L * np.round((rc_odd_re[:,np.newaxis] - rc_odd_re[np.newaxis])/L) )**2, axis=2))

# drs[ik,jk,il,jl] = distance between ik chain at jk monomer and il chain at jl monomer
drs = drs_0.reshape(Nchain,natm_odd,Nchain,natm_odd)
print(drs[0,0,0])
