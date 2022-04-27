import numpy as np
import mdtraj as md
from scipy.spatial import distance_matrix
import sys

def distance_matrix(rvec, lbox, n1, n2):
    # distance matrix for k-l pair, k = ik chain at jk monomer, l = il chain at jl monomer
    d_0 = np.sqrt(np.sum(((rvec[:,np.newaxis] - rvec[np.newaxis]) - lbox * np.round((rvec[:,np.newaxis] - rvec[np.newaxis])/lbox) )**2, axis=2))
    # drs[ik,jk,il,jl] = distance between ik chain at jk monomer and il chain at jl monomer
    d = d_0.reshape(n1,n2,n1,n2)
    return d_0, d

def make_neighborlist(dists_0, r_threshold, n_c, n_atm):
    nlist = [np.where(dists_0[i] <= r_threshold) for i in range(n_c*n_atm)]
    return nlist 

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    vector = np.array(vector)
    if np.linalg.norm(vector) <= 0.00010:
        normv = 1.0
    else:
        normv = np.linalg.norm(vector)
    return vector / normv

def calc_p2(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    p2 = 0.50*(3.0E0*np.dot(v1_u,v2_u)**2.0 - 1.0)
    return p2


# parameters
rd0 = 2.0
rd = 1.0

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
# drs[ik,jk,il,jl] = distance between ik chain at jk monomer and il chain at jl monomer
drs_0, drs = distance_matrix(rc_odd_re, L, Nchain, natm_odd)

# calculate distance between k-l pair and make neighborlist
neighborlist = make_neighborlist(drs_0, rd0, Nchain, natm_odd)

P2 = np.zeros((Nchain,natm_odd))
nk = np.zeros((Nchain,natm_odd))
x = Nchain*natm_odd
for i in range(x):
    nk = 0
    ik = i // natm_odd
    jk = i % natm_odd
    if (jk == natm_odd-1):
        continue # vc_odd[natm_odd-1] does not exist.
    for l in neighborlist[i][0]:
        if i == l:
            continue
        il = l // natm_odd
        jl = l % natm_odd
        if (jl == natm_odd-1):
            continue # vc_odd[natm_odd-1] does not exist.
        #print("index2=",l, il, jl)
        dr = rc_odd[ik,jk] - rc_odd[il,jl]
        dr -= L * np.round(dr/L)
        dr = np.sqrt(np.sum(dr**2))
        #print(dr)
        if (dr < rd ):
            nk[ik, jk] += 1
            vi = vc_odd[ik,jk] 
            vl = vc_odd[il,jl]
            P2[ik,jk] += calc_p2(vi, vl)
            #print("p2:", calc_p2(vi, vl))
    P2[ik, jk] /= nk[ik, jk]
    print("index:", ik, jk, "P2:", P2[ik, jk], nk[ik, jk])
drs_0 = np.sqrt(np.sum(((rc_odd_re[:,np.newaxis] - rc_odd_re[np.newaxis]) - L * np.round((rc_odd_re[:,np.newaxis] - rc_odd_re[np.newaxis])/L) )**2, axis=2))

# drs[ik,jk,il,jl] = distance between ik chain at jk monomer and il chain at jl monomer
drs = drs_0.reshape(Nchain,natm_odd,Nchain,natm_odd)
print(drs[0,0,0])
