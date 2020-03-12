import mdtraj as md
sysgro = 'SYS/system0001.gro'
traj = 'MD/refs0011_01.xtc'
t = md.load(traj, top=sysgro)

N = 500
for i in range(N):
    print(t.openmm_positions(i)[0])
