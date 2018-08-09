from joblib import Parallel, delayed
from multiprocessing import Pool, Process
import threading
import numpy as np
from myopenmm import *
import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys
import os

def mdrun(obj, *args, **kwargs):
    print('mdrun is global function')
    return obj.mdrun(*args, **kwargs)

def wrapper(*args):
    return mdrun(*args)

args = sys.argv
stage = args[1]
index = args[2]

gpu_indexs = get_gpu_info()[1:]
print(gpu_indexs)
N = 2 
gpus = gpu_indexs*N

mds = []
sims = []
for i,index in enumerate(['001', '002']):
    platform = {'Precision': 'single', 'DeviceIndex': gpus[i]}
    print('i:', i, 'platform:', platform)

    md_nvt = MDConductor()
    inpf = '../inpdir/' + stage + '/nvt' + index + '.inp'
    md_nvt.loadFile(inpf)
    sysgro, systop = md_nvt.preparation(sysdir='SYS/', mddir='MD')
    simulation_0 = md_nvt.setup(sysgro, systop, deviceindex=gpus[i])

    # EM
    simulation_em = md_nvt.minimize(simulation_0, 'em', index,  mddir='MD/')
    mds.append(md_nvt)
    sims.append(simulation_em)

# NVT
# serialize
#for i,md_ in enumerate(mds) :
#    index = '{0:03d}'.format(i+1)
#    nvtlog = 'MD/nvt' + index + '.log'
#    mdrun(md_, sims[i], 'nvt', index, mddir='MD/')


# check
#mdrun(mds[0], sims[0], 'nvt', '001')


# parallel
arglist = []
for i, md_ in enumerate(mds):
    index = '{0:03d}'.format(i+1)
    args = [md_, sims[i], 'nvt', index]
    arglist.append(args)

# NG (pickle.PicklingError: Could not pickle the task to send it to the workers.)
# because md_, sim_[i] is SwigPyObject.
#callback = Parallel(n_jobs=-1, verbose=1)( [delayed(mdrun)(i) for i in arglist] )

# NG (TypeError: can't pickle SwigPyObject objects)
# because md_, sim_[i] is SwigPyObject.
#p = Pool(2)
#output = p.map(wrapper, arglist)
#p.close()

# NG (Exception: Error invoking kernel: CUDA_ERROR_NOT_INITIALIZED (3))
# CUDA_ERROR_NOT_INITIALIZED = 3:
# This indicates that the CUDA driver has not been initialized with cuInit() or that initialization has failed.

#jobs = []
#for args in arglist:
#    print('args:', args[0].pltform)
#    job = Process(target=mdrun, args=(args[0], args[1], args[2], args[3]))
#    jobs.append(job)
#    job.start()

#[job.join() for job in jobs]

#print('Finish')

T = []
for i,args in enumerate(arglist):
    print(args)
    print('thread-id:',i)
    T.append(threading.Thread(target=mdrun, args=(args[0], args[1], args[2], args[3])))

for i,t in enumerate(T):
    print('thread-id:', i)
    t.start()

for t in T:
    t.join()
