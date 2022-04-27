from pytrr import GroTrrReader
import sys
import numpy as np

args = sys.argv
fname = args[1]

# Mass Parameters
mO = 16.00
mH =  1.08

N = 3000
Rs = []
with GroTrrReader(fname) as trrfile:
    for frame in trrfile:
        print(frame['step'])
        frame_data = trrfile.get_data()
        for i in range(1, 3000+1):
            iO  = 3*i - 3
            iH1 = 3*i - 2
            iH2 = 3*i - 1
            xO  = np.array(frame_data['x'][iO])
            xH1 = np.array(frame_data['x'][iH1])
            xH2 = np.array(frame_data['x'][iH2])
#            for x in [xO, xH1, xH2]:
#                print(x)
            M = mO + mH + mH
            R = (mO * xO + mH * xH1 + mH * xH2)/ M
            print(i, R)
