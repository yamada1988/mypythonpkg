import numpy as np

# parameters
pi = np.pi
tN =  512**2
dt =  0.00010
tmax = tN*dt

# create discretized t and w
ws = 2.0*pi/dt   # sampling frequency
dw = 2.0*pi/tN # delta omega 
t = np.arange(0, tmax, dt)
w = np.fft.fftfreq(tN, dt/(2.0*pi))
# Set data
gamma = 2.0e0 # relaxation constant
A = np.exp(-1.0e0*gamma*np.abs(t)) # exp(-g|t|)

#print('Check A...')
#for i in range(100):
#    print(i,A[i])

# FFT calculation
Aw = np.fft.fft(A, tN) *dt

print(Aw)
# Lorentzian function
Bw = 1.0/(1.0e0j*w+gamma)

# IFFT calculation
iAw = np.fft.ifft(Aw, tN) /dt

# shift frequency
#w_ = np.fft.ifftshift(w)
#print(w_)

#Cw = 2.0/(w_**2+gamma**2)
#iCw = np.fft.ifft(Bw,tN)/dt

print('Check FFT...')
with open('Aw.dat', 'wt') as f:
    for i in range(100):
        f.write('{0:8.5f}\t{1:8.5f}\t{2:8.5f}\t{3:8.5f}\t{4:8.5f}\n'.format(w[i],Aw[i].real,Aw[i].imag, Bw[i].real, Bw[i].imag))

print('Check inverse-inverse translation...')
print('time\tA(t)\tIFFT(A)(t)\tdiff')
for i in range(100): 
    print(i,A[i],iAw[i].real, A[i]-iAw[i].real)
