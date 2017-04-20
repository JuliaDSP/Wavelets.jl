
# in python3 using PyWavelets
import numpy as np
import pywt

def to_one_array(cA, cD):
    m = cA.shape[0]
    n = cA.shape[1]
    A = np.zeros((2*m,2*n))
    A[0:m,0:n] = cA
    A[m:2*m,0:n] = -cD[0]
    A[0:m,n:2*n] = -cD[1]
    A[m:2*m,n:2*n] = cD[2]
    return A

n1 = 4
n2 = 8

np.random.seed(123)
x = np.random.randn(n1,n2)
np.savetxt('filter2d_nonsquare_data.txt', x)

wavelets = ['haar']
wnames = ['Haar0']

for j in range(0,len(wavelets)):
    xA, xD = pywt.dwt2(x, wavelets[j], mode='periodic')
    A = to_one_array(xA, xD)
    filename = 'filter2d_nonsquare_' + wnames[j] + '.txt'
    np.savetxt(filename, A)
