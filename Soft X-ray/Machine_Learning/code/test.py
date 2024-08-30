import scipy.io
import numpy as np

def FindCircle(L):
    R = np.zeros((2*L))
    for i in range(2*L):
        for j in range(2*L):
            R[i, j] = np.sqrt((L-i+0.5)**2 + (j-L-0.5)**2)
    k = np.where(R < L)
    return k





# Check the shape to confirm it is 1976x8281
print(gm2d1.shape)  # Should output (1976, 8281)