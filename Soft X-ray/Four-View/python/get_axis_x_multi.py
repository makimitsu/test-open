import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

def get_axix_x_multi(grid2D, data2D):
    trange = data2D['trange']
    psi = data2D['psi']
    rq = grid2D['rq']
    zq = grid2D['zq']
    
    rqList = np.tile(rq[:, :, np.newaxis], (1, 1, len(trange)))
    zqList = np.tile(zq[:, :, np.newaxis], (1, 1, len(trange)))
    
    psiRidge = np.max(psi, axis=0)
    psiRidgeIdx = np.argmax(psi, axis=0)
    
    axisCandidate = np.array([find_peaks(psiRidge[:, i], distance=1, width=1)[0] for i in range(psiRidge.shape[1])])

    magAxisList_r = np.full((2, len(trange)), np.nan)
    magAxisList_z = np.full((2, len(trange)), np.nan)
    magAxisList_psi = np.full((2, len(trange)), np.nan)

    xPointList_r = np.full(len(trange), np.nan)
    xPointList_z = np.full(len(trange), np.nan)
    xPointList_psi = np.full(len(trange), np.nan)
    
    rLim = [np.min(rq), np.max(rq)]
    
    for i in range(len(trange)):
        psiRidge_t = psiRidge[:, :, i]
        psiRidgeIdx_t = psiRidgeIdx[:, :, i]
        axisCandidate_t = axisCandidate[:, :, i]
        
        if np.any(psiRidgeIdx_t[axisCandidate_t]):
            magAxisList_r[:, i] = rqList[psiRidgeIdx_t[axisCandidate_t]]
            magAxisList_z[:, i] = zqList[psiRidgeIdx_t[axisCandidate_t]]
            magAxisList_psi[:, i] = psi[psiRidgeIdx_t[axisCandidate_t]]
        
        smoothed_psiRidge_t = smooth(psiRidge_t)
        xpointIdx = find_peaks(-smoothed_psiRidge_t, distance=1)[0]  # find local minima



    