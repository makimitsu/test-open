from scipy.interpolate import Rbf

def bz_rbfinterp(rpos, zpos, grid2D, bz, ok, t):
    smoothval = 0.05
    func = 'multiquadric'
    const = 0.015
    
    x = zpos(ok)
    y = rpos(ok)
    z = bz[t, ok].astype(float)
    
    # Create the RBF interpolator
    rbf_interpolator = Rbf(x, y, z, function=func, epsilon=const, smooth=smoothval)

    # Interpolation points
    zq_flat = grid2D['zq'].ravel()
    rq_flat = grid2D['rq'].ravel()

    # Perform the interpolation
    vq_flat = rbf_interpolator(zq_flat, rq_flat)

    # Reshape the result to match the grid shape
    vq = vq_flat.reshape(grid2D['zq'].shape)