import numpy as np

def elliptic_E(m):
    N = 200
    theta = np.linspace(0, np.pi / 2, N)
    d_theta = np.diff(theta, prepend=0)
    ellipE = np.zeros_like(m)
    rows, columns = m.shape
    for i in range(rows):
        for j in range(columns):
            integrand = np.sqrt(1 - m[i, j] * np.sin(theta)**2)
            ellipE[i, j] = np.sum(integrand * d_theta)
    return ellipE

def elliptic_K(m):
    N = 200
    theta = np.linspace(0, np.pi / 2, N)
    d_theta = np.diff(theta, prepend=0)  # Prepend 0 to align with MATLAB indexing
    ellipK = np.zeros_like(m)

    rows, columns = m.shape
    for i in range(rows):
        for j in range(columns):
            integrand = (1 - m[i, j] * np.sin(theta)**2)**(-0.5)
            ellipK[i, j] = np.sum(integrand * d_theta)

    return ellipK

def B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,probe_mesh_r,probe_mesh_z,Br_on):
    mu0 = 4*np.pi*1e-7
    alpha = probe_mesh_r/r_EF
    beta1 = (probe_mesh_z+z1_EF)/r_EF
    beta2 = (probe_mesh_z+z2_EF)/r_EF
    Q1 = (1+alpha)**2+beta1**2
    Q2 = (1+alpha)**2+beta2**2
    k1 = ((4*alpha)/Q1)**0.5
    k2 = ((4*alpha)/Q2)**0.5
    m1 = k1**2
    m2 = k2**2
    B0 = i_EF*n_EF*mu0/(2*r_EF)
    
    Bz1 = B0 / np.sqrt(Q1 * np.pi) * (elliptic_E(m1) * (1 - alpha**2 - beta1**2) / (Q1 - 4 * alpha) + elliptic_K(m1))
    Bz2 = B0 / np.sqrt(Q2 * np.pi) * (elliptic_E(m2) * (1 - alpha**2 - beta2**2) / (Q2 - 4 * alpha) + elliptic_K(m2))
    Bz_EF = Bz1 + Bz2
    
    if Br_on:
        gamma1 = (probe_mesh_z + z1_EF) / probe_mesh_r
        gamma2 = (probe_mesh_z + z2_EF) / probe_mesh_r
        Br1 = B0 * gamma1 / np.sqrt(Q1 * np.pi) * (elliptic_E(m1) * (1 + alpha**2 + beta1**2) / (Q1 - 4 * alpha) - elliptic_K(m1))
        Br2 = B0 * gamma2 / np.sqrt(Q2 * np.pi) * (elliptic_E(m2) * (1 + alpha**2 + beta2**2) / (Q2 - 4 * alpha) - elliptic_K(m2))
        Br_EF = Br1 + Br2
    else:
        Br_EF = np.zeros_like(Bz_EF)

