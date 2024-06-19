import os
import pandas as pd
import numpy as np
import scipy.io as sio
from scipy.io import loadmat
from scipy.signal import lfilter
from scipy.interpolate import griddata
from scipy.constants import mu_0, pi
from scipy.integrate import cumtrapz
from datetime import datetime 
from get_mds_data import get_mds_data
from B_EF import B_EF
from get_TF_current import get_TF_current
from bz_rbfinterp import bz_rbfinterp
from sympy.physics.vector import curl

def process_PCBdata_200ch(PCB, pathname):
    date = PCB.date
    shot = PCB.shot
    tfshot = PCB.tfshot
    n = PCB.n
    i_EF = PCB.i_EF
    trange = PCB.trange
    idx = PCB.idx
    
    filename = os.path.join(pathname['pre_processed_directory_path'], f'{date}{idx:03d}_200ch.mat')
    if not os.path.exists(filename):
        doCalculation = True
        print('No processed data -- start calculation')
    else:
        doCalculation = False
        print('Loading processed data')
    
    if doCalculation:
        xls = pd.ExcelFile('coeff200ch.xlsx')
        sheet_names = xls.sheet_names
        
        sheet_dates = [pd.to_datetime(name) for name in sheet_names]
        today = datetime.today()
        sheet_date = max(date for date in sheet_dates if date <= today)
        C = pd.read_excel('coeff200ch.xlsx', sheet_name=str(sheet_date.date()))
        C_200ch = C.iloc[0:192, :].to_numpy()
        r_shift = 0.00
        ok = C_200ch[:, 13].astype(bool)
        dtacq_num_list = C_200ch[:, 0]
        dtaq_ch = C_200ch[:, 1]
        polarity = C_200ch[:, 12]
        coeff = C_200ch[:, 11]
        zpos = C_200ch[:, 8]
        rpos = C_200ch[:, 9] + r_shift
        ch = C_200ch[:, 6]
        
        if 39 in dtacq_num_list:
            filename1 = os.path.join(os.environ['rawdata_path'], f'mag_probe/dtacq39/shot{shot[0]}_tfshot{tfshot[0]}.mat')
            if not os.path.exists(filename1):
                print('No rawdata file of a039 -- Start generating!')
                rawdata_wTF, rawdata_woTF = get_mds_data(39, shot[0], tfshot[0])  # save_dtacq_dataの代わり
                sio.savemat(filename1, {'rawdata_wTF': rawdata_wTF, 'rawdata_woTF': rawdata_woTF})
            rawdata = loadmat(filename1)
            a039_raw = rawdata.get('rawdata_woTF')
            if a039_raw is None:
                print(f"'rawdata_woTF' not found in {filename1}")
            else:
                print('Data loaded successfully')
        
        raw = np.zeros(1000,len(dtaq_ch))
        for i in range(len(dtaq_ch)):
            if dtacq_num_list(i) == 39:
                raw[:, i] = a039_raw[:, dtaq_ch(i)]
            elif dtacq_num_list(i) == 40:
                print('a040 data cannot be used!')
                return
        b = np.multiply(raw,coeff)
        b = np.multiply(b, polarity)
        
        bz = np.zeros(1000,100)
        bt = bz
        ok_bz = np.zeros(100, dtype=bool)
        ok_bt = ok_bz
        zpos_bz = np.zeros(100,1)
        rpos_bz = zpos_bz
        zpos_bt = zpos_bz
        rpos_bt = zpos_bz
            
        windowSize = 8
        bb = (1/windowSize)*np.ones(1,windowSize)
        aa = 1
        
        for i in range(len(ch)):
            b[:, i] = lfilter(bb, aa, b[:, i])
            b[:, i] = b[:, i] - np.mean(b[1:40, i])
            if ch[i] % 2 == 1:
                bz[:, (ch[i]) // 2] = b[:, i]
                ok_bz[(ch[i]) // 2] = ok[i]
                zpos_bz[(ch[i]) // 2] = zpos[i]
                rpos_bz[(ch[i]) // 2] = rpos[i]
            elif ch[i] % 2 == 0:
                bt[:, ch[i] // 2] = b[:, i]
                ok_bt[ch[i] // 2] = ok[i]
                zpos_bt[ch[i] // 2] = zpos[i]
                rpos_bt[ch[i] // 2] = rpos[i]
        
        zprobepcb = np.array([-0.17, -0.1275, -0.0850, -0.0315, -0.0105, 0.0105, 0.0315, 0.0850, 0.1275, 0.17])
        rprobepcb = np.array([0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.30, 0.33]) + r_shift
        rprobepcb_t = np.array([0.07, 0.10, 0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34]) + r_shift

        zq, rq = np.meshgrid(np.linspace(np.min(zpos_bz), np.max(zpos_bz), n), np.linspace(np.min(rpos_bz), np.max(rpos_bz), n))
        zq_probepcb, rq_probepcb = np.meshgrid(zprobepcb, rprobepcb)
        ok_bt_matrix = np.zeros((len(rprobepcb), len(zprobepcb)), dtype=bool)
        ok_bz_matrix = np.zeros((len(rprobepcb), len(zprobepcb)), dtype=bool)

        for i in range(len(ok_bt)):
            if rpos_bt[i] > r_shift:
                index_r = np.abs(rpos_bt[i] - rprobepcb_t) < 0.001
                index_z = zpos_bt[i] == zprobepcb
                ok_bt_matrix += np.rot90(index_r, -1)[:, None] * index_z[None, :] * ok_bt[i]
            index_r = np.abs(rpos_bz[i] - rprobepcb) < 0.001
            index_z = zpos_bz[i] == zprobepcb
            ok_bz_matrix += np.dot(np.dot(np.rot90(index_r, -1), index_z), ok_bz[i])
        
        grid2D = {
            'zq': zq,
            'rq': rq,
            'zprobepcb': zprobepcb,
            'rprobepcb': rprobepcb,
            'rprobepcb_t': rprobepcb_t,
            'ok_bz_matrix': ok_bz_matrix,
            'ok_bt_matrix': ok_bt_matrix
        }
        grid2D_probe = {
            'zq': zq_probepcb, 
            'rq': rq_probepcb, 
            'rq_t': rprobepcb_t
        }
        
        del zq, rq, zprobepcb, rprobepcb, zq_probepcb, rq_probepcb, rprobepcb_t, ok_bz_matrix, ok_bt_matrix

        r_EF = 0.5
        n_EF = 234.0
        
        if date < 221119:
            z1_EF = 0.68
            z2_EF = -0.68
        else:
            z1_EF = 0.78
            z2_EF = -0.78

        Bz_EF, _ = B_EF(z1_EF, z2_EF, r_EF, i_EF, n_EF, grid2D['rq'], grid2D['zq'], False)
        del EF, r_EF, n_EF, i_EF, z_EF
        
        data2D = {
            'psi': np.zeros((grid2D['rq'].shape[0], grid2D['rq'].shape[1], len(trange))),
            'Bz': np.zeros((grid2D['rq'].shape[0], grid2D['rq'].shape[1], len(trange))),
            'Bt': np.zeros((grid2D['rq'].shape[0], grid2D['rq'].shape[1], len(trange))),
            'Bt_th': np.zeros((grid2D['rq'].shape[0], grid2D['rq'].shape[1], len(trange))),
            'Br': np.zeros((grid2D['rq'].shape[0], grid2D['rq'].shape[1], len(trange))),
            'Jt': np.zeros((grid2D['rq'].shape[0], grid2D['rq'].shape[1], len(trange))),
            'Jz': np.zeros((grid2D['rq'].shape[0], grid2D['rq'].shape[1], len(trange))),
            'Jr': np.zeros((grid2D['rq'].shape[0], grid2D['rq'].shape[1], len(trange))),
            'Et': np.zeros((grid2D['rq'].shape[0], grid2D['rq'].shape[1], len(trange))),
            'Lambda': np.zeros((grid2D['rq'].shape[0], grid2D['rq'].shape[1], len(trange))),
            'trange': trange
        }
        
        if 'fourier' in pathname:
            directory_rogo = os.path.join(pathname['fourier'], 'rogowski/')
            current_folder = os.path.join(directory_rogo, str(date))
            rgwfile = os.path.join(current_folder, f"{date}{PCB.idx:03d}.rgw")
            
            if os.path.isfile(rgwfile):
                I_TF, x, aquisition_rate = get_TF_current(PCB, pathname) 
                m0 = 4 * np.pi * 10**-7
                rgwflag = True
            else:
                print(f'No rgw file at shot {PCB.idx}')
                rgwflag = False
        else:
            print('Path to fourier does not exist')
            rgwflag = False
        
        
        ########################### no angle correction ###############################
    for i in range(trange.shape[1]):
        t = trange[i]    
        # Bz二次元補間(線形fit)
        vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t)
        B_z = -Bz_EF + vq
        B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t)
        
        # PSI計算
        data2D['psi'][:, :, i] = cumtrapz(2 * np.pi * B_z * grid2D['rq'][:, 0], grid2D['rq'][:, 0], axis=0)
        
        # Assuming grid2D['rq'] and grid2D['zq'] are 2D arrays
        data2D['Br'][:, :, i], data2D['Bz'][:, :, i] = np.gradient(data2D['psi'][:, :, i], grid2D['zq'][0, :], grid2D['rq'][:, 0])
        data2D['Br'][:, :, i] = -data2D['Br'][:, :, i] / (2.0 * np.pi * grid2D['rq'])
        data2D['Bz'][:, :, i] = data2D['Bz'][:, :, i] / (2.0 * np.pi * grid2D['rq'])
        data2D['Bt'][:, :, i] = B_t
        data2D['Jt'][:, :, i] = curl(grid2D['zq'][0, :], grid2D['rq'][:, 0], data2D['Bz'][:, :, i], data2D['Br'][:, :, i]) / (4 * np.pi * 1e-7)
        
        if rgwflag:
            timing = (x / aquisition_rate) == t
            data2D['Bt_th'][:, :, i] = m0 * I_TF[timing] * 1e3 * 12 / (2 * np.pi * grid2D['rq'])
        
        if i > 0:
            data2D['Et'][:, :, i] = -1 * (data2D['psi'][:, :, i] - data2D['psi'][:, :, i - 1]) / (2 * np.pi * grid2D['rq'] * 1e-6)
        
        mat_contents = sio.loadmat(filename)
        
        data2D = mat_contents['data2D']
        grid2D = mat_contents['grid2D']
    if doCalculation:
        sio.savemat(filename, {'data2D': data2D, 'grid2D': grid2D, 'shot': shot, 'pathname': pathname, 'filename': filename}) 
    return grid2D, data2D