import os
import pandas as pd
import numpy as np

def directory_generation_Rogowski(date, shot, directory_rogo):
    date_str = str(date)
    if shot < 10:
        shot_str = f'00{shot}'
    elif shot < 100:
        shot_str = f'0{shot}'
    elif shot < 1000:
        shot_str = str(shot)
    else:
        print('More than 999 shots! You need some rest!!!')
        return None, None, None
    data_dir = os.path.join(directory_rogo, date_str, f'{date_str}{shot_str}.txt')
    return date_str, shot_str, data_dir

def get_TF_current(PCB, pathname):
    directory_rogo = os.path.join(pathname['fourier'], 'rogowski/')
    shot = PCB['idx']
    date = PCB['date']
    aquisition_rate = 10
    offset = 0
    
    if date == 230920:
        date = 230929
        shot = 8
    
    t_start = 1 # us
    t_end = 1000 # us
    
    t_start -= offset
    t_end -= offset
    time_step = 1 # us; time step of plot
    
    calibration = 116.6647 # calibration factor for TF coil
    
    date_str, shot_str, path = directory_generation_Rogowski(date, shot, directory_rogo)

    if date_str is None or shot_str is None or path is None:
        return
    
    current_folder = os.path.join(directory_rogo, date_str, '/')
    path = os.path.join(current_folder, f'{date_str}{shot_str}.rgw')
    
    if not os.path.isfile(path):
        print(f'No such file: {path}')
        return
    elif t_start < 0:
        print('offset should be smaller than t_start!')
        return
    elif aquisition_rate * time_step < 1:
        print('time resolution should be < aquisition rate!')
        return
    
    data = pd.read_csv(path, delim_whitespace=True, header=None).to_numpy()
    step = aquisition_rate * time_step
    x = np.arange(t_start * aquisition_rate, t_end * aquisition_rate + 1, step, dtype=int)
    I_TF = data[x, 2] * calibration