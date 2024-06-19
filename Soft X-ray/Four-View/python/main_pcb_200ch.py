from spreadsheet import getTS6log, searchlog
from check_signal import check_signal
from plot_psi200ch import plot_psi200ch
import numpy as np
import os

#環境変数設定
OnoLab = '/Users/shohgookazaki/Documents/UTokyo/OnoLab'
Documents = '/Users/shohgookazaki/Documents'
GoogleDrive = '/Users/shohgookazaki/Library/CloudStorage/GoogleDrive-shohgo-okazaki@g.ecc.u-tokyo.ac.jp/My Drive/OnoLab/data'

os.environ['ts3u'] = OnoLab +'/koala/mnt/koala-experiment/results/ts-3u'
os.environ['fourier_path'] = OnoLab + '/koala/mnt/fourier'
os.environ['NIFS_path'] = OnoLab + '/koala/home/pub/mnt/data/X-ray'
# os.onviron['savedata_path']
# os.environ['woTFdata_path']
os.environ['rawdata_path'] = GoogleDrive + '/probedata/rawdata'
os.environ['pre_processed_directory_path'] = GoogleDrive + '/probedata/processed'

pathname = {
    'ts3u': os.environ['ts3u'],
    'fourier': os.environ['fourier_path'],
    'NIFS': os.environ['NIFS_path'],
    'save': os.environ['savedata_path'],
    'woTFdata': os.environ['woTFdata_path'],
    'rawdata': os.environ['rawdata_path'],
    'pre_processed_directory_path': os.environ['pre_processed_directory_path']
}

class PCBclass:
  def __init__(self, trange, n, start, date, idx, shot, tfshot, dataType, i_EF, TF):
    self.trange = trange
    self.start = start
    self.n = n
    self.date = date
    self.idx = idx
    self.shot = shot
    self.tfshot = tfshot
    self.dataType = dataType
    self.i_EF = i_EF
    self.TF = TF

date = 240501
IDXlist = [12,13]
doCheck = 0
dataType = 1
'''
date = input("date:")
IDXlist = input("shot number:")
doCheck = input("doCheck(0 or 1):")
dataType = input("dataType('1:Bz', '2: psi', '3: Bt', '4: Jt', '5: Et'):")
'''

T = getTS6log()
T = searchlog(T, node = 'date', pat = date)

n_data = len(IDXlist)
shotlist = T[['a039', 'a040']].iloc[IDXlist]
tfshotlist = T[['a039_TF', 'a040_TF']].iloc[IDXlist]
EFlist = T[['EF[A]']].iloc[IDXlist]
TFlist = T[['TF[kV]']].iloc[IDXlist]
dtacqlist =39 * np.ones(n_data)

value = shotlist.iloc[0]
print(shotlist)
print()
trange = 400#:800
n = 40
start =  40

for i in range(n_data):
    print(i)
    PCB = PCBclass(trange = trange, 
                   n = n, start = start, 
                   date = date, idx = IDXlist[i], 
                   shot = shotlist['a039'].loc[shotlist.index[i]], 
                   tfshot = shotlist['a039'].loc[shotlist.index[i]], 
                   dataType = dataType, 
                   i_EF = EFlist, 
                   TF = TFlist)
    if PCB.shot == PCB.tfshot:
        PCB.tfshot = [0,0]
    if doCheck:
        check_signal(PCB)
    else:
        plot_psi200ch(PCB, pathname)

