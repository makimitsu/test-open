import os
from process_PCBdata_200ch import process_PCBdata_200ch
from get_axis_x_multi import get_axis_x_multi

def plot_psi200ch(PCB, pathname):
    shot = PCB.shot
    trange = PCB.trange
    start = PCB.start
    
    grid2D, data2D = process_PCBdata_200ch(PCB, pathname)
    if not isinstance(grid2D, dict):  # もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
        return
    magAxisList, xPointList = get_axis_x_multi(grid2D, data2D)