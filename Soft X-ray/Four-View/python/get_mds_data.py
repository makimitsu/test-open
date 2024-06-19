import numpy as np
import MDSplus as mds

def get_mds_data(dtacq_num, shot, tfshot):

    rawdata_wTF = 0
    rawdata_TF = 0
    
    if dtacq_num  == 38:
        ch_num = 128
    else:
        ch_num = 291
    
    post = 1000
    dtacq = 'a{:03d}'.format(dtacq_num)
    rawdata_wTF = np.zeros((post, ch_num))
    rawdata_TF=rawdata_wTF
    
    mds.Connection('192.168.1.140')
    mds.Tree(dtacq, shot).open()
    
    for i in range(1,ch_num+1):
        # 各チャンネルにおいて「.AI:CHXXX」というノードを指定するためのノード名を作る
        chname = f'.AI:CH{i:03d}'
        # num2strで数値データをstrデータに変換。この時'%03i'で左側を(0)で埋めた(3)桁の整数(i)という形を指定できる。
        data = mds.Tree.getNode(chname).data()
        # データがとれていないときエラーメッセージが多分237文字で帰ってくるので、1000以下の要素はデータなしとしてリターンする
        if len(data) < 1000:
            break
        # Offset adjustment
        data = data - data[0]
        rawdata_wTF.append(data)
    
    rawdata_woTF = rawdata_wTF
    
    if tfshot > 0:
        mds.Tree(dtacq, tfshot).open()
        for i in range(1, ch_num+1):
            chname = f'.AI:CH{i:03d}'
            data = mds.Tree.getNode(chname).data()
            if len(data) < 1000:
                data = [0] * 1000
            data = data - data[0]
            rawdata_TF.append(data)
        rawdata_woTF = [w - t for w, t in zip(rawdata_wTF, rawdata_TF)]

    return  rawdata_wTF, rawdata_woTF