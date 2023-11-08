%物理定数
Vc = 299792.458;%光速(km/s)
mp = 1.67e-27;%陽子質量(kg)
kB = 1.60e-19;%ボルツマン定数(J/eV)

%アルゴンの時
A = 40;%原子量
lambda0 = 480.602;%使用スペクトル(nm)
lambda1 = 480.7019;%校正ランプスペクトル(nm)
lambda2 = 479.2619;%校正ランプスペクトル(nm)
% %水素の時
% A = 1;%原子量
% lambda0 = 486.135;%使用スペクトル(nm)

%装置変数
numMea = 1;%計測点数
numSight = 4;%計測視線数
numCh = numSight * numMea;%チャンネル数 = 視線数 * 計測点数
Angle = [0 30 150];%視線角度[度](0~180)
Theta = Angle*pi/180;%視線角度[rad]に変換
numTheta = numel(Theta);%視線角度数

%解析変数
numLambda = 61;%切り取り波長数(奇数にする)
width = 3;%Y方向足し合わせ幅(前後)
movlen = 3;%移動平均幅
