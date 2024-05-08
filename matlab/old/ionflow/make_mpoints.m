%計測点配列を生成
function [r_measured,z_measured] = make_mpoints(NofCH,min_r,int_r,nz,min_z,int_z)
nr = NofCH/4;
r_measured = zeros(nr,nz);%ベクトルプロットr座標1列目data1、2列目data2
z_measured = zeros(nr,nz);%ベクトルプロットz座標1列目data1、2列目data2
for i = 1:nz
    r_measured(:,i) = linspace(min_r,min_r+(nr-1)*int_r,nr);%dataiのr座標
    z_measured(:,i) = min_z+int_z*(i-1);%dataiのz座標
end
