%�v���_�z��𐶐�
function [r_measured,z_measured] = make_mpoints(NofCH,min_r,int_r,nz,min_z,int_z)
nr = NofCH/4;
r_measured = zeros(nr,nz);%�x�N�g���v���b�gr���W1���data1�A2���data2
z_measured = zeros(nr,nz);%�x�N�g���v���b�gz���W1���data1�A2���data2
for i = 1:nz
    r_measured(:,i) = linspace(min_r,min_r+(nr-1)*int_r,nr);%datai��r���W
    z_measured(:,i) = min_z+int_z*(i-1);%datai��z���W
end
