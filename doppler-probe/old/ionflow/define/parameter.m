Vc = 299792.458;%����(km/s)

if gas == 1
    %�A���S���̎�
    A = 40;%���q��
    lambda0 = 480.602;%�g�p�X�y�N�g��(nm)
    lambda1 = 480.7019;%�Z�������v�X�y�N�g��(nm)
    lambda2 = 479.2619;%�Z�������v�X�y�N�g��(nm)
elseif gas == 2
    %���f�̎�
    A = 1;%���q��
    lambda0 = 486.135;%�g�p�X�y�N�g��(nm)
end

NofCH = 28;%�`�����l����
Angle = 30;%���ˊp(�x)
width = 3;%Y�����������킹��
l = 7;%�ړ����ϒ���
