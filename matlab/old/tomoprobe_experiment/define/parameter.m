%�����萔
Vc = 299792.458;%����(km/s)
mp = 1.67e-27;%�z�q����(kg)
kB = 1.60e-19;%�{���c�}���萔(J/eV)

%�A���S���̎�
A = 40;%���q��
lambda0 = 480.602;%�g�p�X�y�N�g��(nm)
lambda1 = 480.7019;%�Z�������v�X�y�N�g��(nm)
lambda2 = 479.2619;%�Z�������v�X�y�N�g��(nm)
% %���f�̎�
% A = 1;%���q��
% lambda0 = 486.135;%�g�p�X�y�N�g��(nm)

%���u�ϐ�
numMea = 1;%�v���_��
numSight = 4;%�v��������
numCh = numSight * numMea;%�`�����l���� = ������ * �v���_��
Angle = [0 30 150];%�����p�x[�x](0~180)
Theta = Angle*pi/180;%�����p�x[rad]�ɕϊ�
numTheta = numel(Theta);%�����p�x��

%��͕ϐ�
numLambda = 61;%�؂���g����(��ɂ���)
width = 3;%Y�����������킹��(�O��)
movlen = 3;%�ړ����ϕ�
