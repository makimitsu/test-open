%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�V���b�g�ԍ��A�B�e�p�����[�^������͂���
%�h�b�v���[�v���[�u�ɂ��C�I�����x�A�t���[���v���b�g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = main_ionflow_manual(show_offset,plot_fit,save_flow,plot_flow,save_fig)
%����offset��\��/�K�E�X�t�B�b�e�B���O��\��/�����f�[�^��ۑ�/�������v���b�g/����fig��ۑ�
%�S��true or false
%���s��)main_ionflow_manual(false,false,false,true,false)

%%%%%�������ePC�̃p�X
%�y���R�[�h���g�p����O�Ɂz���ϐ���ݒ肵�Ă������Amatlab���̃R�}���h����setenv('�p�X��','�A�h���X')�Ŏw�肵�Ă��瓮����
pathname.ts3u=getenv('ts3u_path');%old-koala��ts-3u�܂ł̃p�X�imrd�Ȃǁj
pathname.fourier=getenv('fourier_path');%fourier��md0�i�f�[�^�b�N�̃V���b�g�������Ă�j�܂ł�path
% setenv("NIFS_path","/Volumes/experiment/results")
pathname.NIFS=getenv('NIFS_path');%results�܂ł�path�i�h�b�v���[�ASXR�j
pathname.save=getenv('savedata_path');%output�f�[�^�ۑ���
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038��rawdata�̕ۊǏꏊ
pathname.woTFdata=getenv('woTFdata_path');%rawdata�iTFoffset�������j�̕ۊǏꏊ
% setenv("rsGdrive","/Users/rsomeya/Library/CloudStorage/GoogleDrive-rsomeya2016@g.ecc.u-tokyo.ac.jp/�}�C�h���C�u")
pathname.rawdata=[getenv('rsGdrive') '/pcb'];%dtacq��rawdata�̕ۊǏꏊ

%------�yinput�z-------
date = 230313;%�yinput�z������
ICCD.shot = 3;%�yinput�z�V���b�g�ԍ�
ICCD.trg = 474;%�yinput�zICCD�g���K����
ICCD.exp_w = 2;%�yinput�zICCD�I������
ICCD.gain = 4095;%�yinput�zICCD gain
ICCD.line = 'Ar';%�yinput�z�h�b�v���[�������C��('Ar')
min_r = 12.5;%�yinput�z�h�b�v���[�v���[�u�v���_�ŏ�r���W[cm]
int_r = 2.5;%�yinput�z�h�b�v���[�v���[�u�v���_r�����Ԋu[cm]
min_z = 2.1;%�yinput�z�h�b�v���[�v���[�u�v���_�ŏ�z���W[cm]
int_z = 4.2;%�yinput�z�h�b�v���[�v���[�u�v���_z�����Ԋu[cm]
n_CH = 28;%�yinput�z�h�b�v���[�v���[�u�t�@�C�o�[CH��(28)
n_z = 1;%�yinput�z�h�b�v���[�v���[�uz�����f�[�^��(���l)(1)
factor = 0.05;%�yinput�z�C�I���t���[���T�C�Y(���l:0.05�Ȃ�)

%�v���_�z��𐶐�
mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

%�C�I�����x�A�t���[���v�Z�A�v���b�g
[V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_flow);
if plot_flow
    plot_ionflow(V_i,absV,T_i,date,ICCD,factor,mpoints,false,save_fig)
end
end