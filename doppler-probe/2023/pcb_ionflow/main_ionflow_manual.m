%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�V���b�g�ԍ��A�B�e�p�����[�^������͂���
%�h�b�v���[�v���[�u�ɂ��C�I�����x�A�t���[���v���b�g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%�������ePC�̃p�X
%�y���R�[�h���g�p����O�Ɂz���ϐ���ݒ肵�Ă������Amatlab���̃R�}���h����setenv('�p�X��','�A�h���X')�Ŏw�肵�Ă��瓮����
setenv("NIFS_path","/Volumes/experiment/results")
setenv("rsGdrive","/Users/rsomeya/Library/CloudStorage/GoogleDrive-rsomeya2016@g.ecc.u-tokyo.ac.jp/�}�C�h���C�u/lab")
pathname.ts3u=getenv('ts3u_path');%old-koala��ts-3u�܂ł̃p�X�imrd�Ȃǁj
pathname.fourier=getenv('fourier_path');%fourier��md0�i�f�[�^�b�N�̃V���b�g�������Ă�j�܂ł�path
pathname.NIFS=getenv('NIFS_path');%results�܂ł�path�i�h�b�v���[�ASXR�j
pathname.save=[getenv('rsGdrive') '/save'];%output�f�[�^�ۑ���
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038��rawdata�̕ۊǏꏊ
pathname.woTFdata=getenv('woTFdata_path');%rawdata�iTFoffset�������j�̕ۊǏꏊ
pathname.fig=[getenv('rsGdrive') '/figure'];%figure�ۑ���
pathname.mat=[getenv('rsGdrive') '/mat'];%figure�ۑ���
pathname.rawdata=[pathname.mat,'/pcb'];%dtacq��rawdata�̕ۊǏꏊ
pathname.flowdata=[pathname.mat,'/ionflow'];%�����f�[�^�̕ۊǏꏊ
pathname.vdistdata=[pathname.mat,'/ionvdist'];%���x���z�f�[�^�̕ۊǏꏊ

%------�yinput�z-------
date = 230315;%�yinput�z������
expval.PF1 = 39;%PF1�d��(kv)
expval.PF2 = 39;%PF2�d��(kv)
expval.TF = 4;%PF2�d��(kv)
expval.EF = 150;%EF�d��
ICCD.shot = 5;%�yinput�z�V���b�g�ԍ�
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
%------�ڍאݒ�yinput�z-------
cal_flow = true;%�yinput�z�������v�Z(true,false)
show_offset = false;%�yinput�z����offset��\��(true,false)
plot_fit = true;%�yinput�z�K�E�X�t�B�b�e�B���O��\��(true,false)
save_fit = false;%�yinput�z�K�E�X�t�B�b�e�B���Opng��ۑ�(true,false)
save_flow = false;%�yinput�z�����f�[�^��ۑ�(true,false)
load_flow = false;%�yinput�z�����f�[�^��ǂݍ���(true,false)
plot_flow = false;%�yinput�z�������v���b�g(true,false)
save_fig = false;%�yinput�z����png��ۑ�(true,false)
factor = 0.1;%�yinput�z�C�I���t���[���T�C�Y(���l:0.1�Ȃ�)

%�v���_�z��𐶐�
mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

%�C�I�����x�A�t���[���v�Z�A�v���b�g
if cal_flow
    %�C�I�����x�A�t���[���v�Z
    [V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_fit,save_flow);
elseif load_flow
    %�ۑ��ς݃C�I�����x�A�t���[��ǂݎ��
    [V_i,absV,T_i] = load_ionflow(date,ICCD,pathname);
end

if plot_flow
    if isempty(V_i)
        return
    else
        plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,false,save_fig,'ionflow')
    end
end