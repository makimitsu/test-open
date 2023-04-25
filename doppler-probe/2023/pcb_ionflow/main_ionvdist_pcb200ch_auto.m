%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�V���b�g�ԍ��A�B�e�p�����[�^�Ȃǂ��������O���玩���擾����
%�h�b�v���[�v���[�u�ɂ��C�I�����x���z�֐��A���x�A�t���[�Ƃ��̏u�Ԃ̎��C�ʂ��v���b�g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

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

%------�yinput�z---------------------------------------------------
date = 230310;%�yinput�z������
begin_cal = 1;%�yinput�z���C��&�t���[�v�Z�n��shot�ԍ�(�������OD��)
end_cal = 0;%�yinput�z���C��&�t���[�v�Z�I���shot�ԍ�(�������OD��)(0�ɂ����begin_cal�ȍ~�̓����̑Sshot�v�Z)
min_r = 12.5;%�yinput�z�h�b�v���[�v���[�u�v���_�ŏ�r���W[mm]
int_r = 2.5;%�yinput�z�h�b�v���[�v���[�u�v���_r�����Ԋu[mm]
min_z = 2.1;%�yinput�z�h�b�v���[�v���[�u�v���_�ŏ�z���W[mm](-2.1,2.1)
int_z = 4.2;%�yinput�z�h�b�v���[�v���[�u�v���_z�����Ԋu[mm](4.2)
ICCD.line = 'Ar';%�yinput�z�h�b�v���[�������C��('Ar')
n_CH = 28;%�yinput�z�h�b�v���[�v���[�u�t�@�C�o�[CH��(28)
n_z = 1;%�yinput�z�h�b�v���[�v���[�uz�����f�[�^��(���l)(1)
%-----------------------�ڍאݒ�yinput�z----------------------------
cal_vdist = true;%�yinput�z���x���z���v�Z(true,false)
plot_spectra = false;%�yinput�z�X�y�N�g�����v���b�g(true,false)
plot_analisis = false;%�yinput�z�t�ϊ���͂��v���b�g(true,false)
plot_vdist = false;%�yinput�z���x���z���v���b�g(true,false)
plot_compare = false;%�yinput�z�č\����r���v���b�g(true,false)
plot_flow = true;%�yinput�z�������v���b�g(true,false)
plot_psi = false;%�yinput�z���C�ʂ��v���b�g(true,false)
overlay_plot = false;%�yinput�z�����Ǝ��C�ʂ��d�˂�(true,false)

save_fit = false;%�yinput�z�K�E�X�t�B�b�e�B���Opng��ۑ�(true,false)
save_fig = true;%�yinput�z���x���z�A�t���[png��ۑ�(true,false)

save_vdist = true;%�yinput�z���x���z�f�[�^��ۑ�(true,false)
load_vdist = false;%�yinput�z���x���z�f�[�^��ǂݍ���(true,false)

plot_type = 'contour';%�yinput�z���x���z�v���b�g���('contour','surf')
Ti_type = 'dispersion';%�yinput�z�C�I�����x�v�Z�@('dispersion')

show_offset = false;%�yinput�z����offset��\��(true,false)
inversion_method = 5;%�yinput�z���x���z�t�ϊ���@(1~6)
factor = 0.1;%�yinput�z�C�I���t���[���T�C�Y(���l:0.1�Ȃ�)
dtacq.num = 39;%�yinput�z���C�v���[�udtacq�ԍ�(39)
mesh_rz = 50;%�yinput�z���C�v���[�urz�����̃��b�V����(50)
trange = 430:590;%�yinput�z���C�v���[�u�v�Z���Ԕ͈�(430:590)
%------------------------------------------------------------------

%�h�b�v���[�v���[�u�v���_�z��𐶐�
mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

%�������O�ǂݎ��
[exp_log,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

%--------���C��&�t���[���v�Z------
start_i = begin_row + begin_cal - 1;
if start_i <= end_row
    if end_cal == 0
        end_i = end_row;%begin_cal�ȍ~�S���v�Z
    elseif end_cal < begin_cal
        error('end_cal must <= begin_cal.')
    elseif begin_row + end_cal - 1 <= end_row
        end_i = begin_row + end_cal - 1;%begin_cal����end_cal�܂Ōv�Z
    else
        error('end_cal must <= %d.', exp_log(end_row,4))
    end
    for i = start_i:end_i
        ICCD.shot = exp_log(i,4);%�V���b�g�ԍ�
        a039shot = exp_log(i,8);%a039�V���b�g�ԍ�
        a039tfshot = exp_log(i,9);%a039TF�V���b�g�ԍ�
        expval.PF1 = exp_log(i,11);%PF1�d��(kv)
        expval.PF2 = exp_log(i,14);%PF2�d��(kv)
        expval.TF = exp_log(i,18);%PF2�d��(kv)
        expval.EF = exp_log(i,23);%EF�d��
        ICCD.trg = exp_log(i,42);%ICCD�g���K����
        ICCD.exp_w = exp_log(i,43);%ICCD�I������
        ICCD.gain = exp_log(i,44);%Andor gain
        time = round(ICCD.trg+ICCD.exp_w/2);%���C�ʃv���b�g����
        if dtacq.num == 39
            dtacq.shot = a039shot;
            dtacq.tfshot = a039tfshot;
        end
        if cal_vdist
            %�C�I�����x���z���v�Z
            [V_i,absV,T_i] = cal_ionvdist(date,expval,ICCD,mpoints,pathname,show_offset,plot_spectra, ...
                inversion_method,plot_analisis,plot_vdist,plot_type,save_fig,plot_compare,save_vdist,Ti_type);
        elseif load_vdist
            %�ۑ��ς݃C�I�����x�A�t���[��ǂݎ��
            [V_i,absV,T_i,F,W,P,Lambda,Vx,Vy,ppoints,Angle] = load_ionvdist(date,ICCD,pathname);
            if plot_vdist
                plot_ionvdist(Vx,Vy,F,date,expval,ICCD,pathname,mpoints,ppoints,plot_type,save_fig)
            end
            if plot_compare
                plot_inversion_compare(F,W,P,Lambda,mpoints,Angle,ICCD)
            end
        else
            V_i = char.empty;
            absV = char.empty;
            T_i = char.empty;
        end
        %���C�ʂ��v���b�g
        if plot_psi
            plot_psi200ch_at_t(time,date,dtacq,pathname,mesh_rz,expval,trange,false);
        end
        %�C�I�����x�A�t���[���v���b�g
        if plot_flow
            if not(isempty(V_i))
                if plot_psi
                    plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,overlay_plot,save_fig,'ionvdist')
                else
                    plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,false,save_fig,'ionvdist')
                end
            end
        end
    end
else
    error('begin_cal must <= %d.', exp_log(end_row,4))
end
