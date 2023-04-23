%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�V���b�g�ԍ��A�B�e�p�����[�^�Ȃǂ��������O���玩���擾����
%�h�b�v���[�v���[�u�ɂ��C�I�����x�A�t���[�Ƃ��̏u�Ԃ̎��C�ʂ��v���b�g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
pathname.flowdata=[getenv('rsGdrive') '/ionflow'];%�����f�[�^�̕ۊǏꏊ

%------�yinput�z-------
date = 230313;%�yinput�z������
begin_cal = 6;%�yinput�z���C��&�t���[�v�Z�n��shot�ԍ�(�������OD��)
end_cal = 6;%�yinput�z���C��&�t���[�v�Z�I���shot�ԍ�(�������OD��)(0�ɂ����begin_cal�ȍ~�̓����̑Sshot�v�Z)
min_r = 12.5;%�yinput�z�h�b�v���[�v���[�u�v���_�ŏ�r���W[mm]
int_r = 2.5;%�yinput�z�h�b�v���[�v���[�u�v���_r�����Ԋu[mm]
min_z = 2.1;%�yinput�z�h�b�v���[�v���[�u�v���_�ŏ�z���W[mm]
int_z = 4.2;%�yinput�z�h�b�v���[�v���[�u�v���_z�����Ԋu[mm]
ICCD.line = 'Ar';%�yinput�z�h�b�v���[�������C��('Ar')
n_CH = 28;%�yinput�z�h�b�v���[�v���[�u�t�@�C�o�[CH��(28)
n_z = 1;%�yinput�z�h�b�v���[�v���[�uz�����f�[�^��(���l)(1)
%------�ڍאݒ�yinput�z-------
factor = 0.05;%�yinput�z�C�I���t���[���T�C�Y(���l:0.05�Ȃ�)
show_offset = false;%�yinput�z����offset��\��(true,false)
plot_fit = false;%�yinput�z�K�E�X�t�B�b�e�B���O��\��(true,false)
cal_flow = true;%�yinput�z�������v�Z(true,false)
save_flow = true;%�yinput�z�����f�[�^��ۑ�(true,false)
load_flow = false;%�yinput�z�����f�[�^��ǂݍ���(true,false)
plot_psi = true;%�yinput�z���C�ʂ��v���b�g(true,false)
plot_flow = true;%�yinput�z�������v���b�g(true,false)
overlay_plot = true;%�yinput�z�����Ǝ��C�ʂ��d�˂�(true,false)
save_fig = false;%�yinput�z����fig��ۑ�(true,false)
dtacq_num = 39;%�yinput�z���C�v���[�udtacq�ԍ�(39)
mesh_rz = 50;%�yinput�z���C�v���[�urz�����̃��b�V����(50)
trange = 430:590;%�yinput�z���C�v���[�u�v�Z���Ԕ͈�(430:590)

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
        i_EF = exp_log(i,23);%EF�d��
        ICCD.trg = exp_log(i,42);%ICCD�g���K����
        ICCD.exp_w = exp_log(i,43);%ICCD�I������
        ICCD.gain = exp_log(i,44);%Andor gain
        time = round(ICCD.trg+ICCD.exp_w/2);%���C�ʃv���b�g����
        if dtacq_num == 39
            dtacq_shot = a039shot;
            dtacq_tfshot = a039tfshot;
        end
        if cal_flow
            %�C�I�����x�A�t���[���v�Z
            [V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_flow);
        elseif load_flow
            %�ۑ��ς݃C�I�����x�A�t���[��ǂݎ��
            [V_i,absV,T_i] = load_ionflow(date,ICCD,pathname);
        end
        %���C�ʂ��v���b�g
        if plot_psi
            plot_psi200ch_at_t(time,date,dtacq_num,dtacq_shot,dtacq_tfshot,pathname,mesh_rz,i_EF,trange,false);
        end
        %�C�I�����x�A�t���[���v���b�g
        if plot_flow
            if isempty(V_i)
                return
            else
                if plot_psi
                    plot_ionflow(V_i,absV,T_i,date,ICCD,factor,mpoints,overlay_plot,save_fig)
                else
                    plot_ionflow(V_i,absV,T_i,date,ICCD,factor,mpoints,false,save_fig)
                end
            end
        end
    end
else
    error('begin_cal must <= %d.', exp_log(end_row,4))
end

