%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�V���b�g�ԍ��A�B�e�p�����[�^�Ȃǂ��������O���玩���擾����
%�h�b�v���[�v���[�u�ɂ��C�I�����x�A�t���[�Ƃ��̏u�Ԃ̎��C�ʂ��v���b�g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = main_pcb_ionflow_auto(show_offset,plot_fit,plot_flow,save_flow,save_fig,plot_psi)
%offset��\��/�K�E�X�t�B�b�e�B���O��\��/�������v���b�g/�����f�[�^��ۑ�/����fig��ۑ�/���C�ʂ��v���b�g
%�S��true or false

%%%%%�������ePC�̃p�X
%�y���R�[�h���g�p����O�Ɂz���ϐ���ݒ肵�Ă������Amatlab���̃R�}���h����setenv('�p�X��','�A�h���X')�Ŏw�肵�Ă��瓮����
pathname.ts3u=getenv('ts3u_path');%old-koala��ts-3u�܂ł̃p�X�imrd�Ȃǁj
pathname.fourier=getenv('fourier_path');%fourier��md0�i�f�[�^�b�N�̃V���b�g�������Ă�j�܂ł�path
pathname.NIFS=getenv('NIFS_path');%results�܂ł�path�i�h�b�v���[�ASXR�j
pathname.save=getenv('savedata_path');%output�f�[�^�ۑ���
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038��rawdata�̕ۊǏꏊ
pathname.woTFdata=getenv('woTFdata_path');%rawdata�iTFoffset�������j�̕ۊǏꏊ
% setenv("rsGdrive","/Users/rsomeya/Library/CloudStorage/GoogleDrive-rsomeya2016@g.ecc.u-tokyo.ac.jp/�}�C�h���C�u")
pathname.rawdata=[getenv('rsGdrive') '/pcb'];%dtacq��rawdata�̕ۊǏꏊ

%------�yinput�z-------
date = 230313;%�yinput�z������
begin_cal = 3;%�yinput�z���C��&�t���[�v�Z�n��shot�ԍ�(�������OD��)
end_cal = 3;%�yinput�z���C��&�t���[�v�Z�I���shot�ԍ�(�������OD��)(0�ɂ����begin_cal�ȍ~�̓����̑Sshot�v�Z)
min_r = 12.5;%�yinput�z�h�b�v���[�v���[�u�v���_�ŏ�r���W[mm]
int_r = 2.5;%�yinput�z�h�b�v���[�v���[�u�v���_r�����Ԋu[mm]
min_z = 2.1;%�yinput�z�h�b�v���[�v���[�u�v���_�ŏ�z���W[mm]
int_z = 4.2;%�yinput�z�h�b�v���[�v���[�u�v���_z�����Ԋu[mm]
gas = 'Ar';%�yinput�z�K�X��('Ar')
NofCH = 28;%�yinput�z�t�@�C�o�[CH��(28)
nz = 1;%�yinput�zz�����f�[�^��(���l)(1)
factor = 0.05;%�yinput�z���T�C�Y(���l:0.05�Ȃ�)
dtacq_num = 39;%�yinput�z���C�v���[�udtacq�ԍ�
mesh_rz = 50;%�yinput�z���C�v���[�urz�����̃��b�V����(�قڌŒ�)
trange = 430:590;%�yinput�z���C�v���[�u�v�Z���Ԕ͈�(�قڌŒ�)

%�v���_�z��𐶐�
[r_measured,z_measured] = make_mpoints(NofCH,min_r,int_r,nz,min_z,int_z);

%--------------�������O�ǂݎ��---------------
%���Z�p(�Œ�l�ŗǂ�)
load_s = 5300;%�������O�ǂݎn�ߍs�ԍ�(230309~)
load_f = 10000;%�������O�ǂݏI���s�ԍ�(�Œ�)
%�������O��cd�ɂȂ���΍ŐV�̎������O��web����擾
if exist('exp_log.xlsx','file')
    FileInfo = dir('exp_log.xlsx');
    DateNum = FileInfo.datenum;
    FileDate = datetime(DateNum, 'ConvertFrom', 'datenum', 'Format', 'yyMMdd');
    if str2double(string(FileDate)) > date
        disp(append(string(FileDate), '�X�V�̎������O���g�p���܂��B'))
    else
        run save_log.m
        disp('�ŐV�̎������O(exp_log.xlsx)��ۑ����܂����B')
    end
else
    run save_log.m
    disp('�ŐV�̎������O(exp_log.xlsx)��ۑ����܂����B')
end
%�������O���̎������ɑΉ�����͈͂����
exp_log = readmatrix('exp_log.xlsx','Sheet','log','Range', ['A' num2str(load_s) ':AR' num2str(load_f)]);
[n_row,n_col] = size(exp_log);
begin_row = find(exp_log(:,3) == date);%�������̍ŏ���shot�̍s�ԍ����擾
if isempty(begin_row)
    warning('���������������O���ɑ��݂��܂���B')
    return
end
end_row = begin_row;
while end_row<n_row && isnan(exp_log(end_row+1,3)) && exp_log(end_row+1,4)%���t��NaN&&shot�ԍ����L����=��������shot
    end_row = end_row+1;
end%�������̍Ō��shot�̍s�ԍ����擾

%--------���C��&�t���[���v�Z------
start_i = begin_row + begin_cal - 1;
if start_i <= end_row
    if end_cal == 0
        end_i = end_row;%begin_cal�ȍ~�S���v�Z
    elseif end_cal < begin_cal
        warning('end_cal must <= begin_cal.')
        return
    elseif begin_row + end_cal - 1 <= end_row
        end_i = begin_row + end_cal - 1;%begin_cal����end_cal�܂Ōv�Z
    else
        warning('end_cal must <= %d.', exp_log(end_row,4))
        return
    end
    for i = start_i:end_i
        shot = exp_log(i,4);%�V���b�g�ԍ�
        a039shot = exp_log(i,8);%a039�V���b�g�ԍ�
        a039tfshot = exp_log(i,9);%a039TF�V���b�g�ԍ�
        i_EF = exp_log(i,23);%EF�d��
        trg = exp_log(i,42);%ICCD�g���K����
        exp_w = exp_log(i,43);%ICCD�I������
        gain = exp_log(i,44);%Andor gain
        if dtacq_num == 39
            dtacq_shot = a039shot;
            dtacq_tfshot = a039tfshot;
        end
        if plot_psi
            plot_psi200ch_at_t(round(trg+exp_w/2), date, dtacq_num, dtacq_shot, dtacq_tfshot, pathname,mesh_rz,i_EF,trange,false);
        end
        plot_ionflow(date,shot,trg,exp_w,gain,gas,NofCH,nz,show_offset,plot_fit,plot_flow,save_flow,save_fig,factor,r_measured,z_measured);
    end
else
    warning('begin_cal must <= %d.', exp_log(end_row,4))
    return
end

