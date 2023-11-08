function [] = auto_main(gas,NofCH,nz,plot_fitting,cal_flow,save_flow,plot_flow,save_fig,factor)
%�K�X��(Ar:1,H:2)/�t�@�C�o�[CH��/z�����f�[�^��(���l)/�t�B�b�e�B���O��\��(TorF)/�������v�Z(TorF)/������ۑ�(TorF)/������\��(TorF)/fig��ۑ�(TorF)/���T�C�Y(���l:0.05�Ȃ�)

%���͕ϐ�
date = 230314;%������
begin_flow = 2;%�t���[�v�Z�n��shot�ԍ�
end_flow = 2;%�t���[�v�Z�I���shot�ԍ�(0�ɂ����begin_flow�ȍ~�Sshot�v�Z)
min_r = 12.5;%�v���[�u�v���_�ŏ�r���W[mm]
int_r = 2.5;%�v���[�u�v���_r�����Ԋu[mm]
min_z = -2.1;%�v���[�u�v���_�ŏ�z���W[mm]
int_z = 4.2;%�v���[�u�v���_z�����Ԋu[mm]

%�v���_�z��𐶐�
[r_measured,z_measured] = make_mpoints(NofCH,min_r,int_r,nz,min_z,int_z);

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
    disp('���������������O���ɑ��݂��܂���B')
    return
end
end_row = begin_row;
while end_row<n_row && isnan(exp_log(end_row+1,3)) && exp_log(end_row+1,4)%���t��NaN&&shot�ԍ����L����=��������shot
    end_row = end_row+1;
end%�������̍Ō��shot�̍s�ԍ����擾

%�������O����ICCD�t�@�C��������肵�A�t���[���v�Z
if end_flow == 0
    %begin_flow�ȍ~�S���v�Z
    for i = begin_row + begin_flow - 1:end_row
        shot = exp_log(i,4);%�V���b�g�ԍ�
        a039shot = exp_log(i,8);%a039�V���b�g�ԍ�
        a039tfshot = exp_log(i,9);%a039TF�V���b�g�ԍ�
        trg = exp_log(i,42);%ICCD�g���K����
        exp_w = exp_log(i,43);%ICCD�I������
        gain = exp_log(i,44);%Andor gain
        manual_main(true,date,shot,trg,exp_w,gain,gas,NofCH,nz,plot_fitting,cal_flow,save_flow,plot_flow,save_fig,factor,r_measured,z_measured)
    end
else
    %begin_flow����end_flow�܂Ōv�Z
    for i = begin_row + begin_flow - 1:begin_row + end_flow - 1
        shot = exp_log(i,4);%�V���b�g�ԍ�
        trg = exp_log(i,42);%ICCD�g���K����
        exp_w = exp_log(i,43);%ICCD�I������
        gain = exp_log(i,44);%Andor gain
        manual_main(true,date,shot,trg,exp_w,gain,gas,NofCH,nz,plot_fitting,cal_flow,save_flow,plot_flow,save_fig,factor,r_measured,z_measured)
    end
end
