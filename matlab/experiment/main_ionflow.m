%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�h�b�v���[�v���[�u�ɂ��C�I�����x�A�t���[���v���b�g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

%------�yinput�z-------
% FC���́AX�_R=0.2m�AExB�A�E�g�t���[���BIDSP->230828,230829(delay=480,484,488us)
IDSP.date = 230828;%�yinput�zIDSP������
IDSPshotlist = [5 6 8:12 15:17 19:23 25:27 31:61 63];
cal_time = 482;%�yinput�z�v���b�g����[us]482,486

% % SEP���́AX�_R=0.26m�AExB�A�E�g�t���[��BIDSP->230830,230831(delay=468, 472, 476us)
% IDSP.date = 230830;%�yinput�zIDSP������
% IDSPshotlist = [11 13 14 16 17 20 22:26 28 29 32 34 37 41:45 47 49 51 54 55 57 59 60];%[13 14 16 17 20 22 23 25 26 32 37 41:45 47 49 51 54 55 57];%�yinput�zIDSPshot�ԍ����X�g
% cal_time = 470;%�yinput�z�v���b�g����[us]470,474


IDSP.n_CH = 28;%�yinput�z�h�b�v���[�v���[�u�t�@�C�o�[CH��(28)
IDSP.n_r = 7;%�yinput�z�h�b�v���[�v���[�ur�����f�[�^��(���l)(7)
IDSP.n_z = 1;%�yinput�z�h�b�v���[�v���[�uz�����f�[�^��(���l)(1)

%------�ڍאݒ�yinput�z------
plot_fit = true;%�yinput�z�K�E�X�t�B�b�e�B���O��\��(true,false)
save_fit = true;%�yinput�z�K�E�X�t�B�b�e�B���Opng��ۑ�(true,false)
save_fig = true;%�yinput�z����png��ۑ�(true,false)
show_offset = false;%�yinput�z����offset��\��(true,false)
factor = 0.001;%�yinput�z�C�I���t���[���T�C�Y(���l:0.001�Ȃ�)

%--------�������O�擾---------
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%�X�v���b�h�V�[�g��ID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,IDSP.date);
IDSPshotlist(ismember(IDSPshotlist,T.shot)==0) = [];
n_data=numel(IDSPshotlist);%�v���f�[�^��
PF1list=T.CB1_kV_(IDSPshotlist);
PF2list=T.CB2_kV_(IDSPshotlist);
TFlist=T.TF_kV_(IDSPshotlist);
EFlist=T.EF_A_(IDSPshotlist);
Gaslist=T.gas(IDSPshotlist);
IDSPminZlist=T.IDSPZ_cm_(IDSPshotlist);
IDSPminRlist=T.IDSPMinR_cm_(IDSPshotlist);
IDSPdelaylist=T.IDSPDelay_us_(IDSPshotlist);
IDSPwidthlist=T.IDSPWidth_us_(IDSPshotlist);
IDSPgainlist=T.IDSPGain(IDSPshotlist);
IDSPtimelist=round(IDSPdelaylist+IDSPwidthlist/2);

for i=1:n_data
    IDSP.shot=IDSPshotlist(i);
    IDSP.line=Gaslist(i);
    IDSP.delay=IDSPdelaylist(i);
    IDSP.width=IDSPwidthlist(i);
    IDSP.gain=IDSPgainlist(i);
    IDSP.time=IDSPtimelist(i);
    IDSP.z = IDSPminZlist(i)*ones(7,1)*1E-2;
    IDSP.r = ((IDSPminRlist(i):2.5:IDSPminRlist(i)+(IDSP.n_r-1)*2.5)*1E-2)';
    EXP.PF1=PF1list(i);
    EXP.PF2=PF2list(i);
    EXP.TF=TFlist(i);
    EXP.EF=EFlist(i);
    FIG.start = IDSP.time;%�v���b�g�J�n����[us]
    if not(isnan(IDSP.delay))
        if IDSP.time == cal_time
            %IDSP�v�Z
            IDSPdata = cal_ionflow(IDSP,pathname,show_offset,plot_fit,save_fit);
            %�C�I�������v���b�g
            plot_ionflow(IDSPdata,EXP,IDSP,pathname,factor,false,save_fig,'ionflow')
        end
    end
end

