%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�h�b�v���[�v���[�u�ɂ��C�I�����x���z�֐��A���x�A�t���[�Ƃ��̏u�Ԃ̎��C�ʂ��v���b�g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

%------�yinput�z---------------------------------------------------
IDSP.date = 230830;%�yinput�zIDSP������
IDSPshotlist = [13 14 16 17];% 20 22 23 25 26 28 29 32 37 41:45 47 49 51 54 55 57];%�yinput�zIDSPshot�ԍ����X�g
IDSP.n_CH = 28;%�yinput�z�h�b�v���[�v���[�u�t�@�C�o�[CH��(28)
IDSP.n_r = 7;%�yinput�z�h�b�v���[�v���[�ur�����f�[�^��(���l)(7)
IDSP.n_z = 1;%�yinput�z�h�b�v���[�v���[�uz�����f�[�^��(���l)(1)
%-----------------------�ڍאݒ�yinput�z----------------------------
plot_spectra = false;%�yinput�z�X�y�N�g�����v���b�g(true,false)
plot_analisis = false;%�yinput�z�t�ϊ���͂��v���b�g(true,false)
plot_vdist = true;%�yinput�z���x���z���v���b�g(true,false)
plot_compare = false;%�yinput�z�č\����r���v���b�g(true,false)
save_fig = true;%�yinput�z���x���z�A�t���[png��ۑ�(true,false)
plot_type = 'contour';%�yinput�z���x���z�v���b�g���('contour','surf')
Ti_type = 'dispersion';%�yinput�z�C�I�����x�v�Z�@('dispersion')
show_offset = true;%�yinput�z����offset��\��(true,false)
inversion_method = 5;%�yinput�z���x���z�t�ϊ���@(1~6)
factor = 0.001;%�yinput�z�C�I���t���[���T�C�Y(���l:0.1�Ȃ�)

colorplot = 'none';%�yinput�z�J���[�v���b�g���('psi','Et','Bz','Br','Bt_ext','Bt_plasma','Jt')
PCB.mesh = 40; %�yinput�zpsi��rz�������b�V����
PCB.trange = 400:800;%�yinput�zpsi�v�Z���Ԕ͈�
%------------------------------------------------------------------

%--------�������O�擾---------
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%�X�v���b�h�V�[�g��ID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,IDSP.date);
IDSPshotlist(ismember(IDSPshotlist,T.shot)==0) = [];
n_data=numel(IDSPshotlist);%�v���f�[�^��
PCBshotlist_a039 =T.a039(IDSPshotlist);
PCBshotlist_a040 = T.a040(IDSPshotlist);
PCBshotlist = [PCBshotlist_a039, PCBshotlist_a040];
PCBtfshotlist_a039 =T.a039_TF(IDSPshotlist);
PCBtfshotlist_a040 =T.a040_TF(IDSPshotlist);
PCBtfshotlist = [PCBtfshotlist_a039, PCBtfshotlist_a040];
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

FIG.tate = 1;%�v���b�g����(�c)
FIG.yoko = 1;%�v���b�g����(��)
FIG.dt = 0;%�v���b�g���ԊԊu[us]

%--------���C��&�t���[���v�Z------
for i=1:n_data
    PCB.shot=PCBshotlist(i,:);
    PCB.tfshot=PCBtfshotlist(i,:);
    PCB.date = IDSP.date;%�d�˂鎥�C�ʌv����
    PCB.IDX = IDSPshotlist(i);%�d�˂鎥�C��shot�ԍ�
    PCB.i_EF=EFlist(i);
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
    if not(isnan(IDSP.delay)) && not(isnan(PCB.shot(1,1))) && not(isnan(PCB.tfshot(1,1))) && not(isnan(PCB.shot(1,2))) && not(isnan(PCB.tfshot(1,2)))
        [IDSPdata] = cal_ionvdist(IDSP,pathname,show_offset,plot_spectra,inversion_method,plot_analisis,Ti_type);
        %���x���z���v���b�g
        if plot_vdist
            plot_ionvdist(IDSPdata,IDSP,EXP,pathname,plot_type,save_fig)
        end
        %�č\�����ʂ��瓾����X�y�N�g�����v�Z�A�X�y�N�g���f�[�^�Ɣ�r
        if plot_compare
            plot_inversion_compare(IDSPdata,IDSP)
        end
        %���C�v���[�u�v�Z
        [PCBgrid2D,PCBdata2D] = cal_psi(PCB,pathname);
        %���C�ʃv���b�g
        plot_psi(PCBgrid2D,PCBdata2D,IDSP,FIG,colorplot)
        %�C�I�������v���b�g
        plot_ionflow(IDSPdata,EXP,IDSP,pathname,factor,true,save_fig,'ionvdist')
    end
end
