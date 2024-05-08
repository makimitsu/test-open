%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�V���b�g�ԍ��A�B�e�p�����[�^�Ȃǂ��������O���玩���擾����
%�h�b�v���[�v���[�u�ɂ��C�I�����x�A�t���[�Ƃ��̏u�Ԃ̎��C�ʂ��v���b�g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

%------�yinput�z-------
IDSP.date = 230828;%�yinput�zIDSP������
IDSPshotlist = [5 6 8:12 15:17 19:23 25:27 31:61 63];%�yinput�zIDSPshot�ԍ����X�g
IDSP.n_CH = 28;%�yinput�z�h�b�v���[�v���[�u�t�@�C�o�[CH��(28)
IDSP.n_r = 7;%�yinput�z�h�b�v���[�v���[�ur�����f�[�^��(���l)(7)
IDSP.n_z = 1;%�yinput�z�h�b�v���[�v���[�uz�����f�[�^��(���l)(1)

%------�ڍאݒ�yinput�z------
plot_fit = true;%�yinput�z�K�E�X�t�B�b�e�B���O��\��(true,false)
save_fit = true;%�yinput�z�K�E�X�t�B�b�e�B���Opng��ۑ�(true,false)
save_fig = true;%�yinput�z����png��ۑ�(true,false)
show_offset = false;%�yinput�z����offset��\��(true,false)
factor = 0.001;%�yinput�z�C�I���t���[���T�C�Y(���l:0.1�Ȃ�)
colorplot = 'none';%�yinput�z�J���[�v���b�g���('psi','Et','Bz','Br','Bt_ext','Bt_plasma','Jt')
PCB.mesh = 40; %�yinput�zpsi��rz�������b�V����
PCB.trange = 400:800;%�yinput�zpsi�v�Z���Ԕ͈�

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
        %IDSP�v�Z
        IDSPdata = cal_ionflow(IDSP,pathname,show_offset,plot_fit,save_fit);
        %���C�v���[�u�v�Z
        [PCBgrid2D,PCBdata2D] = cal_psi(PCB,pathname);
        %���C�ʃv���b�g
        plot_psi(PCBgrid2D,PCBdata2D,IDSP,FIG,colorplot)
        %�C�I�������v���b�g
        plot_ionflow(IDSPdata,EXP,IDSP,pathname,factor,true,save_fig,'ionflow')
    end
end

