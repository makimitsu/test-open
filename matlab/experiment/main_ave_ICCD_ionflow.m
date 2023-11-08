%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�V���b�g�ԍ��A�B�e�p�����[�^�Ȃǂ��������O���玩���擾����
%���������ICCD�f�[�^�̕��ς���
%�h�b�v���[�v���[�u�ɂ��C�I�����x�A�t���[�Ƃ��̏u�Ԃ̎��C�ʂ��v���b�g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

%------�yinput�z-------
IDSP.date = 230830;%�yinput�zIDSP������
IDSPshotlist = [11 13 14 16 17 20 22:26 28 29 32 34 37 41:45 47 49 51 54 55 57 59 60];%�yinput�zIDSPshot�ԍ����X�g
IDSP.n_CH = 28;%�yinput�z�h�b�v���[�v���[�u�t�@�C�o�[CH��(28)
IDSP.n_r = 7;%�yinput�z�h�b�v���[�v���[�ur�����f�[�^��(���l)(7)
IDSP.n_z = 1;%�yinput�z�h�b�v���[�v���[�uz�����f�[�^��(���l)(1)
cal_time = 470;

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

shot_group=zeros(1,1);
delay_group=zeros(1,1);
width_group=zeros(1,1);
gain_group=zeros(1,1);
time_group=zeros(1,1);
minZ_group=zeros(1,1);
minR_group=zeros(1,1);

cnt = 0;
for i=1:n_data
    if i > 1
        for j=1:i-1
            %���������shot���~���ŒT��
            k = i-j;
            if (IDSPdelaylist(i) == IDSPdelaylist(k))  && (IDSPwidthlist(i) == IDSPwidthlist(k))...
                    && (IDSPgainlist(i) == IDSPgainlist(k)) && (IDSPtimelist(i) == IDSPtimelist(k))...
                    && (IDSPminZlist(i) == IDSPminZlist(k)) && (IDSPminRlist(i) == IDSPminRlist(k))
                I = find(shot_group == k);
                col_I = mod(I,size(shot_group,1));
                if col_I == 0
                    shot_group = [shot_group;zeros(1,cnt)];
                    I = find(shot_group == k);
                    col_I = mod(I,size(shot_group,1));
                end
                raw_I = idivide(int8(I),int8(size(shot_group,1)))+1;
                shot_group(col_I+1,raw_I) = i;
                break
            end
        end
    end
    if isempty(find(shot_group == i, 1))
        cnt = cnt+1;
        shot_group(1,cnt) = i;
        delay_group(1,cnt)=IDSPdelaylist(i);
        width_group(1,cnt)=IDSPwidthlist(i);
        gain_group(1,cnt)=IDSPgainlist(i);
        time_group(1,cnt)=IDSPtimelist(i);
        minZ_group(1,cnt)=IDSPminZlist(i);
        minR_group(1,cnt)=IDSPminRlist(i);
    end
end

for k=1:size(shot_group,2)
    if time_group(1,k) == cal_time
        i = shot_group(1,k);%�O���[�v�̑�\shot
        PCB.shot=PCBshotlist(i,:);
        PCB.tfshot=PCBtfshotlist(i,:);
        PCB.date = IDSP.date;%�d�˂鎥�C�ʌv����
        PCB.IDX = IDSPshotlist(i);%�d�˂鎥�C��shot�ԍ�
        PCB.i_EF=EFlist(i);
        buf_k = shot_group(:,k);
        buf_k = buf_k(buf_k~=0);
        IDSP.shot = 0;%����p
        IDSP.shotlist=IDSPshotlist(buf_k);
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
            IDSPdata = cal_ave_ICCD_ionflow(IDSP,pathname,show_offset,plot_fit,save_fit);
            %���C�v���[�u�v�Z
            [PCBgrid2D,PCBdata2D] = cal_psi(PCB,pathname);
            %���C�ʃv���b�g
            plot_psi(PCBgrid2D,PCBdata2D,IDSP,FIG,colorplot)
            %�C�I�������v���b�g
            plot_ionflow(IDSPdata,EXP,IDSP,pathname,factor,true,save_fig,'ionflow')
        end
    end
end

