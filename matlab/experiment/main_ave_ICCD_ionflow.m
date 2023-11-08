%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータなどを実験ログから自動取得して
%同一条件のICCDデータの平均から
%ドップラープローブによるイオン温度、フローとその瞬間の磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

%------【input】-------
IDSP.date = 230830;%【input】IDSP実験日
IDSPshotlist = [11 13 14 16 17 20 22:26 28 29 32 34 37 41:45 47 49 51 54 55 57 59 60];%【input】IDSPshot番号リスト
IDSP.n_CH = 28;%【input】ドップラープローブファイバーCH数(28)
IDSP.n_r = 7;%【input】ドップラープローブr方向データ数(数値)(7)
IDSP.n_z = 1;%【input】ドップラープローブz方向データ数(数値)(1)
cal_time = 470;

%------詳細設定【input】------
plot_fit = true;%【input】ガウスフィッティングを表示(true,false)
save_fit = true;%【input】ガウスフィッティングpngを保存(true,false)
save_fig = true;%【input】流速pngを保存(true,false)
show_offset = false;%【input】分光offsetを表示(true,false)
factor = 0.001;%【input】イオンフロー矢印サイズ(数値:0.1など)
colorplot = 'none';%【input】カラープロット種類('psi','Et','Bz','Br','Bt_ext','Bt_plasma','Jt')
PCB.mesh = 40; %【input】psiのrz方向メッシュ数
PCB.trange = 400:800;%【input】psi計算時間範囲

%--------実験ログ取得---------
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,IDSP.date);
IDSPshotlist(ismember(IDSPshotlist,T.shot)==0) = [];
n_data=numel(IDSPshotlist);%計測データ数
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

FIG.tate = 1;%プロット枚数(縦)
FIG.yoko = 1;%プロット枚数(横)
FIG.dt = 0;%プロット時間間隔[us]

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
            %同一条件のshotを降順で探す
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
        i = shot_group(1,k);%グループの代表shot
        PCB.shot=PCBshotlist(i,:);
        PCB.tfshot=PCBtfshotlist(i,:);
        PCB.date = IDSP.date;%重ねる磁気面計測日
        PCB.IDX = IDSPshotlist(i);%重ねる磁気面shot番号
        PCB.i_EF=EFlist(i);
        buf_k = shot_group(:,k);
        buf_k = buf_k(buf_k~=0);
        IDSP.shot = 0;%判定用
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
        FIG.start = IDSP.time;%プロット開始時刻[us]
        if not(isnan(IDSP.delay)) && not(isnan(PCB.shot(1,1))) && not(isnan(PCB.tfshot(1,1))) && not(isnan(PCB.shot(1,2))) && not(isnan(PCB.tfshot(1,2)))
            %IDSP計算
            IDSPdata = cal_ave_ICCD_ionflow(IDSP,pathname,show_offset,plot_fit,save_fit);
            %磁気プローブ計算
            [PCBgrid2D,PCBdata2D] = cal_psi(PCB,pathname);
            %磁気面プロット
            plot_psi(PCBgrid2D,PCBdata2D,IDSP,FIG,colorplot)
            %イオン流速プロット
            plot_ionflow(IDSPdata,EXP,IDSP,pathname,factor,true,save_fig,'ionflow')
        end
    end
end

