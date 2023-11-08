%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータなどを実験ログから自動取得して
%ドップラープローブによるイオン温度、フローとその瞬間の磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

%------【input】-------
IDSP.date = 230828;%【input】IDSP実験日
IDSPshotlist = [5 6 8:12 15:17 19:23 25:27 31:61 63];%【input】IDSPshot番号リスト
IDSP.n_CH = 28;%【input】ドップラープローブファイバーCH数(28)
IDSP.n_r = 7;%【input】ドップラープローブr方向データ数(数値)(7)
IDSP.n_z = 1;%【input】ドップラープローブz方向データ数(数値)(1)

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

for i=1:n_data
    PCB.shot=PCBshotlist(i,:);
    PCB.tfshot=PCBtfshotlist(i,:);
    PCB.date = IDSP.date;%重ねる磁気面計測日
    PCB.IDX = IDSPshotlist(i);%重ねる磁気面shot番号
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
    FIG.start = IDSP.time;%プロット開始時刻[us]
    if not(isnan(IDSP.delay)) && not(isnan(PCB.shot(1,1))) && not(isnan(PCB.tfshot(1,1))) && not(isnan(PCB.shot(1,2))) && not(isnan(PCB.tfshot(1,2)))
        %IDSP計算
        IDSPdata = cal_ionflow(IDSP,pathname,show_offset,plot_fit,save_fit);
        %磁気プローブ計算
        [PCBgrid2D,PCBdata2D] = cal_psi(PCB,pathname);
        %磁気面プロット
        plot_psi(PCBgrid2D,PCBdata2D,IDSP,FIG,colorplot)
        %イオン流速プロット
        plot_ionflow(IDSPdata,EXP,IDSP,pathname,factor,true,save_fig,'ionflow')
    end
end

