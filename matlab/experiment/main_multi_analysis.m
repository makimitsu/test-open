%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 磁気プローブ、静電プローブによる
% 磁気面、静電ポテンシャル、ExBドリフトをプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

% %FC合体、X点R=0.2m、ExBアウトフロー小。IDSP->230828,230829(delay=480,484,488us)
% ESP.date = 230828;%【input】静電プローブ計測日
% ESP.shotlist = [5 6 8:12 14:17 19:23 25:27 31:61 63];
% FIG.start = 482;%【input】プロット開始時刻[us]
% PCB.date = 230828;%【input】重ねる磁気面計測日
% PCB.IDX = 15;

%SEP合体、X点R=0.26m、ExBアウトフロー大。IDSP->230830,230831(delay=468, 472, 476us)
ESP.date = 230830;%【input】静電プローブ計測日
ESP.shotlist = [11 13 14 16 17 20 22:26 28 29 32 34 37 41:45 47 49:51 54:57 59 60];%【input】静電プローブ解析shotlist(同一オペレーション)
FIG.start = 470;%【input】プロット開始時刻[us]
PCB.date = 230830;%【input】重ねる磁気面計測日
PCB.IDX = 22;%【input】重ねる磁気面shot番号

FIG.dt = 0;%【input】プロット時間間隔[us]
FIG.tate = 1;%【input】プロット枚数(縦)
FIG.yoko = 1;%【input】プロット枚数(横)

ESP.mesh = 40;%【input】静電プローブ補間メッシュ数
ESP.trange = 460:0.1:500;%【input】計算時間範囲(0.1刻み)
ESP.vector = true;%【input】電場ベクトルをプロット

colorplotlist = ["Bz","Br","Bt_ext","Ez","Er","Et","psi","phi","Jt"];
% ;%【input】カラープロット種類('phi','psi','Ez','Er','Et',...
% 'Bz','Br','Bt_ext','Bt_plasma','absB','absB2','Jt','VExBr','VExBz','|VExB|')
n_clist = size(colorplotlist,2);

PCB.mesh = 40; %【input】psiのrz方向メッシュ数
PCB.trange = 400:800;%【input】psi計算時間範囲

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,ESP.date);
ESP.rlist=T.ESProbeRPosition_mm_(ESP.shotlist);%静電プローブr座標[mm]
shot_a039 =T.a039(PCB.IDX);
shot_a040 = T.a040(PCB.IDX);
PCB.shot = [shot_a039, shot_a040];
tfshot_a039 =T.a039_TF(PCB.IDX);
tfshot_a040 =T.a040_TF(PCB.IDX);
PCB.tfshot = [tfshot_a039, tfshot_a040];
PCB.i_EF=T.EF_A_(PCB.IDX);
IDSPminZlist=T.IDSPZ_cm_(PCB.IDX);
IDSPminRlist=T.IDSPMinR_cm_(PCB.IDX);

IDSP.z = IDSPminZlist*ones(7,1)*1E-2;
IDSP.r = (IDSPminRlist:2.5:IDSPminRlist+6*2.5)*1E-2;

%静電プローブ計算
ESPdata2D = cal_ESP(pathname,ESP);
%磁気プローブ計算
[PCBgrid2D,PCBdata2D] = cal_psi(PCB,pathname);
%ExBドリフト計算
[ExBdata2D,newPCBdata2D] = cal_ExB(pathname,PCBgrid2D,PCBdata2D,ESPdata2D,ESP,PCB,FIG);

%磁気面、ExBドリフト2次元プロット
figure('Position', [0 0 1500 1500],'visible','on')
tiledlayout(3,3)
for i=1:n_clist
    ax = nexttile;
    % subplot(3,3,i)
    plot_ExB(PCBdata2D,ESPdata2D,ExBdata2D,newPCBdata2D,IDSP,FIG,colorplotlist(i),true)
    switch colorplotlist(i)
        case {'phi','Ez','Er','Et','Jt'}
            colormap(ax,redblue(3000));
        case {'psi','Bz','Br','Bt_ext','Bt_plasma','absB','absB2','VExBr','VExBz','|VExB|'}
            colormap(ax,jet)
    end
end
sgtitle([num2str(FIG.start) 'us'])
fontsize(14,"points")
