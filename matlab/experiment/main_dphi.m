%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 磁気プローブ、静電プローブによる
% 磁気面、静電ポテンシャル、ExBドリフトをプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
addpath '/Users/shinjirotakeda/Documents/GitHub/test-open/matlab/common';
addpath '/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/Four-View';

% run define_path.m
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.pre_processed_directory = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）
pathname.ESP=getenv('NIFS_ESP');%smb接続
pathname.fig=[getenv('PROBE_DATA_DIR') '/figure'];%figure保存先
pathname.mat=[getenv('PROBE_DATA_DIR') '/mat'];%mat保存先

% ESP.date = 230826;%計測日
% ESP.shotlist = [6 8 10 11 12 16 17 19 20 22 24 25 26 29 32 34 35 37 38 43 45 47 49 51 53 55];%shot番号
% PCB.date = 230826;%【input】重ねる磁気面計測日
% PCB.IDX = 22;

% ESP.date = 230826;%計測日
% ESP.shotlist = [7 9 13 15 18 21 23 27 28 30 31 33 36 39 40 41 42 44 46 48 50 52 54 56:60];%shot番号
% PCB.date = 230826;%【input】重ねる磁気面計測日
% PCB.IDX = 21;

% % FC合体、X点R=0.2m、ExBアウトフロー小。IDSP->230828,230829(delay=480,484,488us)
% ESP.date = 230828;%【input】静電プローブ計測日
% ESP.shotlist = [5 6 8:12 15:17 19:23 25:27 31:61 63];
% FIG.start = 480;%【input】プロット開始時刻[us]
% FIG.dt = 1;%【input】プロット時間間隔[us]
% PCB.date = 230828;%【input】重ねる磁気面計測日
% PCB.IDX = 15;

% %SEP合体、X点R=0.26m、ExBアウトフロー大。IDSP->230830,230831(delay=468, 472, 476us)
ESP.date = 230830;%【input】静電プローブ計測日
ESP.shotlist = [11 13 14 16 17 20 22:26 28 29 32 34 37 41:45 47 49 51 54 55 57 59 60];%【input】静電プローブ解析shotlist(同一オペレーション)
% FIG.start = 470;%【input】プロット開始時刻[us]
% FIG.dt = 1;%【input】プロット時間間隔[us]
% 
% FIG.tate = 1;%【input】プロット枚数(縦)
% FIG.yoko = 1;%【input】プロット枚数(横)

ESP.mesh = 21;%【input】静電プローブ補間メッシュ数(21)
% ESP.mesh = 40;%【input】静電プローブ補間メッシュ数
ESP.trange = 440:0.1:500;%【input】計算時間範囲(0.1刻み)
ESP.vector = false;%【input】電場ベクトルをプロット
colorplot = 'phi';%【input】カラープロット種類('phi','Er','Ez')
ESP.tate = 1;%【input】プロット枚数(縦)
ESP.yoko = 1;%【input】プロット枚数(横)
ESP.dt = 2;%【input】プロット時間間隔[us]
ESP.start = 474;%【input】プロット開始時刻[us]

% 
% colorplot = 'none';%【input】カラープロット種類('phi','psi','Ez','Er','Et',...
% % 'Bz','Br','Bt_ext','Bt_plasma','absB','absB2','Jt','VExBr','VExBz','|VExB|')
% profileplot = 'VExBz';%【input】一次元プロット種類('VExBr','VExBz','|VExB|')

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,ESP.date);
ESP.rlist=T.ESProbeRPosition_mm_(ESP.shotlist);%静電プローブr座標[mm]
% shot_a039 =T.a039(PCB.IDX);
% shot_a040 = T.a040(PCB.IDX);
% PCB.shot = [shot_a039, shot_a040];
% tfshot_a039 =T.a039_TF(PCB.IDX);
% tfshot_a040 =T.a040_TF(PCB.IDX);
% PCB.tfshot = [tfshot_a039, tfshot_a040];
% PCB.i_EF=T.EF_A_(PCB.IDX);

%静電プローブ計算
ESPdata2D = cal_ESP(pathname,ESP);
d_phi = max(squeeze(max(ESPdata2D.phi,[],2)),[],2) - min(squeeze(min(ESPdata2D.phi,[],2)),[],2);
f = figure;
plot(ESPdata2D.trange,d_phi,'LineWidth',3);
xlabel('time [us]');
ylabel('\Delta\Phi [V]');
ax = gca;
ax.FontSize = 18;