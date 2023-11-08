%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 静電プローブによる
% 静電ポテンシャル、電場ベクトルをプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';

run define_path.m

% FC合体、X点R=0.2m、ExBアウトフロー小。IDSP->230828,230829(delay=480,484,488us)
ESP.date = 230828;%【input】静電プローブ計測日
ESP.shotlist = [5 6 8:12 15:17 19:23 25:27 31:61 63];
ESP.start = 482;%【input】プロット開始時刻[us]

% %SEP合体、X点R=0.26m、ExBアウトフロー大。IDSP->230830,230831(delay=468, 472, 476us)
% ESP.date = 230830;%【input】静電プローブ計測日
% ESP.shotlist = [11 13 14 16 17 20 22:26 28 29 32 34 37 41:45 47 49 51 54 55 57 59 60];%【input】静電プローブ解析shotlist(同一オペレーション)
% ESP.start = 470;%【input】プロット開始時刻[us]

ESP.mesh = 21;%【input】静電プローブ補間メッシュ数
ESP.trange = 460:0.1:490;%【input】計算時間範囲(0.1刻み)
ESP.tate = 1;%【input】プロット枚数(縦)
ESP.yoko = 1;%【input】プロット枚数(横)
ESP.dt = 2;%【input】プロット時間間隔[us]
ESP.vector = false;%【input】電場ベクトルをプロット
colorplot = 'Ez';%【input】カラープロット種類('phi','Er','Ez')

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,ESP.date);
ESP.rlist=T.ESProbeRPosition_mm_(ESP.shotlist);%静電プローブr座標[mm]

ESPdata2D = cal_ESP(pathname,ESP);

% plot_ESP(ESP,ESPdata2D,colorplot)
movie_ESP(ESP,ESPdata2D,colorplot)