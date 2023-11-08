%%% cal & plot ExB drift velocity %%%
% clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';

run define_path.m

ESP.date = 230830;%【input】静電プローブ計測日
ESP.shotlist = [11 13 14 16 17 20 22 23];%【input】静電プローブ解析shotlist(同一オペレーション)
ESP.mesh = 50;%【input】静電プローブ補間メッシュ数
ESP.trange = 460:0.1:530;%【input】計算時間範囲(0.1刻み)
ESP.tate = 2;%【input】プロット枚数(縦)
ESP.yoko = 3;%【input】プロット枚数(横)
ESP.start_t = 460;%【input】プロット開始時刻[us]
ESP.dt = 2;%【input】プロット時間間隔[us]
ESP.vector = false;%【input】電場ベクトルをプロット

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,ESP.date);
ESP.rlist=T.ESProbeRPosition_mm_(ESP.shotlist);%静電プローブr座標[mm]

ESPdata2D = cal_ESP(pathname,ESP);
plot_ESP(ESP,ESPdata2D)
% movie_ESP(plot_Efield,trange,data2D)