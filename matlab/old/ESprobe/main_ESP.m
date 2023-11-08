%%% cal & plot ExB drift velocity %%%
% clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';

run define_path.m

% %直接入力の場合
date = 230817;%計測日
n_mesh = 50;%メッシュ数
trange = 470:0.1:490;%【input】計算時間範囲(0.1刻み)
IDXlist = [5:7 9 11:12 14];%shot番号

tate = 4;%プロット枚数(縦)
yoko = 2;%プロット枚数(横)
start_t = 470;%プロット開始時刻[us]
dt = 2;%プロット時間間隔[us]
plot_Efield = false;%電場ベクトルをプロット

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
n_data=numel(IDXlist);%計測データ数
rlist=T.ESProbeRPosition_mm_(IDXlist);%静電プローブr座標[mm]

% data2D = cal_ESP(pathname,date,IDXlist,rlist,trange,n_mesh);
plot_ESP(tate,yoko,start_t,dt,plot_Efield,trange,data2D)
% movie_ESP(plot_Efield,trange,data2D)