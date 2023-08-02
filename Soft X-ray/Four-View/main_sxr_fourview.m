%%%%%%%%%%%%%%%%%%%%%%%%
% Top-level file for calculating and plotting SXR emission for four-view
% experimental setup
%土井プローブ(Bz100ch)の磁気面
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%
clear
%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.NIFS=getenv('NIFS_path');%192.168.1.111
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所;
pathname.pre_processed_directory = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）

%【※コードを使用する前に】Change the addpath() command below to access the
% test-open folder on your computer.
addpath(genpath('/Users/shinjirotakeda/Documents/GitHub/test-open'));
% addpath(genpath('/Users/saadayub/Desktop/test-open'));                     

%%%%実験オペレーションの取得
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
pat = 230721;
date = pat;
T=searchlog(T,node,pat);
% IDXlist = [33,35:40];
IDXlist = 6;
n_data=numel(IDXlist);%計測データ数
shotlist=T.a039(IDXlist);
tfshotlist=T.a039_TF(IDXlist);
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);
startlist = T.SXRStart(IDXlist);
intervallist = T.SXRInterval(IDXlist);

% % %直接入力の場合【注意】全て同じサイズの行列になるように記入
% dtacqlist=39;
% shotlist=1118;%【input】実験ログのa039の番号
% tfshotlist=1106;%【input】実験ログのa039_TFの番号
% date = 230315;%【input】計測日
% n_data=numel(shotlist);%計測データ数
% EFlist = 150;%【input】EF電流

trange=440:500;%【input】計算時間範囲
n=50; %【input】rz方向のメッシュ数

t = 475;
show_xpoint = false;
show_localmax = false;
% start = 450;
% interval = 5;
save = true;
filter = false;
NL = false;

for i=1:n_data
    dtacq_num=dtacqlist(i);
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    i_EF=EFlist(i);
    TF=TFlist(i);
    start = startlist(i);
    interval = intervallist(i);
    % TF=4;
%     plot_psi200ch(date, dtacq_num, shot, tfshot, pathname,n,i_EF,trange,TF); 
    % [grid2D,data2D] = process_PCBdata(date, dtacq_num, shot, tfshot, pathname, n,i_EF,trange);
    [grid2D,data2D] = process_PCBdata_280ch(date, shot, tfshot, pathname, n,i_EF,trange);
    % grid2D = NaN;
    % data2D = NaN;
    shot_SXR = IDXlist(i);
    % shot_SXR = 14;
    % SXRfilename = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/SXR_Images/',num2str(date),'/shot',num2str(shot_SXR,'%03i'),'.tif');
    SXRfilename = strcat(getenv('SXR_IMAGE_DIR'),'/',num2str(date),'/shot',num2str(shot_SXR,'%03i'),'.tif');
    % [EE_high,EE_low] = plot_SXR_at_t(grid2D,data2D,date,shot_SXR,t,show_xpoint,show_localmax,start,interval,save,SXRfilename,filter,NL);
    plot_sxr_multi(grid2D,data2D,date,shot_SXR,show_xpoint,show_localmax,start,interval,save,SXRfilename,filter,NL);
    % Brec = clc_Breconnection(grid2D,data2D);
end
% figure;plot(data2D.trange,Brec);
