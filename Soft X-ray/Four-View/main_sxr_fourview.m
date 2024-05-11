%%%%%%%%%%%%%%%%%%%%%%%%
% Top-level file for calculating and plotting SXR emission for four-view
% experimental setup
%%%%%%%%%%%%%%%%%%%%%%%%
% clear
% close all
clearvars -except date IDXlist doSave doFilter doNLR
addpath '/Users/shohgookazaki/Documents/GitHub/test-open/pcb_experiment'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
% ~/Documents/MATLAB にてstartup.mを作って、その中でsetenv('パス名','アドレス')していくと自動になる。
pathname.NIFS=getenv('NIFS_path');%192.168.1.111
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所;
pathname.pre_processed_directory = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）

%%%%実験オペレーションの取得
prompt = {'Date:','Shot number:','doSave:','doFilter:','doNLR:'};
definput = {'','','','',''};
if exist('date','var')
    definput{1} = num2str(date);
end
if exist('IDXlist','var')
    definput{2} = num2str(IDXlist);
end
if exist('doSave','var')
    definput{3} = num2str(doSave);
end
if exist('doFilter','var')
    definput{4} = num2str(doFilter);
end
if exist('doNLR','var')
    definput{5} = num2str(doNLR);
end
dlgtitle = 'Input';
dims = [1 35];
% if exist('date','var') && exist('IDXlist','var')
%     definput = {num2str(date),num2str(IDXlist)};
% else
%     definput = {'',''};
% end
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    return
end
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2)));
doSave = logical(str2num(cell2mat(answer(3))));
doFilter = logical(str2num(cell2mat(answer(4))));
doNLR = logical(str2num(cell2mat(answer(5))));

SXR.doSave = doSave;
SXR.doFilter = doFilter;
SXR.doNLR = doNLR;

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
n_data=numel(IDXlist);%計測データ数
% shotlist_a039 = T.a039(IDXlist);
% shotlist_a040 = T.a040(IDXlist);
% shotlist = [shotlist_a039, shotlist_a040];
shotlist = [T.a039(IDXlist), T.a040(IDXlist)];
% tfshotlist_a039 = T.a039_TF(IDXlist);
% tfshotlist_a040 = T.a040_TF(IDXlist);
% tfshotlist = [tfshotlist_a039, tfshotlist_a040];
tfshotlist = [T.a039_TF(IDXlist), T.a040_TF(IDXlist)];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1); % 39が計測データ数だけ縦に並ぶ。
startlist = T.SXRStart(IDXlist);
intervallist = T.SXRInterval(IDXlist);

% % %直接入力の場合【注意】全て同じサイズの行列になるように記入
% dtacqlist=39;
% shotlist=1118;%【input】実験ログのa039の番号
% tfshotlist=1106;%【input】実験ログのa039_TFの番号
% date = 230315;%【input】計測日
% n_data=numel(shotlist);%計測データ数
% EFlist = 150;%【input】EF電流

PCB.trange=400:800;%【input】計算時間範囲
PCB.n=50; %【input】rz方向のメッシュ数

t = 470;
SXR.show_xpoint = false;
SXR.show_localmax = false;
% doSave = false;
% doFilter = true;
% doNLR = false; %do non-linear reconstruction

copyFolderIfNotExist(strcat(getenv('SXR_IMAGE_DIR'),'/',num2str(date)), strcat(getenv('NIFS_path'),'/',num2str(date)));

for i=1:n_data
    disp(strcat('(',num2str(i),'/',num2str(n_data),')'));
    % dtacq_num=dtacqlist;
    PCB.idx = IDXlist(i);
    PCB.shot=shotlist(i,:);
    PCB.tfshot=tfshotlist(i,:);
    if PCB.shot == PCB.tfshot
        PCB.tfshot = [0,0];
    end
    PCB.i_EF=EFlist(i);
    PCB.date = date;
    TF=TFlist(i);
    SXR.start = startlist(i);
    SXR.interval = intervallist(i);
    % [PCBdata.grid2D,PCBdata.data2D] = process_PCBdata_280ch(date, shot, tfshot, pathname, n,i_EF,trange);
    % [PCBdata.grid2D,PCBdata.data2D] = process_PCBdata_280ch(PCB,pathname);
    [PCBdata.grid2D,PCBdata.data2D] = process_PCBdata_200ch(PCB,pathname);
    % grid2D = NaN;
    % data2D = NaN;
    SXR.date = date;
    SXR.shot = IDXlist(i);
    % SXR.SXRfilename = strcat(getenv('SXR_IMAGE_DIR'),'/',num2str(date),'/shot',num2str(shot_SXR,'%03i'),'.tif');
    SXR.SXRfilename = strcat(getenv('SXR_IMAGE_DIR'),'/',num2str(date),'/shot',num2str(SXR.shot,'%03i'),'.tif');
    % plot_sxr_multi(grid2D,data2D,date,shot_SXR,show_xpoint,show_localmax,start,interval,doSave,SXRfilename,doFilter,doNLR);
    plot_sxr_multi(PCBdata,SXR);
    % plot_sxr_at_t(grid2D,data2D,date,shot,t,show_xpoint,show_localmax,start,interval,doSave,SXRfilename,doFilter,NL)
end