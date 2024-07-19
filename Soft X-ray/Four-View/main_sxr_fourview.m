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
pathname.pre_processed_directory_path = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）

PCB.restart = 0; % 今だけ

%%%%実験オペレーションの取得
prompt = {'Date:','Shot number:','a039(not necessary):','doSave:','doFilter:','doNLR:'};
definput = {'','','','','',''};
if exist('date','var')
    definput{1} = num2str(date);
end
if exist('IDXlist','var')
    definput{2} = num2str(IDXlist);
end
if exist('a039','var')
    definput{3} = num2str(a039);
end
if exist('doSave','var')
    definput{4} = num2str(doSave);
end
if exist('doFilter','var')
    definput{5} = num2str(doFilter);
end
if exist('doNLR','var')
    definput{6} = num2str(doNLR);
end
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    return
end
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2)));
a039 = str2num(cell2mat(answer(3)));
doSave = logical(str2num(cell2mat(answer(4))));
doFilter = logical(str2num(cell2mat(answer(5))));
doNLR = logical(str2num(cell2mat(answer(6))));

SXR.doSave = doSave;
SXR.doFilter = doFilter;
SXR.doNLR = doNLR;

%-----------スプレッドシートからデータ抜き取り--------------------%
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);

if ~isempty(date) && ~isempty(IDXlist)% 日付とショット入力の場合
    T=searchlog(T,'date',date);
    if isnan(T.shot(1))
        T(1, :) = [];
    end
    n_data=numel(IDXlist);%計測データ
    shotlist = [T.a039(IDXlist), T.a040(IDXlist)];
    tfshotlist = [T.a039_TF(IDXlist), T.a040_TF(IDXlist)];
    EFlist=T.EF_A_(IDXlist);
    TFlist=T.TF_kV_(IDXlist);
    dtacqlist=39.*ones(n_data,1); % 39が計測データ数だけ縦に並ぶ。
    startlist = T.SXRStart(IDXlist);
    intervallist = T.SXRInterval(IDXlist);
elseif ~isempty(a039)% a039入力の場合
    T=searchlog(T,'a039',a039);
    if isnan(T.shot(1))
        T(1, :) = [];
    end
    n_data = numel(a039);
    shotlist = [T.a039, T.a040];
    tfshotlist = [T.a039_TF,T.a040_TF];
    EFlist = T.EF_A_;
    TFlist = T.TF_kV_;
    dtacqlist=39.*ones(n_data,1);
    startlist = T.SXRStart;
    intervallist = T.SXRInterval;
    date = T.date;
    IDXlist = T.shot;
end
%-------------------------------------------------%

PCB.trange=400:800;%【input】計算時間範囲
PCB.n=50; %【input】rz方向のメッシュ数

t = 470;
SXR.show_xpoint = false;
SXR.show_localmax = false;
% doSave = false;
% doFilter = true;
% doNLR = false; %do non-linear reconstruction

% NIFSの軟X線データをドライブにコピーする
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
    [PCBdata.grid2D,PCBdata.data2D] = process_PCBdata_200ch(PCB,pathname); %process_PCBdata_200ch.mに行く
    SXR.date = date;
    SXR.shot = IDXlist(i);
    SXR.SXRfilename = strcat(getenv('SXR_IMAGE_DIR'),'/',num2str(date),'/shot',num2str(SXR.shot,'%03i'),'.tif');
    plot_sxr_multi(PCBdata,SXR);
end