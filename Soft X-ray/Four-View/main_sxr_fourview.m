% Top-level file for 4 SXR emission plot
% このファイルは，SXRと磁気面を重ねて表示します．

clear
% パス関係
pathname.NIFS=getenv('NIFS_path');% 192.168.1.111
pathname.rawdata=getenv('rawdata_path');% dtacqのrawdataの保管場所;
pathname.pre_processed_directory = getenv('pre_processed_directory_path');% 計算結果の保存先
addpath(genpath('/Users/yuleo/Documents/GitHub/test-open'));

% --
trange=440:500;% 計算時間範囲
n=50; % rzメッシュ数
SXR_Time = 475; % 見たい時間
locate_Xpoint = false;
locate_Localmax = false;
Save = true;
ApplyFilter = false;
ApplyNonLiner = false;
% --

prompt = {'Date:','Shot number:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'',''};
answer = inputdlg(prompt,dlgtitle,dims,definput);
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2)));

% 実験オペレーション取得
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';% 実験ログ固有ID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
n_data=numel(IDXlist);% 計測データ数
shotlist=T.a039(IDXlist);
tfshotlist=T.a039_TF(IDXlist);
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);
startlist = T.SXRStart(IDXlist);
intervallist = T.SXRInterval(IDXlist);

for i=1:n_data
    dtacq_num=dtacqlist(i);
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    i_EF=EFlist(i);
    TF=TFlist(i);
    start = startlist(i);
    interval = intervallist(i);
    [grid2D,data2D] = process_PCBdata_280ch(date, shot, tfshot, pathname, n,i_EF,trange);
    Shotnumber_SXR = IDXlist(i);
    SXRFilePath = strcat(getenv('SXR_IMAGE_DIR'),'/',num2str(date),'/shot',num2str(Shotnumber_SXR,'%03i'),'.tif');
    plot_SXR_at_t(grid2D,data2D,date,Shotnumber_SXR,SXR_Time,locate_Xpoint,locate_Localmax,start,interval,Save,SXRFilePath,ApplyFilter,ApplyNonLiner);
    % plot_sxr_multi(grid2D,data2D,date,shot_SXR,show_xpoint,show_localmax,start,interval,save,SXRFilePath,filter,NL);
    % Brec = clc_Breconnection(grid2D,data2D);
end
% figure;plot(data2D.trange,Brec);
