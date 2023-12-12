%%%%%%%%%%%%%%%%%%%%%%%%
% Top-level file for calculating and plotting SXR emission for four-view
% experimental setup
%%%%%%%%%%%%%%%%%%%%%%%%
% clear
% close all
clearvars -except date IDXlist doSave doFilter doNLR time
addpath '/Users/shinjirotakeda/Documents/GitHub/test-open/pcb_experiment'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.NIFS=getenv('NIFS_path');%192.168.1.111
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所;
pathname.pre_processed_directory = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）

%%%%実験オペレーションの取得
prompt = {'Date:','Shot number:','Time','doSave:','doFilter:','doNLR:'};
definput = {'','','','','',''};
if exist('date','var')
    definput{1} = num2str(date);
end
if exist('IDXlist','var')
    definput{2} = num2str(IDXlist);
end
if exist('time','var')
    definput{3} = num2str(time);
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
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2)));
time = str2num(cell2mat(answer(2)));
doSave = logical(str2num(cell2mat(answer(4))));
doFilter = logical(str2num(cell2mat(answer(5))));
doNLR = logical(str2num(cell2mat(answer(6))));

SXR.t = time;
SXR.doSave = doSave;
SXR.doFilter = doFilter;
SXR.doNLR = doNLR;

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
n_data=numel(IDXlist);%計測データ数
shotlist_a039 = T.a039(IDXlist);
shotlist_a040 = T.a040(IDXlist);
shotlist = [shotlist_a039, shotlist_a040];
tfshotlist_a039 = T.a039_TF(IDXlist);
tfshotlist_a040 = T.a040_TF(IDXlist);
tfshotlist = [tfshotlist_a039, tfshotlist_a040];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);
startlist = T.SXRStart(IDXlist);
intervallist = T.SXRInterval(IDXlist);

PCB.trange=400:800;%【input】計算時間範囲
PCB.n=50; %【input】rz方向のメッシュ数

SXR.show_xpoint = false;
SXR.show_localmax = false;

for i=1:n_data
    disp(strcat('(',num2str(i),'/',num2str(n_data),')'));
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
    [PCBdata.grid2D,PCBdata.data2D] = process_PCBdata_280ch(PCB,pathname);
    SXR.date = date;
    SXR.shot = IDXlist(i);
    SXR.SXRfilename = strcat(getenv('SXR_IMAGE_DIR'),'/',num2str(date),'/shot',num2str(SXR.shot,'%03i'),'.tif');
    plot_sxr_at_t(PCBdata,SXR);
end