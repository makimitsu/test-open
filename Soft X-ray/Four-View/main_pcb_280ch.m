% clear
% close all
clearvars -except date IDXlist doCheck
addpath '/Users/shinjirotakeda/Documents/GitHub/test-open/pcb_experiment'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス

%%%%%%%%%%%%%%%%%%%%%%%%
%280ch用新規pcbプローブのみでの磁気面（Bz）
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパスx
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.pre_processed_directory = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）

%%%%実験オペレーションの取得
prompt = {'Date:','Shot number:','doCheck:'};
definput = {'','',''};
if exist('date','var')
    definput{1} = num2str(date);
end
if exist('IDXlist','var')
    definput{2} = num2str(IDXlist);
end
if exist('doCheck','var')
    definput{3} = num2str(doCheck);
end
dlgtitle = 'Input';
dims = [1 35];
% if exist('date','var') && exist('IDXlist','var') && exist('doCheck','var')
%     definput = {num2str(date),num2str(IDXlist),num2str(doCheck)};
% else
%     definput = {'','',''};
% end
% definput = {'','',''};
% definput = {num2str(date),num2str(IDXlist),num2str(doCheck)};
answer = inputdlg(prompt,dlgtitle,dims,definput);
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2)));
doCheck = logical(str2num(cell2mat(answer(3))));

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
% date=230714;
T=searchlog(T,node,date);
% IDXlist= 1; %[5:50 52:55 58:59];%[4:6 8:11 13 15:19 21:23 24:30 33:37 39:40 42:51 53:59 61:63 65:69 71:74];
n_data=numel(IDXlist);%計測データ数
shotlist_a039 =T.a039(IDXlist);
shotlist_a040 = T.a040(IDXlist);
shotlist = [shotlist_a039, shotlist_a040];
tfshotlist_a039 =T.a039_TF(IDXlist);
tfshotlist_a040 =T.a040_TF(IDXlist);
tfshotlist = [tfshotlist_a039, tfshotlist_a040];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);

trange=400:800;%【input】計算時間範囲
n=40; %【input】rz方向のメッシュ数

% doCheck = false;
% doCheck = true;
% if ~doCheck
%     figure('Position', [0 0 1500 1500],'visible','on');
% end

for i=1:n_data
    % dtacq_num=dtacqlist;
    shot=shotlist(i,:);
    tfshot=tfshotlist(i,:);
    if shot == tfshot
        tfshot = [0,0];
    end
    i_EF=EFlist(i);
    TF=TFlist(i);
    if doCheck
        check_signal(date, shot, tfshot, pathname, n);
    else
        plot_psi280ch(date, shot, tfshot, pathname,n,i_EF,trange);
    end
end