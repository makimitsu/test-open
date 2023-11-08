%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dtacqショット番号を指定して
%磁気プローブデータを保存
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%各PCのパスを定義
run define_path.m

%デジタイザ別の保存(2022/11/17)
%個別に環境変数a038_path, a039_path, a040_pathを設定する必要あり
clear all
setenv('MDSPLUS_DIR','/usr/local/mdsplus');
addpath(fullfile(getenv('MDSPLUS_DIR'), 'matlab'));
pathname.rawdata=[getenv("rsOnedrive") '/mat/pcb_raw']; %保存先

dtacqlist=39;
shotlist=2178;%【input】dtacqの保存番号
% tfshotlist=zeros(size(shotlist));
tfshotlist = 2177;
date = 230816;%【input】計測日
n=numel(shotlist);%計測データ数

%RC係数読み込み

for i=1:n
    dtacq_num=dtacqlist;
    shot=shotlist(i);
    % tfshot=tfshotlist(i);%変数の場合
    tfshot=tfshotlist;%固定値の場合
    [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
    if not(exist([pathname.rawdata,'/',num2str(date)],'dir'))
        mkdir(sprintf("%s", pathname.rawdata), sprintf("%s", num2str(date)));
    end
    save(strcat(pathname.rawdata,'/',num2str(date),'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
    if tfshot>0
        [rawdata0]=getMDSdata(dtacq_num,shot,0);
        save(strcat(pathname.rawdata,'/',num2str(date),'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
    end
end
