%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dtacqショット番号を実験ログから
%自動取得して磁気プローブデータを保存
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%個別に環境変数a038_path, a039_path, a040_pathを設定する必要あり
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス
setenv('MDSPLUS_DIR','/usr/local/mdsplus');
addpath(fullfile(getenv('MDSPLUS_DIR'), 'matlab'));
%各PCのパスを定義
run define_path.m

cal_begin = 1;%計算開始shot番号

date = 230831;%【input】計測日
IDXlist= 2:71;

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
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

%dtacqdataをmat形式で保存
for i=1:n_data
    dtacq_num=dtacqlist;
    shot=shotlist(i,:);
    tfshot=tfshotlist(i,:);
    if shot == tfshot
        tfshot = [0,0];
    end
    i_EF=EFlist(i);
    TF=TFlist(i);
    filename1 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
    filename3 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot0.mat');
    if exist(filename1,"file")==0
        disp('No rawdata file of a039 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(39, shot(1), tfshot(1),rawdataPath)
        % disp(['File:',filename1,' does not exit']);
        % return
    end
    if exist(filename3,"file")==0
        disp('No rawdata0 file of a039 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(39, shot(1), 0,rawdataPath)
        % disp(['File:',filename1,' does not exit']);
        % return
    end
    filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
    filename4 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot0.mat');
    if exist(filename2,"file")==0
        disp('No rawdata file of a040 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(40, shot(2), tfshot(2),rawdataPath)
        % disp(['File:',filename2,' does not exit']);
        % return
    end
    if exist(filename4,"file")==0
        disp('No rawdata0 file of a040 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(40, shot(2), 0,rawdataPath)
        % disp(['File:',filename1,' does not exit']);
        % return
    end
end

