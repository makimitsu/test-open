%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dtacqショット番号を実験ログから
%自動取得して磁気プローブデータを保存
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%個別に環境変数a038_path, a039_path, a040_pathを設定する必要あり
clear all
setenv('MDSPLUS_DIR','/usr/local/mdsplus');
addpath(fullfile(getenv('MDSPLUS_DIR'), 'matlab'));
%各PCのパスを定義
run define_path.m
dtacq_num=39;

cal_begin = 1;%計算開始shot番号

date = 230721;%【input】計測日

%実験ログ読み取り
[exp_log,index,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

if not(exist([pathname.rawdata,'/',num2str(date)],'dir'))
    mkdir(sprintf("%s", pathname.rawdata), sprintf("%s", num2str(date)));
end

%dtacqdataをmat形式で保存
for i = begin_row+(cal_begin-1) : end_row
    a039shot = exp_log(i,index.a039);%a039ショット番号
    a039tfshot = exp_log(i,index.a039_TF);%a039TFショット番号
    %RC係数読み込み
    if dtacq_num == 39
        shot = a039shot;
        tfshot = a039tfshot;
    end
    if shot>0 && tfshot>0
        [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
        if not(exist([pathname.rawdata,'/',num2str(date)],'dir'))
            mkdir(sprintf("%s", pathname.rawdata), sprintf("%s", num2str(date)));
        end
        save(strcat(pathname.rawdata,'/',num2str(date),'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
        if tfshot>0
            [rawdata0]=getMDSdata(dtacq_num,shot,0);
            save(strcat(pathname.rawdata,'/',num2str(date),'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
        end
    end
end

