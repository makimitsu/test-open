%デジタイザ別の保存(2022/11/17)
%個別に環境変数a038_path, a039_path, a040_pathを設定する必要あり
clear all
setenv('MDSPLUS_DIR','/usr/local/mdsplus');
addpath(fullfile(getenv('MDSPLUS_DIR'), 'matlab'));
pathname.rawdata=[getenv("rsGdrive") '/pcb']; %保存先
dtacq_num=39;

date = 230317;%【input】計測日

%時短用(固定値で良い)
load_s = 5300;%実験ログ読み始め行番号(230309~)
load_f = 10000;%実験ログ読み終わり行番号(固定)

%実験ログがcdになければ最新の実験ログをwebから取得
if exist('exp_log.xlsx','file')
    FileInfo = dir('exp_log.xlsx');
    DateNum = FileInfo.datenum;
    FileDate = datetime(DateNum, 'ConvertFrom', 'datenum', 'Format', 'yyMMdd');
    if str2double(string(FileDate)) > date
        disp(append(string(FileDate), '更新の実験ログを使用します。'))
    else
        run save_log.m
        disp('最新の実験ログ(exp_log.xlsx)を保存しました。')
    end
else
    run save_log.m
    disp('最新の実験ログ(exp_log.xlsx)を保存しました。')
end

%実験ログ中の実験日に対応する範囲を特定
exp_log = readmatrix('exp_log.xlsx','Sheet','log','Range', ['A' num2str(load_s) ':AR' num2str(load_f)]);
[n_row,n_col] = size(exp_log);
begin_row = find(exp_log(:,3) == date);%実験日の最初のshotの行番号を取得
if isempty(begin_row)
    disp('実験日が実験ログ中に存在しません。')
    return
end
end_row = begin_row;
while end_row<n_row && isnan(exp_log(end_row+1,3)) && exp_log(end_row+1,4)%日付がNaN&&shot番号が記入済=実験日のshot
    end_row = end_row+1;
end%実験日の最後のshotの行番号を取得
for i = begin_row : end_row
    a039shot = exp_log(i,8);%a039ショット番号
    a039tfshot = exp_log(i,9);%a039TFショット番号
    %RC係数読み込み
    if dtacq_num == 39
        shot = a039shot;
        tfshot = a039tfshot;
    end
    if shot>0 && tfshot>0
        [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
        save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
        if tfshot>0
            [rawdata0]=getMDSdata(dtacq_num,shot,0);
            save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
        end
    end
end

